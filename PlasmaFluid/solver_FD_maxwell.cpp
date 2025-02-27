#include <cmath>
#include "solver_FD_maxwell.hpp"
#include "petscsys.h"
#include "Table.hpp"
using namespace std ;
CFDMaxwell::CFDMaxwell()
{
}
void CFDMaxwell::Init(  boost::shared_ptr<CConfig> &config ,boost::shared_ptr<CVariable> &var )
{

	if ( mpi_rank == 0 ){
		cout<<"Creat FD_maxwell"<<endl;
	}
    its=0 ;

	//Warning : Power and current controlling is very tricky, the following result is kind of numerical experiment result,
	//			so there are lots of varible here. It's better to see the ppt file to understand the definition of them.

    json &ICP_simulation_condition    = *(FDMaxwell_coupled_eqs.get_json_input_parameter("ICP_simulation_condition")) ;
	var->Coil_frequency 			= ICP_simulation_condition["Coil_frequency"] ; 			// Coil frequency, (ex : 13.56 MHz)
    var->Coil_Current 				= ICP_simulation_condition["Initial_Coil_current"] ; 	// Coil current, (ex : 20A)
    var->Coil_area 					= ICP_simulation_condition["Coil_area"] ;   			// Coil crosection area (ex : 0.000001 m^2)
    var->Coil_power 				= ICP_simulation_condition["Coil_power"] ;  			// Final coil power (ex : 400 W)
    var->Coil_change_factor 		= ICP_simulation_condition["Coil_change_factor"] ; 		// Adjust coil current when it need change. Not use now.
    var->power_grows_rate 			= ICP_simulation_condition["power_grows_rate"] ; 		// Adjust coil power when it growing. (ex : 10 W/cycle)
    var->Controlled_Coil_power 		= ICP_simulation_condition["Initial_Coil_power"] ; 		// Expection of coil power.
    var->Max_current 				= ICP_simulation_condition["Max_current"] ; 			// Maximum current. (ex : 80A)
    var->omega						= 2 * var->PI * var->Coil_frequency ;
	var->power_inductive			= var->Controlled_Coil_power ; 							// EM_power, which integral EM_power absorption over all computation domain.
	
	
	// Give the permeability for all cells.
    for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++){
    	Cell *cell		=	FDMaxwell_Re.get_cell( cth ) ;
		
		if (	cell->type == MPP_cell_tag[ "DIELECTRIC_FVFD" ]	)
		{ 
				var->eps_FVFD	[ cth ]	=   ICP_simulation_condition["DIELECTRIC_eps_r"] ;	 		
		} else 
		{
				var->eps_FVFD	[ cth ]	=	1	; 				
		}
    }  
    UltraMPPComputeCurrentDenAndSourceTerm( config, var ) ;

	// Read collision table
	FVFD_CollTable.Init( config->CasePath+"5Collision.inp" ) ;
	
	FDMaxwell_Re.apply_linear_solver_setting();
    FDMaxwell_Im.apply_linear_solver_setting();
    FDMaxwell_coupled_eqs.apply_linear_solver_setting();
}
void CFDMaxwell::SOLVE( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{	
	// Count that how many times Maxwell equ. solved. This will output to ICP.dat. 
	var->Maxwell_solver_count++ ;
	
	/*--- Matrix A ---*/
    FDMaxwell_Re.before_matrix_construction() ;
    FDMaxwell_Im.before_matrix_construction() ;
    FDMaxwell_coupled_eqs.before_matrix_construction() ;
//--
    FDMaxwell_Re.add_laplacian_matrix_form_op();
    FDMaxwell_Im.add_laplacian_matrix_form_op();
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

	// Compute plasma conductivity. (Both real & imaginary part)
    for( int cth = 0; cth < plasma.Mesh.cell_number; cth++){
		
    	var->collision_frequency = var->TotalNumberDensity[ cth ]* FVFD_CollTable.GetValue( var->T[ 0 ][ cth ] ) ; //unit
        var->sigma_p_Re_plasma[ cth ] 	=   unit_charge*unit_charge*var->U0[0][ cth ]*var->collision_frequency	/electron_mass/( var->collision_frequency * var->collision_frequency + var->omega * var->omega ) ;
 		var->sigma_p_Im_plasma[ cth ] 	= - unit_charge*unit_charge*var->U0[0][ cth ]*var->omega				/electron_mass/( var->collision_frequency * var->collision_frequency + var->omega * var->omega ) ; 
		
    } 
    
	// Map plasma conductivity from plasma mesh to FVFD mesh.
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Re_FVFD"]		, var->VarTag["sigma_p_Re_plasma"] ) ;     
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Im_FVFD"]		, var->VarTag["sigma_p_Im_plasma"] ) ;  
 
  	    
    for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++){

        Cell *cell = FDMaxwell_Re.get_cell(cth);
        
		// Compute effective, square of effective Helmhotlz constant. (Both real & imaginary part)
        var->k_square_Re[ cth ] 	=	var->omega * vacuum_permeability * var->sigma_p_Im_FVFD[ cth ] + var->omega*var->omega*var->eps_FVFD[ cth ]/vacuum_light_speed/vacuum_light_speed ;
        var->k_square_Im[ cth ] 	= -	var->omega * vacuum_permeability * var->sigma_p_Re_FVFD[ cth ] ;

		// Building matrix for diagonal term.
        FDMaxwell_Re.add_entry_in_matrix(cth, cell->id, -1/cell->r[0]/cell->r[0] + var->k_square_Re[ cth ]);
        FDMaxwell_Im.add_entry_in_matrix(cth, cell->id, -1/cell->r[0]/cell->r[0] + var->k_square_Re[ cth ]);

        //coupled_eqs.add_entry_in_matrix( equation number ,  variable_tag , cth , cell->mat_id, entry_value ) ;
		FDMaxwell_coupled_eqs.add_entry_in_matrix(0,  var->E_phi_Im.tag_current, cth , cell->id, -var->k_square_Im[ cth ] ) ;
		FDMaxwell_coupled_eqs.add_entry_in_matrix(1,  var->E_phi_Re.tag_current, cth , cell->id,  var->k_square_Im[ cth ] ) ;

    }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    FDMaxwell_coupled_eqs.add_matrix( &FDMaxwell_Re );
    FDMaxwell_coupled_eqs.add_matrix( &FDMaxwell_Im );

//--
    FDMaxwell_Re.finish_matrix_construction() ;
    FDMaxwell_Im.finish_matrix_construction() ;
    FDMaxwell_coupled_eqs.finish_matrix_construction() ;

	/*--- Source B ---*/
    FDMaxwell_Re.set_bc_value(MPP_face_tag["GROUND_FVFD"], 0.0,var->E_phi_Re.face );
    FDMaxwell_Im.set_bc_value(MPP_face_tag["GROUND_FVFD"], 0.0,var->E_phi_Im.face );

	// Compute source term, omega * mu0 * J.
    UltraMPPComputeCurrentDenAndSourceTerm( config, var ) ;

    FDMaxwell_Re.before_source_term_construction();
    FDMaxwell_Im.before_source_term_construction();
    FDMaxwell_coupled_eqs.before_source_term_construction();

    FDMaxwell_Re.add_laplacian_source_term_op(var->Re_eq_source, var->E_phi_Re.face);
    FDMaxwell_Im.add_laplacian_source_term_op(var->Im_eq_source, var->E_phi_Im.face);

    FDMaxwell_Re.finish_source_term_construction( );
    FDMaxwell_Im.finish_source_term_construction( );
    FDMaxwell_coupled_eqs.finish_source_term_construction( );

    FDMaxwell_coupled_eqs.get_solution();
    its = FDMaxwell_coupled_eqs.get_iteration_number() ;
	
    UltraMPPComputePowerAbsorptionFromMaxwell( config, var ) ;
	//UltraMPPComputeInstantPowerAbsorptionFromMaxwell( config, var ) ;

}
void CFDMaxwell::UltraMPPComputeCurrentDenAndSourceTerm( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	double current_temp =0.0 , current_diff=0.0;
	/*Only the 'coil' region. */
  for ( int cth=0 ; cth< FDMaxwell_Re.Mesh.cell_number ; cth++ ) {

    Cell *cell = FDMaxwell_Re.get_cell( cth ) ;

    if ( cell->type == MPP_cell_tag[ "coil" ] ) {
    	// J = I/Area
		var->CurrentDen[ cth ] = var->Coil_Current/var->Coil_area ;
	}else{
		var->CurrentDen[ cth ] = 0 ;		
	}

	// source term
	var->Re_eq_source[ cth ] = 0.0 ; 	
	var->Im_eq_source[ cth ] = var->omega * vacuum_permeability * var->CurrentDen[ cth ] ; 

  }//cell loop.
  
  	// This section is for adjusting coil current, please see the algorithm
 	if( var->power_inductive < var->Controlled_Coil_power){	
    	current_temp = var->Coil_Current * var->Controlled_Coil_power / var->power_inductive ;
    	current_diff =  current_temp - var->Coil_Current ;
    	var->Coil_Current = var->Coil_Current + 25/var->Controlled_Coil_power * current_diff;
		 	
 //   	var->Coil_Current = var->Coil_Current * var->Coil_change_factor; // Not use now.
    }else if(var->power_inductive > var->Controlled_Coil_power){
    	current_temp = var->Coil_Current * var->Controlled_Coil_power / var->power_inductive ;
    	current_diff =  current_temp - var->Coil_Current ;
    	var->Coil_Current = var->Coil_Current + 25/var->Controlled_Coil_power * current_diff;
	}
	
	// Limiting the coil current in case it gets too high.
	if (var->Coil_Current > var->Max_current) var->Coil_Current = var->Max_current ; 

  plasma.syn_parallel_cell_data( var->VarTag["CurrentDen"] );   
}
void CFDMaxwell::UltraMPPComputePowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	for( int cth = 0; cth < plasma.Mesh.cell_number; cth++){
	
    	var->collision_frequency = var->TotalNumberDensity[ cth ]* FVFD_CollTable.GetValue( var->T[ 0 ][ cth ] ) ; //unit
        var->sigma_p_Re_plasma[ cth ] 	=   unit_charge*unit_charge*var->U0[0][ cth ]*var->collision_frequency	/electron_mass/( var->collision_frequency * var->collision_frequency + var->omega * var->omega ) ;
 		var->sigma_p_Im_plasma[ cth ] 	= - unit_charge*unit_charge*var->U0[0][ cth ]*var->omega				/electron_mass/( var->collision_frequency * var->collision_frequency + var->omega * var->omega ) ; 
		
    } 

	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Re_FVFD"]		, var->VarTag["sigma_p_Re_plasma"] ) ;     
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Im_FVFD"]		, var->VarTag["sigma_p_Im_plasma"] ) ;  
	
	for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++){		
		var->Power_Absorption_FVFD[ cth ] = 0.5 * ( (var->sigma_p_Re_FVFD[ cth ] * var->E_phi_Re.current[ cth ] + var->sigma_p_Im_FVFD[ cth ] * var->E_phi_Im.current[ cth ] ) * var->E_phi_Re.current[ cth ]    
								             + (var->sigma_p_Re_FVFD[ cth ] * var->E_phi_Im.current[ cth ] - var->sigma_p_Im_FVFD[ cth ] * var->E_phi_Re.current[ cth ] ) * var->E_phi_Im.current[ cth ] ) ;								                								            
	} 
	
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["Power_Absorption_plasma"]	, var->VarTag["Power_Absorption_FVFD"] 	   ) ; 
}
void CFDMaxwell::UltraMPPComputeInstantPowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++)
	{
		var->Power_Absorption_FVFD[ cth ] =  sin(var->omega*var->PhysicalTime) * cos(var->omega*var->PhysicalTime)*
				                ( - 2.0 * var->sigma_p_Re_FVFD[ cth ] * var->E_phi_Re.current[ cth ] * var->E_phi_Im.current[ cth ]
							            - var->sigma_p_Im_FVFD[ cth ] * var->E_phi_Im.current[ cth ] * var->E_phi_Im.current[ cth ]
										+ var->sigma_p_Im_FVFD[ cth ] * var->E_phi_Re.current[ cth ] * var->E_phi_Re.current[ cth ] );
										
		var->Power_Absorption_FVFD[ cth ] +=   var->sigma_p_Re_FVFD[ cth ] * var->E_phi_Re.current[ cth ] * var->E_phi_Re.current[ cth ] * cos(var->omega*var->PhysicalTime) * cos(var->omega*var->PhysicalTime);
		var->Power_Absorption_FVFD[ cth ] +=   var->sigma_p_Im_FVFD[ cth ] * var->E_phi_Re.current[ cth ] * var->E_phi_Im.current[ cth ] * cos(var->omega*var->PhysicalTime) * cos(var->omega*var->PhysicalTime);
		var->Power_Absorption_FVFD[ cth ] +=   var->sigma_p_Re_FVFD[ cth ] * var->E_phi_Im.current[ cth ] * var->E_phi_Im.current[ cth ] * sin(var->omega*var->PhysicalTime) * sin(var->omega*var->PhysicalTime);
		var->Power_Absorption_FVFD[ cth ] += - var->sigma_p_Im_FVFD[ cth ] * var->E_phi_Re.current[ cth ] * var->E_phi_Im.current[ cth ] * sin(var->omega*var->PhysicalTime) * sin(var->omega*var->PhysicalTime);
	}
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["Power_Absorption_plasma"]	, var->VarTag["Power_Absorption_FVFD"] 	   ) ; 	
}
void CFDMaxwell::UltraMPPComputeTotalPower( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	//Calculate power [W]
	var->power_inductive = 0 ;
	var->power_static = 0 ;	
	for ( int cth = 0 ; cth < plasma.Mesh.cell_number ; cth++ )
	{
		Cell *cell	=	plasma.get_cell( cth ) ;
		var->power_inductive	+=	var->Power_Absorption_plasma[ cth ] * cell->volume;
		var->power_static	+=	var->JouleHeating[0][ cth ] * cell->volume ;
	}
	var->power_inductive =	 plasma.parallel_sum( &var->power_inductive ) ;
	var->power_static	 =	 plasma.parallel_sum( &var->power_static ) ;
}


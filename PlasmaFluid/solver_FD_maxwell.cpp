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
		
	FDMaxwell_Re.set_linear_solver_library("PETSC");
    FDMaxwell_Im.set_linear_solver_library("PETSC");
    FDMaxwell_coupled_eqs.set_linear_solver_library("PETSC");
    
    FDMaxwell_Re.set_N_variable_number( 1 ) ;
    FDMaxwell_Im.set_N_variable_number( 1 ) ;
    FDMaxwell_coupled_eqs.set_N_variable_number( 2 ) ;
    
    FDMaxwell_Re.load_mesh( "FD_maxwell_mesh" ) ;
    FDMaxwell_Im.load_mesh( "FD_maxwell_mesh" ) ;
    FDMaxwell_coupled_eqs.load_mesh( "FD_maxwell_mesh" ) ;
    
    json &ICP_simlation_condition = *(FDMaxwell_coupled_eqs.get_json_input_parameter("ICP_simlation_condition") );
    
    var->Coil_frequency = ICP_simlation_condition["Coil_frequency"] ;
    var->Coil_Current 	= ICP_simlation_condition["Coil_Current"] ;
    var->omega			= 2 * var->PI * var->Coil_frequency ;
    
    for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++){
    	Cell *cell		=	FDMaxwell_Re.get_cell( cth ) ;
		
		if (	cell->type == MPP_cell_tag[ "Quartz" ]	)
		{ 
				var->eps_FVFD	[ cth ]	=   ICP_simlation_condition["Quartz_eps_r"]	;				
		} else 
		{
				var->eps_FVFD	[ cth ]	=	1	; 				
		}
    }   
    
    UltraMPPComputeCurrentDanAndSourceTerm( config, var ) ;
}
void CFDMaxwell::SOLVE( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{	

	double collision_frequency = 0.0 ;
	/*--- Matrix A ---*/
    FDMaxwell_Re.before_matrix_construction() ;
    FDMaxwell_Im.before_matrix_construction() ;
    FDMaxwell_coupled_eqs.before_matrix_construction() ;
//--
    FDMaxwell_Re.add_laplacian_matrix_form_op();
    FDMaxwell_Im.add_laplacian_matrix_form_op();
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

    for( int cth = 0; cth < plasma.Mesh.cell_number; cth++){
    	collision_frequency = var->TotalNumberDensity[ cth ]*CollTable.GetValue( var->T[ 0 ][ cth ] ) ; //unit
        
        var->sigma_p_Re_plasma[ cth ] 	=   unit_charge*unit_charge*var->U0[0][ cth ]*var->collision_frequency	/electron_mass/( var->collision_frequency * var->collision_frequency + var->omega * var->omega ) ;
 		var->sigma_p_Im_plasma[ cth ] 	= - unit_charge*unit_charge*var->U0[0][ cth ]*var->omega				/electron_mass/( var->collision_frequency * var->collision_frequency + var->omega * var->omega ) ; 
    }
    
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Re_FVFD"]		, var->VarTag["sigma_p_Re_plasma"] ) ;     
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Im_FVFD"]		, var->VarTag["sigma_p_Im_plasma"] ) ;    
		    
    for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++){

        Cell *cell = FDMaxwell_Re.get_cell(cth);
        
      
        var->k_square_Re[ cth ] 	=	var->omega * vacuum_permeability * var->sigma_p_Im_FVFD[ cth ] + var->omega*var->omega*var->eps_FVFD[ cth ]/vacuum_light_speed/vacuum_light_speed ;
        var->k_square_Im[ cth ] 	= -	var->omega * vacuum_permeability * var->sigma_p_Re_FVFD[ cth ] ;
				
        FDMaxwell_Re.add_entry_in_matrix(cth, cell->mat_id, -1/cell->r[0]/cell->r[0] + var->k_square_Re[ cth ]);
        FDMaxwell_Im.add_entry_in_matrix(cth, cell->mat_id, -1/cell->r[0]/cell->r[0] + var->k_square_Re[ cth ]);
        
        //coupled_eqs.add_entry_in_matrix( equation number ,  variable_tag , cth , cell->mat_id, entry_value ) ;
		FDMaxwell_coupled_eqs.add_entry_in_matrix(0,  var->E_phi_Im.tag_current, cth , cell->mat_id, -var->k_square_Im[ cth ] ) ;
		FDMaxwell_coupled_eqs.add_entry_in_matrix(1,  var->E_phi_Re.tag_current, cth , cell->mat_id,  var->k_square_Im[ cth ] ) ;

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
    FDMaxwell_Re.set_bc_value(MPP_face_tag["GROUND"], 0.0,var->E_phi_Re.face );
    FDMaxwell_Im.set_bc_value(MPP_face_tag["GROUND"], 0.0,var->E_phi_Im.face );

    FDMaxwell_Re.before_source_term_construction();
    FDMaxwell_Im.before_source_term_construction();
    FDMaxwell_coupled_eqs.before_source_term_construction();

    FDMaxwell_Re.add_laplacian_source_term_op(var->Re_eq_source, var->E_phi_Re.face);
    FDMaxwell_Im.add_laplacian_source_term_op(var->Im_eq_source, var->E_phi_Im.face);

    FDMaxwell_Re.finish_source_term_construction( );
    FDMaxwell_Im.finish_source_term_construction( );
    FDMaxwell_coupled_eqs.finish_source_term_construction( );

    FDMaxwell_coupled_eqs.get_solution();
    
    UltraMPPComputePowerAbsorptionFromMaxwell( config, var ) ;
	//UltraMPPComputeInstantPowerAbsorptionFromMaxwell( config, var ) ;
}
void CFDMaxwell::UltraMPPComputeCurrentDanAndSourceTerm( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*Only the 'coil' region. */
  for ( int cth=0 ; cth< FDMaxwell_Re.Mesh.cell_number ; cth++ ) {

    Cell *cell = plasma.get_cell( cth ) ;

    if ( cell->type == MPP_cell_tag[ "COIL" ] ) {
		var->CurrentDen[ cth ] = var->Coil_Current/cell->volume ;
	}else{
		var->CurrentDen[ cth ] = 0 ;		
	}//if charged species.
	
	var->Re_eq_source[ cth ] = 0.0 ;
	var->Im_eq_source[ cth ]	= var->omega * vacuum_permeability * var->CurrentDen[ cth ] ;

  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["CurrentDen"] );
}
void CFDMaxwell::UltraMPPComputePowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
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


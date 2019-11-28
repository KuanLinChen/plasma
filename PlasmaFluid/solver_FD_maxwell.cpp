#include <cmath>
#include "solver_poisson.hpp"
#include "petscsys.h"
using namespace std ;
CFDMaxwell::CFDMaxwell()
{
}
void CFDMaxwell::Init( boost::shared_ptr<CConfig> &config )
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
    
    Coil_frequency 	= ICP_simlation_condition["Coil_frequency"] ;
    Coil_Current 	= ICP_simlation_condition["Coil_Current"] ;
    omega			= 2 * PI * Coil_frequency ;
    
    for( int cth = 0; cth < Plasma.Mesh.cell_number; cth++){
    	Cell *cell		=	PDE_solver.get_cell( cth ) ;
		cpiID	[ cth ]	=	cell->mpi_id ;
		
		if (	cell->type == CellType[ "Quartz" ]	)
		{ 
				eps_FVFD	[ cth ]	=   ICP_simlation_condition["Quartz_eps_r"]	;				
		} else 
		{
				eps_FVFD	[ cth ]	=	1	; 				
		}
    }   
    
    UltraMPPComputeCurrentDan_And_SourceTerm( config, var ) ;
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

    for( int cth = 0; cth < Plasma.Mesh.cell_number; cth++){
    	collision_frequency = [cth]*CollTable.GetValue( T[ 0 ][ cth ] ) ; //unit
        
        sigma_p_Re_plasma[ cth ] 	=   unit_charge*unit_charge*variable->U0[0][ cth ]*collision_frequency	/electron_mass/( collision_frequency * collision_frequency + omega * omega ) ;
 		sigma_p_Im_plasma[ cth ] 	= - unit_charge*unit_charge*variable->U0[0][ cth ]*omega				/electron_mass/( collision_frequency * collision_frequency + omega * omega ) ; 
    }
    
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Re_FVFD"]		, var->VarTag["sigma_p_Re_plasma"] ) ;     
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["sigma_p_Im_FVFD"]		, var->VarTag["sigma_p_Im_plasma"] ) ;    
		    
    for( int cth = 0; cth < FDMaxwell_Re.Mesh.cell_number; cth++){

        Cell *cell = FDMaxwell_Re.get_cell(cth);
        
      
        k_square_Re[ cth ] 	=	omega * vacuum_permeability * sigma_p_Im_FVFD[ cth ] + omega*omega*eps_FVFD[ cth ]/vacuum_light_speed/vacuum_light_speed ;
        k_square_Im[ cth ] 	= -	omega * vacuum_permeability * sigma_p_Re_FVFD[ cth ] ;
				
        FDMaxwell_Re.add_entry_in_matrix(cth, cell->mat_id, -1/cell->r[0]/cell->r[0] + k_square_Re[ cth ]);
        FDMaxwell_Im.add_entry_in_matrix(cth, cell->mat_id, -1/cell->r[0]/cell->r[0] + k_square_Re[ cth ]);
        
        //coupled_eqs.add_entry_in_matrix( equation number ,  variable_tag , cth , cell->mat_id, entry_value ) ;
		FDMaxwell_coupled_eqs.add_entry_in_matrix(0,  E_phi_Im.tag_current, cth , cell->mat_id, -k_square_Im[ cth ] ) ;
		FDMaxwell_coupled_eqs.add_entry_in_matrix(1,  E_phi_Re.tag_current, cth , cell->mat_id,  k_square_Im[ cth ] ) ;

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
    FDMaxwell_Re.set_bc_value(MPP_face_tag["GROUND"], 0.0,E_phi_R.face );
    FDMaxwell_Im.set_bc_value(MPP_face_tag["GROUND"], 0.0,E_phi_I.face );

    FDMaxwell_Re.before_source_term_construction();
    FDMaxwell_Im.before_source_term_construction();
    FDMaxwell_coupled_eqs.before_source_term_construction();

    FDMaxwell_Re.add_laplacian_source_term_op(Re_eq_source, E_phi_Re.face);
    FDMaxwell_Im.add_laplacian_source_term_op(Im_eq_source, E_phi_Im.face);

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
  for ( int cth=0 ; cth<plasma.Mesh.cell_number ; cth++ ) {

    Cell *cell = plasma.get_cell( cth ) ;

    if ( cell_type[ cell->type ] == COIL ) {
		var->CurrentDen[ cth ] = Coil_Current/cell->volume ;
	}else{
		var->CurrentDen[ cth ] = 0 ;		
	}//if charged species.
	
	Re_eq_source[ cth ] = 0.0 ;
	Im_eq_source[ cth ]	= omega * mu_0 * CurrentDen[ cth ] ;

  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["CurrentDen"] );
}
void CFDMaxwell::UltraMPPComputePowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	for( int cth = 0; cth < eq1.Mesh.cell_number; cth++){		
		Power_Absorption_FVFD[ cth ] = 0.5 * ( (sigma_p_Re_FVFD[ cth ] * E_phi_Re.current[ cth ] + sigma_p_Im_FVFD[ cth ] * E_phi_Im.current[ cth ] ) * E_phi_Re.current[ cth ]    
								             + (sigma_p_Re_FVFD[ cth ] * E_phi_Im.current[ cth ] - sigma_p_Im_FVFD[ cth ] * E_phi_Re.current[ cth ] ) * E_phi_Im.current[ cth ] ) ;								                								            

	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["Power_Absorption_plasma"]	, var->VarTag["Power_Absorption_FVFD"] 	   ) ;  
}
void CFDMaxwell::UltraMPPComputeInstantPowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	for( int cth = 0; cth < eq1.Mesh.cell_number; cth++){
		Power_Absorption_FVFD[ cth ] =  sin(omega*PhysicalTime) * cos(omega*PhysicalTime)*
				                ( - 2.0 * sigma_p_Re_FVFD[ cth ] * E_phi_Re.current[ cth ] * E_phi_Im.current[ cth ]
							            - sigma_p_Im_FVFD[ cth ] * E_phi_Im.current[ cth ] * E_phi_Im.current[ cth ]
										+ sigma_p_Im_FVFD[ cth ] * E_phi_Re.current[ cth ] * E_phi_Re.current[ cth ] );
										
		Power_Absorption_FVFD[ cth ] +=   sigma_p_Re_FVFD[ cth ] * E_phi_Re.current[ cth ] * E_phi_Re.current[ cth ] * cos(omega*PhysicalTime) * cos(omega*PhysicalTime);
		Power_Absorption_FVFD[ cth ] +=   sigma_p_Im_FVFD[ cth ] * E_phi_Re.current[ cth ] * E_phi_Im.current[ cth ] * cos(omega*PhysicalTime) * cos(omega*PhysicalTime);
		Power_Absorption_FVFD[ cth ] +=   sigma_p_Re_FVFD[ cth ] * E_phi_Im.current[ cth ] * E_phi_Im.current[ cth ] * sin(omega*PhysicalTime) * sin(omega*PhysicalTime);
		Power_Absorption_FVFD[ cth ] += - sigma_p_Im_FVFD[ cth ] * E_phi_Re.current[ cth ] * E_phi_Im.current[ cth ] * sin(omega*PhysicalTime) * sin(omega*PhysicalTime);
	}
	FDMaxwell_Re.syn_parallel_cell_data( var->VarTag["Power_Absorption_plasma"]	, var->VarTag["Power_Absorption_FVFD"] 	   ) ;  
}


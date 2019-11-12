
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
		cout<<"Correction: "<<config->Equation[ POISSON ].Correction<<endl;
	}
	its=0 ;
		
	FDMaxwell_Re.set_linear_solver_library("PETSC");
    FDMaxwell_Im.set_linear_solver_library("PETSC");
    FDMaxwell_coupled_eqs.set_linear_solver_library("PETSC");
    
    FDMaxwell_Re.set_N_variable_number( 1 ) ;
    FDMaxwell_Im.set_N_variable_number( 1 ) ;
    FDMaxwell_coupled_eqs.set_N_variable_number( 2 ) ;
    
    FDMaxwell_Re.load_mesh( "FD_maxwell_mesh" ) ;
    FDMaxwell_Im.load_mesh( "FD_maxwell_mesh" ) ;
    FDMaxwell_coupled_eqs.load_mesh( "FD_maxwell_mesh" ) ;
    
    UltraMPPComputeCurrentDan_And_SourceTerm( config, var ) ;
}
void CFDMaxwell::SOLVE( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{	

	/*--- Matrix A ---*/
    FDMaxwell_Re.before_matrix_construction() ;
    FDMaxwell_Im.before_matrix_construction() ;
    FDMaxwell_coupled_eqs.before_matrix_construction() ;
//--
    FDMaxwell_Re.add_laplacian_matrix_form_op();
    FDMaxwell_Im.add_laplacian_matrix_form_op();
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    for( int cth = 0; cth < eq1.Mesh.cell_number; cth++){

        Cell *cell = FDMaxwell_Re.get_cell(cth);
		
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

}
void CFDMaxwell::UltraMPPComputeCurrentDan_And_SourceTerm( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
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
  plasma.syn_parallel_cell_data( var->VarTag["ChargeDen"] );
}


#include "solver_photoionization.hpp"

using namespace std;

CHelmholtz::CHelmholtz()
{
}
void CHelmholtz::Init( boost::shared_ptr<CConfig> &config )
{

	numFitTerm = 3 ;
	A      = new double [ numFitTerm ] ;
	lambda = new double [ numFitTerm ] ;

	if ( numFitTerm == 3 ) {

		/* Three exponential fit coefficients. */
		A[0] = 1.9860 ; // unit # m^-2*Torr^-2
		A[1] = 4886.0 ; // unit # m^-2*Torr^-2
		A[2] = 51.000 ; // unit # m^-2*Torr^-2

		lambda[0] = 5.530 ; // unit # m^-1*Torr^-1
		lambda[1] = 89.00 ; // unit # m^-1*Torr^-1
		lambda[2] = 14.60 ; // unit # m^-1*Torr^-1

	} else if ( numFitTerm == 6 )	{

	}

	O2PartialPressure = 150.0 ;


	/* Initial ultraMPP object. */
	helmholtz.set_linear_solver_library("PETSC");
	helmholtz.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;

	string input_json = "helmholtz3term.json" ;
	helmholtz.load_mesh( config->CasePath + input_json ) ;
	helmholtz.apply_linear_solver_setting();

	/* Extract the boundary & cell setting from input json file. */
	json &json_bc_setting    = *(helmholtz.get_json_input_parameter("boundary_setting")) ;
	json &json_cell_setting  = *(helmholtz.get_json_input_parameter("volume_setting"  )) ;

	/* Face */
	for ( int ibc = 0; ibc < json_bc_setting["name"].size(); ibc++){
		helmholtz_face_parameter[ json_bc_setting["name"][ ibc ] ] = json_bc_setting["values"][ ibc ] ;
		helmholtz_face_tag      [ json_bc_setting["name"][ ibc ] ] = helmholtz.set_bc_mapping(  json_bc_setting["name"][ ibc ],  json_bc_setting["boundary_type"][ ibc ]  );
	}

	/* Cell */
	for ( int icc = 0; icc < json_cell_setting["name"].size(); icc++ ) {
		helmholtz_cell_parameter[ json_cell_setting["name"][ icc ] ] = json_cell_setting["permittivity"][ icc ] ;
		helmholtz_cell_tag      [ json_cell_setting["name"][ icc ] ] = helmholtz.get_cell_type_mapping( json_cell_setting["name"][ icc ] ) ;
	}


	/* Allocate the memory. */
	helmholtz.set_parallel_variable( &photoionization, "photoionization" ) ;
	Tag_RHS    = helmholtz.set_parallel_cell_data( &RHS, "RHS" ) ;
	Tag_SphSum = helmholtz.set_parallel_cell_data( &SphSum, "SphSum" ) ;


	// set bc value
	for (int i=0 ; i < helmholtz.Mesh.face_number ; i++ ) {
		photoionization.face[ i ] = 0.0 ;
	}

}
void CHelmholtz::SOLVE( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	Cell *Cell_i ;
	double I_function=0.0 ;
	/* Reset zero for photoionization source term */
	for (int i = 0 ; i < helmholtz.Mesh.cell_number ; i++ ) SphSum[ i ] = 0.0 ;

	for ( int nTerm = 0 ; nTerm < numFitTerm ; nTerm++) {

		/*--- Matrix A ---*/
		helmholtz.before_matrix_construction() ;
		helmholtz.add_laplacian_matrix_form_op() ;

		for ( int i = 0 ; i < helmholtz.Mesh.cell_number ; i++ ) {
			Cell_i = helmholtz.get_cell(i) ;
    	helmholtz.add_entry_in_matrix( i, Cell_i->id, - pow( lambda[nTerm]*O2PartialPressure, 2.0) ) ;
		}
		helmholtz.finish_matrix_construction() ;

		/* Source B */
		for ( int i = 0 ; i < helmholtz.Mesh.cell_number ; i++ ) {
			Cell_i = helmholtz.get_cell(i) ;
			I_function = 3.4884e28*exp(-pow((Cell_i->r[1]-0.005)/0.0001,2.0)-pow(Cell_i->r[0]/0.0001,2.0)) ; // unit # m^-3*s^-1
			RHS[ i ]   = -A[ nTerm ] * pow( O2PartialPressure, 2.0 ) * I_function ;
		}
		helmholtz.before_source_term_construction() ;
		helmholtz.add_laplacian_source_term_op( RHS, photoionization.face ) ;
		helmholtz.finish_source_term_construction() ;
		helmholtz.get_solution( photoionization.current ) ;

		/* add photoionization source term to SphSum array.*/
		for ( int i = 0 ; i < helmholtz.Mesh.cell_number ; i++ ) 
		{
			SphSum[ i ] += photoionization.current[i] ;
		}

	} // End nTerm


	/*outout (will be delete in the future. */
	 helmholtz.set_output("helmholtz3term_test") ;
	 helmholtz.set_output( Tag_RHS    ) ;
	 helmholtz.set_output( Tag_SphSum ) ;
	 helmholtz.write_output("1") ;
}
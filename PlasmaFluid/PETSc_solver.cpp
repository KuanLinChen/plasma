#include "PETSc_solver.h"
#include <petscksp.h>
#include <petscsnes.h>
#include <petsc.h>
#include <petscsys.h>
#include "main.h"
#include <iostream>
#include <boost/assert.hpp>

using namespace std;

bool PETScSolver::flagPetscInitialize = false ;

PETScSolver::PETScSolver ( )
{
	local_unknown_number = -1 ;
	flag_local_number_initialization = false ;
	comm = MPI_COMM_WORLD ;
	PETSC_COMM_WORLD = comm ;
}

PETScSolver::~PETScSolver ( )
{

}

void PETScSolver::set_local_unknown_number ( int _u )
{
	local_unknown_number = _u ;
	flag_local_number_initialization = true ;
}

void PETScSolver::set_nz(int _MaxNZ, int *d_nnz, int *o_nnz){

	dia_nz = _MaxNZ;
	off_nz = _MaxNZ;

	dia_iNum_nz = new int[local_unknown_number];
	off_iNum_nz = new int[local_unknown_number];

	for(int i = 0; i < local_unknown_number; i++ )dia_iNum_nz[i] = d_nnz[i];
	for(int i = 0; i < local_unknown_number; i++ )off_iNum_nz[i] = o_nnz[i];

};

void PETScSolver::init()
{
	init( 1.e-8 , PETSC_PCASM ) ;
}

void PETScSolver::close()
{
	VecRestoreArray( X, &x ) ;
	//VecRestoreArray( B, &b ) ;
	VecDestroy( &X );
	VecDestroy( &B );
	MatDestroy( &A ) ;
	KSPDestroy( &ksp ) ;

}

void PETScSolver::init( int PC_TYPE )
{

	init( 1.e-8 , PC_TYPE ) ;
}

void PETScSolver::init( double tol, int PC_TYPE )
{
	int	i ;
	PetscErrorCode		ierr;
	int				size;


	BOOST_ASSERT_MSG ( local_unknown_number > 0, "Uninitialized variable local_unknown_number." ) ;

	MPI_Comm_rank( PETSC_COMM_WORLD, &mpi_id) ;
	MPI_Comm_size( PETSC_COMM_WORLD, &size ) ;

	asm_overlap = 1 ;
	total_it = 0;

	VecCreate( PETSC_COMM_WORLD, &X ) ;

	VecSetSizes( X, local_unknown_number, PETSC_DECIDE ) ;
	VecSetFromOptions( X ) ;
	VecDuplicate ( X, &B ) ;

	VecGetSize ( X, &global_unknown_number) ;
	VecGetOwnershipRange( X, &low, &high) ;
/*
MatCreateAIJ(MPI_Comm comm,PetscInt m,PetscInt n,PetscInt M,PetscInt N,PetscInt d_nz,const PetscInt d_nnz[],PetscInt o_nz,const PetscInt o_nnz[],Mat *A)
Input Parameters
comm	- MPI communicator
m	- number of local rows (or PETSC_DECIDE to have calculated if M is given) This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
n	- This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax. (or PETSC_DECIDE to have calculated if N is given) For square matrices n is almost always m.
M	- number of global rows (or PETSC_DETERMINE to have calculated if m is given)
N	- number of global columns (or PETSC_DETERMINE to have calculated if n is given)
d_nz	- number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
d_nnz	- array containing the number of nonzeros in the various rows of the DIAGONAL portion of the local submatrix (possibly different for each row) or NULL, if d_nz is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e 'm'.
o_nz	- number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows).
o_nnz	- array containing the number of nonzeros in the various rows of the OFF-DIAGONAL portion of the local submatrix (possibly different for each row) or NULL, if o_nz is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e 'm'.
*/

	MatCreateAIJ( PETSC_COMM_WORLD, local_unknown_number, local_unknown_number, global_unknown_number, global_unknown_number, dia_nz, dia_iNum_nz, off_nz, off_iNum_nz, &A )  ;

	KSPCreate( PETSC_COMM_WORLD, &ksp ) ;
	// PETSc 3.4.4
	KSPSetOperators( ksp, A, A ) ;
	// Set PC
	KSPGetPC( ksp , &pc );
	if ( size > 1 )
	{
		if ( PC_TYPE == PETSC_PCASM )
		{
			PCSetType ( pc , PCASM ) ;
			PCASMSetOverlap( pc, asm_overlap );
			PCASMSetType ( pc, PC_ASM_BASIC ) ;

		}
		else if ( PC_TYPE == PETSC_PCAMG )
		{
			//PCGAMGSetNlevels( pc, 3) ;
			//PCMGSetCycleType( pc,PC_MG_CYCLE_V ) ;
			//PCMGSetCycleType( pc,PC_MG_CYCLE_W ) ;
			//cout<<"Test"<<endl;exit(1);
			PCSetType ( pc, PCGAMG) ;
			PCGAMGSetType ( pc, PCGAMGAGG ) ;
			//PCGAMGSetNSmooths ( pc, 1 ) ;
			PCGAMGSetNSmooths ( pc, 0 ) ;
			PCGAMGSetThresholdScale( pc, 0.08 ) ; // can be 0 to 0.08
			//PCGAMGSetThreshold( pc, 0.04, 1 ) ; // can be 0 to 0.08
			//close by KL
			PetscOptionsSetValue( NULL, "-mg_levels_ksp_type", "richardson" ) ; // gmres, richardson
			PetscOptionsSetValue( NULL, "-mg_levels_pc_type", "sor" ) ; // sor, asm
			PetscOptionsSetValue( NULL, "-mg_levels_ksp_max_it", "1" ) ;
			PetscOptionsSetValue( NULL, "-mg_coarse_ksp_type", "richardson" ) ;
			PetscOptionsSetValue( NULL, "-mg_coarse_pc_type", "sor" ) ; // sor, asm
			PetscOptionsSetValue( NULL, "-mg_coarse_ksp_max_it", "1" ) ;
			PetscOptionsSetValue( NULL, "-pc_gamg_sym_graph", "true" ) ;

		}

	}
	else
	{
		PCSetType ( pc , PCLU ) ;
	}

	KSPSetType( ksp, KSPBCGS );
	//KSPSetType( ksp,  KSPGMRES );
	KSPSetTolerances( ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );

	KSPSetFromOptions( ksp ) ;
	VecGetArray ( X, &x ) ;

}

void PETScSolver::reuse_preconditioner( bool flag )
{
	if ( flag )
		PCSetReusePreconditioner( pc, PETSC_TRUE  ) ;
	else
		PCSetReusePreconditioner( pc, PETSC_FALSE  ) ;
}

int PETScSolver::push_matrix ( int row, int number_of_entry, int *column, double *v )
{
	MatSetValues( A , 1, &row, number_of_entry, column, v,  INSERT_VALUES ) ;
	return 0 ;
}
int PETScSolver::Add_Entries( int row, int column, double v )
{

	MatSetValues( A , 1, &row, 1, &column, &v,  ADD_VALUES ) ;
	return 0 ;
}
int PETScSolver::ZeroEntriesMat()
{

	MatZeroEntries( A ) ;
	return 0 ;
}

int PETScSolver::push_source ( int row, double v )
{
	VecSetValue ( B, row, v, INSERT_VALUES  ) ;
	return 0 ;
}

int PETScSolver::add_source ( int row, double v )
{
	VecSetValue ( B, row, v, ADD_VALUES  ) ;
	return 0 ;
}

int PETScSolver::solve ( )
{
	

	PetscErrorCode ierr;

	VecAssemblyBegin( B ) ;
	VecAssemblyEnd( B ) ;

	MatAssemblyBegin ( A , MAT_FINAL_ASSEMBLY ) ;
	MatAssemblyEnd (  A , MAT_FINAL_ASSEMBLY ) ;

	//MatView (A, PETSC_VIEWER_STDOUT_WORLD );
	//VecView ( B, PETSC_VIEWER_STDOUT_WORLD ) ;

	VecRestoreArray( X, &x ) ;

	ierr = KSPSolve( ksp, B, X );
	//BOOST_ASSERT_MSG( ierr != PETSC_ERR_CONV_FAILED , "PETScSolver converge failed." );
	BOOST_ASSERT_MSG( ierr == 0, "PETScSolver converge failed." );
	//CHKERRQ(ierr);
	KSPGetIterationNumber( ksp, &its );
	total_it += its ;
	KSPSetInitialGuessNonzero( ksp, PETSC_TRUE) ;
	VecGetArray( X, &x ) ;

	return 0 ;
}
void PETScSolver::reuse_initialGuess( bool flag )
{
	if ( flag )
		KSPSetInitialGuessNonzero( ksp, PETSC_TRUE) ;
	else
		KSPSetInitialGuessNonzero( ksp, PETSC_FALSE) ;
}

/*					NONLINEAR SOLVER SNES					*/
void PETScSolverSnes::init( int gargc, char **gargv )
{
	int	i ;

	PetscInitialize( &gargc , &gargv, PETSC_NULL, PETSC_NULL ) ;
	SNESCreate( PETSC_COMM_WORLD, &snes ) ;
	/*	Create Vector for solution and nonlinear function  */
	VecCreate( PETSC_COMM_WORLD, &UnknownX ) ;
	VecSetSizes(UnknownX, PETSC_DECIDE, SnesUnknownNum ) ;
	VecSetFromOptions( UnknownX ) ;
	VecDuplicate( UnknownX, &Residue ) ;

	/*	Create Jacobian matrix data structure  */
	MatCreate( PETSC_COMM_WORLD, &Jacobian ) ;
	MatSetSizes( Jacobian, PETSC_DECIDE, PETSC_DECIDE, SnesUnknownNum, SnesUnknownNum ) ;
	MatSetFromOptions( Jacobian ) ;
	MatSetUp( Jacobian ) ;
}

int PETScSolverSnes::solve ( double tolerance, int maxiternum )
{
	SNESGetKSP( snes, &ksp ) ;
	KSPGetPC( ksp, &pc ) ;
	PCSetType( pc, PCNONE ) ;
	KSPSetTolerances( ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxiternum ) ;
	SNESSolve( snes, NULL, UnknownX ) ;
	return 0 ;
}

void PETScSolver::reset_matrix( )
{
	MatZeroEntries ( A );
}

void PETScSolver::reset_explicit ( )
{
	VecZeroEntries ( B );
}

// void PETScSolver::add_operator_A ( fvm::Laplacian &op )
// {
// 	int i, j ;

// 	int row ;
// 	int *column	;
// 	double *value;

// 	// Implicit part
// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row = op.entry_index[i][0] ;
// 		column = new int [ op.entry_index[i].size() ];
// 		value = new double [ op.entry_index[i].size() ];

// 		for ( j = 0 ; j < op.entry_index[i].size() ; j++ )
// 		{
// 			column[j] = op.entry_index[i][j]  ;
// 			value[j] = op.entry[i][j]  ;
// 			BOOST_ASSERT_MSG( value[j] == value[j], "NAN in matrix A (Laplacian operator)" );
// 		}
// 		MatSetValues( A, 1, &row, op.entry_index[i].size(), column, value, ADD_VALUES ) ;
// 		delete [] column ;
// 		delete [] value ;
// 	}
// }

// void PETScSolver::add_operator_B ( fvm::Laplacian &op )
// {
// 	int i, j ;
// 	int row , cid ;
// 	double source_term ;
// 	int *column	;
// 	double *value;

// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row =  op.entry_index[i][0] ;
// 		source_term = 0.;


// 		for ( j = 0 ; j < op.deferred_local_entry_index[i].size() ; j++ )
// 		{
// 			cid = op.deferred_local_entry_index[i][j] ;
// 			BOOST_ASSERT_MSG ( op.variable[cid] ==  op.variable[cid] , "NaN in deferred variable." );
// 			source_term += op.deferred_entry[i][j] * op.variable[cid] ;
// 		}

// 		for ( j = 0 ; j < op.deferred_BC_face_index[i].size() ; j++ )
// 		{
// 			BOOST_ASSERT_MSG ( op.BC_face_value[ op.deferred_BC_face_index[i][j] ] != NAN , "NaN in deferred boundary condition." );
// 			source_term += op.deferred_BC_face_coefficients[i][j] * op.BC_face_value[ op.deferred_BC_face_index[i][j] ] ;
// 		}

// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}

// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row =  op.entry_index[i][0] ;
// 		source_term = 0.;
// 		for ( j = 0 ; j < op.BC_face_index[i].size() ; j++ )
// 		{
// 			if ( op.BC_face_index[i][j] != -1 )
// 			{
// 				source_term += op.BC_face_coefficients[i][j] * op.BC_face_value[ op.BC_face_index[i][j] ] ;
// 			}

// 			if ( op.discontinued_index[i][j] != -1 )
// 			{
// 				source_term += op.discontinued_coefficient[i][j] * op.discontinued_face_value[ op.discontinued_index[i][j] ] ;
// 			}
// 		}
// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}
// }

// void PETScSolver::add_operator_A ( fvm::Divergence &op )
// {
// 	int i, j ;

// 	int row ;
// 	int *column	;
// 	double *value;

// 	// Implicit part
// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row = op.entry_index[i][0] ;
// 		column = new int [ op.entry_index[i].size() ];
// 		value = new double [ op.entry_index[i].size() ];

// 		for ( j = 0 ; j < op.entry_index[i].size() ; j++ )
// 		{
// 			column[j] = op.entry_index[i][j]  ;
// 			value[j] = op.entry[i][j]  ;
// 		}
// 		MatSetValues( A, 1, &row, op.entry_index[i].size(), column, value, ADD_VALUES ) ;
// 		delete [] column ;
// 		delete [] value ;
// 	}


// 	//MatAssemblyBegin ( A , MAT_FINAL_ASSEMBLY ) ;
// 	//MatAssemblyEnd (  A , MAT_FINAL_ASSEMBLY ) ;
// 	//MatView (A, PETSC_VIEWER_STDOUT_WORLD );
// }

// void PETScSolver::add_operator_B ( fvm::Divergence &op )
// {
// 	int i, j ;
// 	int row , cid ;
// 	double source_term ;
// 	int *column	;
// 	double *value;

// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row =  op.entry_index[i][0] ;
// 		source_term = 0.;

// 		for ( j = 0 ; j < op.explicit_local_entry_index[i].size() ; j++ )
// 		{
// 			cid = op.explicit_local_entry_index[i][j] ;
// 			source_term += op.explicit_entry [i][j] * op.variable[cid] ;
// 		}

// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}

// }

// void PETScSolver::add_operator_A ( TSM &op )
// {
// 	int i, j ;

// 	int row ;
// 	int *column	;
// 	double *value;

// 	if ( op.stepping_type == FORWARD_EULER ) return ;

// 	// Implicit part
// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row = op.entry_index[i][0] ;
// 		column = new int [ op.entry_index[i].size() ];
// 		value = new double [ op.entry_index[i].size() ];

// 		for ( j = 0 ; j < op.entry_index[i].size() ; j++ )
// 		{
// 			column[j] = op.entry_index[i][j]  ;
// 			value[j] = op.entry[i][j]  ;
// 		}
// 		MatSetValues( A, 1, &row, op.entry_index[i].size(), column, value, ADD_VALUES ) ;
// 		delete [] column ;
// 		delete [] value ;
// 	}


// }

// void PETScSolver::add_operator_B ( TSM &op )
// {

// 	int i, j ;
// 	int row , cid ;
// 	double source_term ;
// 	int *column	;
// 	double *value;

// 	if ( op.stepping_type == CRANK_NICOLSON  )
// 	{
// 		for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 		{
// 			row =  op.entry_index[i][0] ;
// 			source_term = 0.;

// 			for ( j = 0 ; j < op.explicit_local_entry_index[i].size() ; j++ )
// 			{
// 				cid = op.explicit_local_entry_index[i][j] ;
// 				source_term += - op.explicit_entry [i][j] * op.pre_variable[cid] ;
// 			}

// 			BOOST_ASSERT_MSG ( source_term != NAN, "NaN in time-stepping source from CRANK NICOLSON scheme." );
// 			VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 		}
// 	}


// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row = op.entry_index[i][0] ;
// 		source_term = 0.;

// 		for ( j = 0 ; j < op.deferred_local_entry_index[i].size() ; j++ )
// 		{
// 			cid = op.deferred_local_entry_index[i][j] ;
// 			source_term += op.deferred_entry[i][j] * op.variable[cid] ;
// 		}

// 		for ( j = 0 ; j < op.deferred_BC_face_index[i].size() ; j++ )
// 		{
// 			source_term += op.deferred_BC_face_coefficients[i][j] * op.BC_face_value[ op.deferred_BC_face_index[i][j] ] ;
// 		}

// 		BOOST_ASSERT_MSG ( source_term != NAN , "NaN in time-stepping deferred source." );
// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}

// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row =  op.entry_index[i][0] ;
// 		source_term = 0.;
// 		for ( j = 0 ; j < op.BC_face_index[i].size() ; j++ )
// 		{
// 			if ( op.BC_face_index[i][j] != -1 )
// 			{
// 				if ( op.stepping_type == CRANK_NICOLSON )
// 				{
// 					source_term -= op.BC_face_coefficients[i][j] * op.BC_face_value[ op.BC_face_index[i][j] ];
// 				} else
// 				{
// 					source_term -= op.BC_face_coefficients[i][j] * op.BC_face_value[ op.BC_face_index[i][j] ] ;
// 				}
// 			}
// 		}
// 		BOOST_ASSERT_MSG ( source_term != NAN , "NaN in time-stepping BC source." );
// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}

// 	// Given source
// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row =  op.entry_index[i][0] ;
// 		source_term = op.dt * op.source[i] ;
// 		BOOST_ASSERT_MSG ( source_term != NAN , "NaN in time-stepping from the given source." );
// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}

// 	// Actually this is the variable at previous time step
// 	for ( i = 0 ; i < op.local_entry_number() ; i++ )
// 	{
// 		row =  op.entry_index[i][0] ;
// 		source_term = op.source_pre_step[i] ;

// 		BOOST_ASSERT_MSG ( source_term != NAN , "NaN in time-stepping source from previous time step." );
// 		VecSetValue ( B, row, source_term, ADD_VALUES ) ;
// 	}

// }

#include <petscksp.h>
#include <petscsnes.h>
using namespace std;

#if !defined(__PETSC_SOLVER_H)
#define __PETSC_SOLVER_H

#define PETSC_PCASM 	1
#define PETSC_PCAMG		2

class PETScSolver
{
	public:
		PETScSolver (  ) ;
		~PETScSolver (  ) ;
		MPI_Comm	comm ;
		int		mpi_id ;
		int		local_unknown_number ;
		int		global_unknown_number ;
		int     dia_nz, off_nz;
		int     *dia_iNum_nz, *off_iNum_nz;
		int		low, high ;
		void	set_local_unknown_number ( int ) ;
		void	init( ) ;
		void	init( int ) ;
		void	init( double , int ) ;
		void    set_nz(int Max_nz, int *d_nnz, int *o_nnz);
		int		push_matrix ( int, int, int * , double* ) ;
		int		push_source ( int, double ) ;
		int 	Add_Entries( int, int, double ) ;
		int		add_source ( int, double ) ;
		int		solve() ;
		void	close() ;
		int 	ZeroEntriesMat() ;
		void	reuse_preconditioner( bool ) ;

		void 	reset_matrix() ;
		void 	reset_explicit() ;
		void 	reuse_initialGuess( bool flag ) ;

		unsigned long total_it ;
		static bool flagPetscInitialize ;

		/* PETSc A x = B */
		int		asm_overlap ;
		Mat		A ;
		Vec		X ;
		Vec		B ;
		KSP		ksp ;
		PC		pc ;
		PetscScalar	*b, *x ;
		int	its;
	private:

		bool flag_local_number_initialization ;

};

class PETScSolverSnes
{
	public:
		int		mpi_id ;
		int		SnesUnknownNum ;
		void	init( int, char ** ) ;
		int		solve( double, int ) ;
		double bb;
		/* PETSc A x = B */
		SNES           snes;         /* nonlinear solver context */
		KSP            ksp;          /* linear solver context */
		PC             pc;           /* preconditioner context */
		Vec            UnknownX,Residue;          /* solution, residual vectors */
		PetscInt       its;
		Mat            Jacobian;            /* Jacobian matrix */
		PetscMPIInt    size;
		PetscBool      flg ;
	private:

};





#endif


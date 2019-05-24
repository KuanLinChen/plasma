#include <mpi.h>
#include "petscsys.h" 
#include "config_structure.hpp"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "solver_poisson.hpp"
#include "solver_drift_diffusion.hpp"
#include "solver_drift_diffusion.hpp"
 #include "post_structure.hpp"
//#include "ultraMPP.h"
// #include "solver_poisson_explicit.hpp"
// #include "solver_energy_density.hpp"
// #include "solver_fluid_model.hpp"
// #include "solver_navier_stokes.hpp"
// #include "variable_structure_NS.hpp"
// #include "PETSc_solver.h"

using namespace std ;
int mpi_size; /*!< \brief The number of processes. */
int	mpi_rank ; /*!< \brief The rank(id) of the process. */
ultraMPP plasma ; 
//map<string,int> var_name ;/*!< \brief map of variable names. It contain the nama of the variable and index associated to this name */

int main( int argc, char * argv[] )
{
	/* the 'plasma' is global variable, it is define in the 'PFM.hpp' file */
	plasma.initial( argc, argv, &mpi_rank, &mpi_size ) ;

	/* read the case input file path */
	if ( argv[1] == NULL ) {
		if( mpi_rank ==0 ) cout<<"Pls using ./main  case_path "<<endl;
		exit(1) ;
	} 

	/*--- Configuration ---*/
		cout<<"Configuration"<<endl;
		boost::shared_ptr<CConfig> Config ;
		Config = boost::shared_ptr<CConfig> ( new CConfig ) ;
		Config->Init( argv[1] ) ;
		//Note: the number of matrix solvers are poisson equation + number of total sepcies + electron energy.
		//TODO: this part should be fix. becaeuse sometime you need solve the ion temperature.


		//TODO: this part should be fix....
		Config->MeshFile = "input.json" ;
		plasma.load_mesh( argv[1] + Config->MeshFile ) ;
		/* potential + numSpecies*Continuity + electron energy */
		int numMatrixSolver = 1 + Config->TotalSpeciesNum  + 1  ;
		plasma.set_matrix_number( numMatrixSolver ) ;

	/*--- Link ultraMPP & Read Mesh file ---*/
	 	boost::shared_ptr<CDomain> mesh ;
	 	mesh = boost::shared_ptr<CDomain> ( new CDomain ) ;
		//cout<<"Link ultraMPP"<<endl;
	 	mesh->link_ultraMPP() ;
		//cout<<"SetupCellFaceType"<<endl;
		mesh->SetupCellFaceType() ;
		//cout<<"BulidCellStructure"<<endl;
	 	mesh->BulidCellStructure() ;

	/*--- Solution Variables ---*/
	 	//cout<<" Solution Variables"<<endl;
		boost::shared_ptr<CVariable> Var ;
	 	Var = boost::shared_ptr<CVariable> ( new CVariable ) ;
	 	Var->Init( mesh , Config ) ;
		Var->Calculate_LSQ_Coeff_Scalar( mesh ) ;

	/*--- Poisson solver ---*/
		//cout<<" Poisson solver"<<endl;
		boost::shared_ptr<CPoisson> poisson_solver ;
		poisson_solver = boost::shared_ptr<CPoisson> ( new CPoisson ) ;
		poisson_solver->Init( mesh, Config ) ;

	/*--- Dirft-diffusion solver ---*/
		int DriftDiffusionNum=0, FullEqnNum=0 ;/*!< \brief Number of drift-diffusion eqn. & Full eqn. */
		int SpeciesType=0 ;
		boost::shared_ptr<CDriftDiffusion> *continuity_solver ;
		for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {
			SpeciesType = Config->Species[ jSpecies ].Type ;
			if ( Config->Equation[ SpeciesType ].Equation == 0 ){
				DriftDiffusionNum++ ;
			} else if ( Config->Equation[ SpeciesType ].Equation > 0 ) {
				FullEqnNum++ ;
			}
		} if(mpi_rank==0) cout<<"Number of D-D eqn: "<< DriftDiffusionNum <<", Number of full eqn: "<<FullEqnNum<<endl;

		if( DriftDiffusionNum > 0 ){			
			continuity_solver = new boost::shared_ptr<CDriftDiffusion> [ DriftDiffusionNum ] ;
			DriftDiffusionNum = 0 ;
			for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {
				SpeciesType = Config->Species[ jSpecies ].Type ;
				if ( Config->Equation[ SpeciesType ].Equation == 0 ){
					continuity_solver[ DriftDiffusionNum ] = boost::shared_ptr<CDriftDiffusion> ( new CDriftDiffusion ) ;
					continuity_solver[ DriftDiffusionNum ]->Init( mesh, Config, Config->Species[ jSpecies ].Index ) ;
					DriftDiffusionNum++ ;
				}
			}
		}

	/*--- full equations solver ---*/
		// boost::shared_ptr<CFluidModel> *fluid_model_solver ;
		// if( FullEqnNum > 0 ){			
		// 	fluid_model_solver = new boost::shared_ptr<CFluidModel> [ FullEqnNum ] ;
		// 	FullEqnNum = 0 ;
		// 	for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {
		// 		SpeciesType = Config->Species[ jSpecies ].Type ;
		// 		if ( Config->Equation[ SpeciesType ].Equation > 0 ){
		// 			fluid_model_solver[ FullEqnNum ] = boost::shared_ptr<CFluidModel> ( new CFluidModel ) ;
		// 			fluid_model_solver[ FullEqnNum ]->Init( mesh, Config, Config->Species[ jSpecies ].Index ) ;
		// 			FullEqnNum++ ;
		// 		}
		// 	}
		// }

	/*--- Energy density Eqn. w/ Drift-diffusion solver ---*/
		// boost::shared_ptr<CEnergyDensity> electron_energy_solver=NULL ;
		// if( DriftDiffusionNum > 0 ){
		// 	electron_energy_solver = boost::shared_ptr<CEnergyDensity> ( new CEnergyDensity ) ;
		// 	electron_energy_solver->Init( mesh, Config, 0 ) ;
		// }
	/*--- Output ---*/
		boost::shared_ptr<CPost> post ;
		post = boost::shared_ptr<CPost> ( new CPost ) ;

// 		map< int, CElectrical>::iterator Iter;

 		bool WRT_CYC_AVG = false, WRT_INS = false ;
 		bool     MON_CYC = false, MON_INS = false ;

		poisson_solver->Solve( mesh, Config, Var ) ;
 		Var->UpdateSolution( mesh ) ; 
 		post->OutputFlow( mesh, Config, Var, 0, 0 ) ;
 		for ( int MainCycle = 1 ; MainCycle < Config->ExitCycle ; MainCycle ++  ) {

 			if( MainCycle%Config->WRT_Cycle_Freq == 0 and Config->Average_Switch ) WRT_CYC_AVG = true ;
			else WRT_CYC_AVG = false ;

			if( MainCycle%Config->MON_Cycle_Freq == 0 ) MON_CYC = true ;
			else MON_CYC = false ;

// 			/*--- If output cycle average, then fist u need to reset zero in average variable. ---*/
			Var->ResetAvgZero_Electrode( mesh, Config ) ;
 			if( WRT_CYC_AVG ) {
 				Var->ResetAvgZero( mesh, Config ) ;
			}

			/*--- Main step loop in eack cycle. ---*/
			/*--- Main step loop in eack cycle. ---*/
			/*--- Main step loop in eack cycle. ---*/
 			for ( int MainStep = 0 ; MainStep < Config->StepPerCycle ; MainStep++ ) {
 			//for ( int MainStep = 0 ; MainStep < 4 ; MainStep++ ) {

 				if( WRT_CYC_AVG and MainStep%(Config->StepPerCycle/Config->WRT_Insta_Freq) == 0  ) WRT_INS = true ;
 				else  WRT_INS 	= false ;

 				if( MON_CYC and MainStep%(Config->StepPerCycle/Config->MON_Insta_Freq) == 0  ) MON_INS = true ;
 				else  MON_INS 	= false ;

 				if( mpi_rank == MASTER_NODE and MON_CYC and MON_INS ){
 					cout<<"MainCycle: "<<MainCycle<<"\t"<<"MainStep: "<<MainStep<<"\t"<<"Voltage: "<<Var->Volt<<"  [V]"<<"\t"<<"PhysicalTime: "<<Var->PhysicalTime<<"  [s]"<<endl ; 
 					//cout<<"Poissn ksp iter : "<< plasma.get_convergence_reason( 0 ) <<endl ;
 					for ( int iEqn = 0 ; iEqn < DriftDiffusionNum ; iEqn++ ) {
 						//cout<<"Continuity["<<iEqn<<"] ksp iter : "<<plasma.get_convergence_reason( iEqn+1 )<<endl ;
						
 					}
 				}
 				/*--- Copy solution (n+1) step -> (n) step ---*/ 
				//cout<<"A"<<endl;
 				Var->UpdateSolution( mesh ) ; 
				//cout<<"B"<<endl;

 				/*--- Update the transport coeffieients and rate constants ---*/ 
				//cout<<"C"<<endl;
 				Var->ChemistryUpdate( mesh, Config ) ; 
				//cout<<"D"<<endl;

				/*--- Solve for potential and electric field. ---*/
				//cout<<"E"<<endl;
 				poisson_solver->Solve( mesh, Config, Var ) ;
				//cout<<"F"<<endl;


				/*--- Solve number density for next time step (n+1). ---*/
				for ( int iEqn = 0 ; iEqn < DriftDiffusionNum ; iEqn++ ) {
				//cout<<"G"<<iEqn<<endl;
					continuity_solver[ iEqn ]->Solve( mesh, Config, Var ) ;
				//cout<<"H"<<iEqn<<endl;
				}

 				if( WRT_INS ){			
 					post->OutputFlow( mesh, Config, Var, MainCycle, MainStep ) ;
 				}
 				//post->OutputFlow( mesh, Config, Var, MainCycle, MainStep ) ;
// 				/*--- Calculate cycle average data. ---*/
// 				if( WRT_CYC_AVG ) Var->AddAverage( mesh, Config ) ;

// 				Var->AddAverage_Electrode( mesh, Config ) ;

 				Var->PhysicalTime += Var->Dt ;
			}//End Main Step

			
// 			/*--- Output cycle average data ---*/
// 			if ( MainCycle % 20 == 0 ){
// 				for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
//     				if ( Iter->first == POWER and fabs(Var->I_AvgPowerElectrode) > 1.E-6 ) {
// 	    				//Iter->second.BiasVoltage += (Var->I_AvgPowerElectrode*Var->Dt*Config->StepPerCycle)/(500*1.0E-12)*0.5 ;
//     					BiasVoltage += (Var->I_AvgPowerElectrode*Var->Dt*Config->StepPerCycle)/(500*1.0E-12)*0.5 ;
//     				}
// 				}
// 			}

// 			if( mpi_rank == MASTER_NODE and MON_CYC ){
// 				cout<<"BiasVoltage: "<< BiasVoltage <<endl ;
// 				cout<<"Power Current[A] : "<< Var->I_AvgPowerElectrode <<endl ;
// 			}

// 			if( WRT_CYC_AVG ){			

// 				post->OutputAverageFlow( Config, Var, MainCycle ) ;

// 				FileOutput.close();
// 				FileOutput.clear();
// 			}

// 			Config->Cycle ++ ;
 		}//End Main Cycle

	// PetscFinalize();
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Finalize() ;
	plasma.end_MPP() ;
	return 0 ;
}






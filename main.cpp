#include <mpi.h>
#include "petscsys.h" 
#include "config_structure.hpp"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "solver_poisson.hpp"
#include "solver_drift_diffusion.hpp"
#include "post_structure.hpp"
#include "solver_energy_density.hpp"
#include "solver_fluid_model.hpp"
// #include "solver_navier_stokes.hpp"
// #include "variable_structure_NS.hpp"
// #include "PETSc_solver.h"
#define Debug false
using namespace std ;
int mpi_size, /*!< \brief The number of processes. */
		mpi_rank ;/*!< \brief The rank(id) of the process. */
int nDim ; /*!< \brief The dimension of the problem. */
bool WRT_CYC_AVG = false ; /*!< \brief Trigger for write cycle averaged data. */
bool WRT_INS = false ; /*!< \brief Trigger for write instantaneous averaged data. */
bool MON_CYC = false ; /*!< \brief Trigger for monitor cycle averaged data. */
bool MON_INS = false ; /*!< \brief Trigger for monitor instantaneous averaged data. */
ultraMPP plasma ; 

int	gargc2 ;
char **gargv2 ;

int main( int argc, char * argv[] )
{
	/* the 'plasma' is global variable, it is define in the 'PFM.hpp' file */
	gargc2	=	argc ;
	gargv2	=	argv ;

	plasma.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;

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
	/*--- This module will be delets after some modification --*/
	boost::shared_ptr<CDomain> mesh ;
	mesh = boost::shared_ptr<CDomain> ( new CDomain ) ;
	mesh->BulidCellStructure() ;

	/*--- Solution Variables ---*/
	boost::shared_ptr<CVariable> Var ;
	Var = boost::shared_ptr<CVariable> ( new CVariable ) ;
	Var->Init( mesh , Config ) ;
	Var->Calculate_LSQ_Coeff_Scalar( mesh ) ;

	/*--- Poisson solver ---*/
	boost::shared_ptr<CPoisson> poisson_solver ;
	poisson_solver = boost::shared_ptr<CPoisson> ( new CPoisson ) ;
	poisson_solver->Init( mesh, Config ) ;

	/*--- Dirft-diffusion solver ---*/
	int DriftDiffusionNum=0, FullEqnNum=0 ;
	int SpeciesType=0 ;
	boost::shared_ptr<CDriftDiffusion> *continuity_solver ;


	/* First: Count the number of D-D equation and full equation. */
	for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {
		SpeciesType = Config->Species[ jSpecies ].Type ;
		if ( Config->Equation[ SpeciesType ].Equation == 0 ){
			DriftDiffusionNum++ ;
		} else if ( Config->Equation[ SpeciesType ].Equation > 0 ) {
			FullEqnNum++ ;
		}
	} if(mpi_rank==0) cout<<"Number of D-D eqn: "<< DriftDiffusionNum <<", Number of full eqn: "<<FullEqnNum<<endl<<endl;


	/* Second: Create the D.-D. equation modules */
	if( DriftDiffusionNum > 0 ){			
		continuity_solver = new boost::shared_ptr<CDriftDiffusion> [ DriftDiffusionNum ] ;
		DriftDiffusionNum = 0 ;
		for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {
			SpeciesType = Config->Species[ jSpecies ].Type ;
			if ( Config->Equation[ SpeciesType ].Equation == 0 ){
				continuity_solver[ DriftDiffusionNum ] = boost::shared_ptr<CDriftDiffusion> ( new CDriftDiffusion ) ;
				continuity_solver[ DriftDiffusionNum ]->Init( mesh, Config, Config->Species[ jSpecies ].Index ) ;
				DriftDiffusionNum++ ;
				if(mpi_rank==0) cout<<"drift-diffustion eqn. index: "<<Config->Species[ jSpecies ].Index<<endl<<endl;
			}
		}//End iSpecies
	}//End Drift-Diffusion

	/*--- full equations solver ---*/
		boost::shared_ptr<CFluidModel> *fluid_model_solver ;
		if( FullEqnNum > 0 ){			
			fluid_model_solver = new boost::shared_ptr<CFluidModel> [ FullEqnNum ] ;
			FullEqnNum = 0 ;
			for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {
				SpeciesType = Config->Species[ jSpecies ].Type ;
				if ( Config->Equation[ SpeciesType ].Equation > 0 ){
					fluid_model_solver[ FullEqnNum ] = boost::shared_ptr<CFluidModel> ( new CFluidModel ) ;
					fluid_model_solver[ FullEqnNum ]->Init( mesh, Config, Config->Species[ jSpecies ].Index ) ;
					FullEqnNum++ ;
				}
			}
		}

	/*--- Energy density Eqn. w/ Drift-diffusion solver ---*/
		boost::shared_ptr<CEnergyDensity> electron_energy_solver ;
		electron_energy_solver = boost::shared_ptr<CEnergyDensity> ( new CEnergyDensity ) ;
		electron_energy_solver->Init( mesh, Config, 0 ) ;

	/*--- post-processing module ---*/
		boost::shared_ptr<CPost> post ;
		post = boost::shared_ptr<CPost> ( new CPost ) ;

 	/* first solve potential and electric field as initial. */
		poisson_solver->Solve( mesh, Config, Var ) ;
 		Var->UpdateSolution( mesh ) ; 
 		Var->ChemistryUpdate( mesh, Config ) ; 
 		post->OutputFlow( mesh, Config, Var, 0, 0 ) ;

		ofstream        FileOutput, PCB_FileOutput ;
		map< int, CElectrical>::iterator Iter;

		if ( mpi_rank == MASTER_NODE ){
			PCB_FileOutput.open( "PCB.dat" , ios::out | ios::trunc ) ;
			PCB_FileOutput<<"VARIABLES=\"Time\", \"PowerAbs [W]\", \"Bias Voltage [V]\" , \"Residue\""<<endl ;
		}
		double BiasVoltage=0.0 ;

 		/*--- Main Cycle ---*/
 		for ( int MainCycle = 1 ; MainCycle < Config->ExitCycle ; MainCycle ++  ) {

 			if( MainCycle%Config->WRT_Cycle_Freq == 0 and Config->Average_Switch ) WRT_CYC_AVG = true ;
			else WRT_CYC_AVG = false ;

			if( MainCycle%Config->MON_Cycle_Freq == 0 ) MON_CYC = true ;
			else MON_CYC = false ;

			Var->ResetAvgZero_Electrode( mesh, Config ) ;
 			Var->ResetAvgZero_PowerAbs( mesh, Config ) ;


			/* If output cycle average, then fist you need to reset zero in average variable. */
			if( WRT_CYC_AVG ) {

				Var->ResetAvgZero( mesh, Config ) ;

				if ( mpi_rank == MASTER_NODE ) {		

					FileOutput.open( "I_V_"+to_string(MainCycle)+".dat" , ios::out | ios::trunc ) ;
					FileOutput<<"VARIABLES=\"step\", \"voltage [V]\", \"I<sub>p</sub>\" , \"Disp_I<sub>p</sub>\"" ;
					for (int iSpecies = 0 ; iSpecies< Config->ChargeSpeciesNum ; iSpecies ++ ){
						FileOutput<<"\"Cond["+to_string(iSpecies)+"]_I<sub>p</sub>"+"\"" ;
					}
					FileOutput<<"\"I<sub>g</sub>\" , \"Disp_I<sub>g</sub>\"" ;
					for (int iSpecies = 0 ; iSpecies< Config->ChargeSpeciesNum ; iSpecies ++ ){
						FileOutput<<"\"Cond["+to_string(iSpecies)+"]_I<sub>g</sub>"+"\"" ;
					}
					FileOutput<<endl;

				}
			}

			/*--- Main step loop in eack cycle. ---*/
 			for ( int MainStep = 0 ; MainStep < Config->StepPerCycle ; MainStep++ ) {

 				if( WRT_CYC_AVG and MainStep%(Config->StepPerCycle/Config->WRT_Insta_Freq) == 0  ) WRT_INS = true ;
 				else  WRT_INS 	= false ;

 				if( MON_CYC and MainStep%(Config->StepPerCycle/Config->MON_Insta_Freq) == 0  ) MON_INS = true ;
 				else  MON_INS 	= false ;

 				if( mpi_rank == MASTER_NODE and MON_CYC and MON_INS ){
 					cout<<"MainCycle: "<<MainCycle<<"\t"<<"MainStep: "<<MainStep<<"\t"<<"Voltage: "<<Var->Volt<<"  [V]"<<"\t"<<"PhysicalTime: "<<Var->PhysicalTime<<"  [s]"<<endl ; 
 				}

 				/* Update solution: Copy solution (n+1) step -> (n) step */ 
 				Var->UpdateSolution( mesh ) ; 
	 				#if (Debug == true )
	 					PetscPrintf( PETSC_COMM_WORLD, "UpdateSolution done...\n" ) ;
	 				#endif
				
 				/* Update the rate constants & transport coefficients. */
 				Var->ChemistryUpdate( mesh, Config ) ; 
	 				#if (Debug == true )
	 					PetscPrintf( PETSC_COMM_WORLD, "ChemistryUpdate done...\n" ) ;
	 				#endif



				/*--- Predtic the ion number density for poisson eqn. (Only for full eqn.) ---*/
				for ( int iEqn = 0 ; iEqn < FullEqnNum ; iEqn++ ) {
					if ( fluid_model_solver[ iEqn ]->SpeciesType == ION ) {
						Var->Alpha_Beta(  mesh, Config ) ;
						fluid_model_solver[ iEqn ]->Solve_Continuity( mesh, Config, Var ) ;
					}
				}
	 				#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, " Predtic the ion number density...\n" ) ;
	 				#endif				


				/* Solve for potential and electric field. */
 				poisson_solver->Solve( mesh, Config, Var ) ;
	 				#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, "poisson_solver done...\n" ) ;
	 				#endif

				/* Solve number density using D-D approximation for next time step (n+1). */
				for ( int iEqn = 0 ; iEqn < DriftDiffusionNum ; iEqn++ ) {
					continuity_solver[ iEqn ]->Solve( mesh, Config, Var ) ;
				}

				/*--- Solve number density and flux for next time step (k+1). ---*/
				for ( int iEqn = 0 ; iEqn < FullEqnNum ; iEqn++ ) {
					if ( fluid_model_solver[ iEqn ]->SpeciesType != ION ) {
						fluid_model_solver[ iEqn ]->Solve_Continuity( mesh, Config, Var ) ;
					}

					fluid_model_solver[ iEqn ]->Solve_Momentum( mesh, Config, Var ) ;
				}

				Var->CalculateElectrodeCurrent( mesh, Config ) ;


				/* Solve energy density for next time step (n+1). */
				if ( Config->PFM_Assumption == "LMEA" and DriftDiffusionNum > 0 ) {
					electron_energy_solver->Solver( mesh, Config, Var ) ;
				}

				/*--- Output I-V data ---*/
				if ( mpi_rank == MASTER_NODE and WRT_CYC_AVG ) {

					FileOutput <<(double)MainStep/(Config->StepPerCycle-1)<<"\t"
							   <<Var->Volt<<"\t"
							   <<(-1.0)*Var->I_PowerElectrode_global_sum<<"\t"
							   <<(-1.0)*Var->Disp_PowerElectrode_global_sum<<"\t" ;
					for (int iSpecies = 0 ; iSpecies< Config->ChargeSpeciesNum ; iSpecies ++ ){
						FileOutput <<(-1.0)*Var->CondI_PowerElectrode_global_sum[iSpecies]<<"\t" ;
					}

					FileOutput<<(-1.0)*Var->I_GroundElectrode_global_sum<<"\t" ;
					FileOutput<<(-1.0)*Var->Disp_GroundElectrode_global_sum<<"\t" ;
					for (int iSpecies = 0 ; iSpecies< Config->ChargeSpeciesNum ; iSpecies ++ ){
						FileOutput <<(-1.0)*Var->CondI_GroundElectrode_global_sum[iSpecies]<<"\t" ;
					}
					FileOutput<<endl ;

				}

				/* Output instantaneous flow field data */
 				if( WRT_INS ) 
 					post->OutputFlow( mesh, Config, Var, MainCycle, MainStep ) ;

				/* Calculate cycle average data. */
				if( WRT_CYC_AVG ) 
					Var->AddAverage( mesh, Config ) ;

 				Var->PhysicalTime += Var->Dt ;

			}//End Main Step


			if ( MainCycle % 50 == 0 ){
				for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
    				if ( Iter->first == POWER and fabs(Var->I_AvgPowerElectrode) > 1.E-6 ) {
	    				//Iter->second.BiasVoltage += (Var->I_AvgPowerElectrode*Var->Dt*Config->StepPerCycle)/(500*1.0E-12)*0.5 ;
    					//BiasVoltage += (Var->I_AvgPowerElectrode*Var->Dt*Config->StepPerCycle)/(500*1.0E-12)*0.5 ;
    				}
				}
			}


			if( mpi_rank == MASTER_NODE and MON_CYC ){
				cout<<"BiasVoltage: "<< BiasVoltage <<endl ;
				cout<<"Power Current[A] : "<< Var->I_AvgPowerElectrode <<endl ;
				cout<<"POWER [W]: "<<Var->AvgPowerAbs <<endl;
			}

			if( mpi_rank == MASTER_NODE )
			PCB_FileOutput<<Var->PhysicalTime<<"\t"<<Var->AvgPowerAbs<<"\t"<<BiasVoltage<<"\t"<<fabs(Var->I_AvgPowerElectrode)<<endl;

			/* Output cycle average data */
			if( WRT_CYC_AVG ) {			
				post->OutputAverageFlow( Config, Var, MainCycle ) ;

				FileOutput.close();
				FileOutput.clear();
			}

// 			Config->Cycle ++ ;
 		}//End Main Cycle

	return 0 ;
}






#include <mpi.h>
#include "petscsys.h" 
#include "config_structure.hpp"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "solver_poisson.hpp"
#include "solver_FD_maxwell.hpp"
#include "solver_drift_diffusion.hpp"
#include "post_structure.hpp"
#include "solver_energy_density.hpp"
#include "solver_fluid_model.hpp"
#define Debug false


#define time_monitor true 
#include "time.h"
#define poisson_test false

using namespace std ;
int mpi_size, /*!< \brief The number of processes. */
		mpi_rank ;/*!< \brief The rank(id) of the process. */
int nDim ; /*!< \brief The dimension of the problem. */
bool WRT_CYC_AVG = false ; /*!< \brief Trigger for write cycle averaged data. */
bool WRT_INS = false ; /*!< \brief Trigger for write instantaneous averaged data. */
bool MON_CYC = false ; /*!< \brief Trigger for monitor cycle averaged data. */
bool MON_INS = false ; /*!< \brief Trigger for monitor instantaneous averaged data. */

/* ultraMPP, these variables are global variable, it is define in the 'PFM.hpp' file  */
ultraMPP plasma  ;/*!< \brief UltraMPP main object. */

// ICP
#if (FDMaxwell == true ) 
ultraMPP FDMaxwell_Re  ;
ultraMPP FDMaxwell_Im  ;
ultraMPP FDMaxwell_coupled_eqs  ;
#endif

json &json_bc_setting ;/*!< \brief Boundary condition informations. */
json &json_cell_setting;/*!< \brief Cell condition informations. */  

map<string,double> cell_parameter ;
map<string,double> face_parameter ;

map<string,int>    MPP_face_tag ;
map<string,int>    MPP_cell_tag ;

#if (FDMaxwell == true )
map<string,int>    FVFD_face_tag ;
map<string,int>    FVFD_cell_tag ;
#endif

map<int,int> face_type ;
map<int,int> cell_type ;

map<int, string>	type_typename ;
map<string, int>	typename_type ;

int	gargc2 ;
char **gargv2 ;

int main( int argc, char * argv[] )
{
	gargc2	=	argc ;	gargv2	=	argv ;

	/* Read the case input file path */
		if ( argv[1] == NULL ) {
			if( mpi_rank ==0 ) cout<<"Pls using ./main  case_path "<<endl; exit(1) ;
		} 

	/* Configuration */
		cout<<"Configuration"<<endl;
		boost::shared_ptr<CConfig> Config ;
		Config = boost::shared_ptr<CConfig> ( new CConfig ) ;
		Config->Init( argv[1] ) ;

	/* Initial the ultraMPP object. */
		plasma.set_linear_solver_library("PETSC");
		//plasma.apply_linear_solver_setting();
		plasma.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;
		string MPPFile = "input.json" ;
		plasma.load_mesh( argv[1] + MPPFile ) ;
		
	/* Initial the FVFD object. */
	#if (FDMaxwell == true )
        FDMaxwell_Re.set_linear_solver_library("PETSC");
        FDMaxwell_Im.set_linear_solver_library("PETSC");
        FDMaxwell_coupled_eqs.set_linear_solver_library("PETSC");
    	FDMaxwell_Re.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;
    	FDMaxwell_Im.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;
    	FDMaxwell_coupled_eqs.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;
    	FDMaxwell_Re.set_N_variable_number( 1 ) ;
    	FDMaxwell_Im.set_N_variable_number( 1 ) ;
    	FDMaxwell_coupled_eqs.set_N_variable_number( 2 ) ;
		string MPPFile_FVFD = "FD_maxwell.json" ;
    	FDMaxwell_Re.load_mesh( argv[1] + MPPFile_FVFD ) ;
    	FDMaxwell_Im.load_mesh( argv[1] + MPPFile_FVFD ) ;
    	FDMaxwell_coupled_eqs.load_mesh( argv[1] + MPPFile_FVFD ) ;
	#endif	
	/*--- This module will be delets after some modification --*/
		boost::shared_ptr<CDomain> mesh ;
		mesh = boost::shared_ptr<CDomain> ( new CDomain ) ;
		mesh->BulidCellStructure() ;
		mesh->Init();


	//Note: the number of matrix solvers are poisson equation + number of total sepcies + electron energy.
	//TODO: this part should be fix. becaeuse sometime you need solve the ion temperature.

	//TODO: this part should be fix....
	/* potential + numSpecies*Continuity + electron energy */
	int numMatrixSolver = 1 + Config->TotalSpeciesNum  + 1  ;


	/*--- Solution Variables ---*/
	boost::shared_ptr<CVariable> Var ;
	Var = boost::shared_ptr<CVariable> ( new CVariable ) ;
	Var->Init( mesh , Config ) ;


	/*--- Poisson solver ---*/
	boost::shared_ptr<CPoisson> poisson_solver ;
	poisson_solver = boost::shared_ptr<CPoisson> ( new CPoisson ) ;
	poisson_solver->Init( Config ) ;
	
	/*--- Frequency domain Maxwell equation solver ---*/
	#if (FDMaxwell == true )
	boost::shared_ptr<CFDMaxwell> FD_maxwell_solver ;
	FD_maxwell_solver = boost::shared_ptr<CFDMaxwell> ( new CFDMaxwell ) ;
	FD_maxwell_solver->Init( Config,Var ) ;
	#endif

	/*--- Dirft-diffusion solver ---*/
	int DriftDiffusionNum=0, FullEqnNum=0 ;
	int SpeciesType=0 ;

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
	boost::shared_ptr<CDriftDiffusion> *continuity_solver ;
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
	 	poisson_solver->SOLVE( Config, Var ) ;
	 	Var->UpdateSolution( mesh ) ; 
 	 	Var->ChemistryUpdate( mesh, Config ) ; 
	 	post->OutputFlow( mesh, Config, Var, 0, 0 ) ;

	 	
		#if ( Debug == true ) 
		PetscPrintf( PETSC_COMM_WORLD, "Init. potential & electric field\n" ) ; 
		#endif
		// if ( poisson_test == true ) {
		// 	poisson_solver->SOLVE_TEST( Config, Var ) ;
		// } else {
		// 	poisson_solver->SOLVE( Config, Var ) ;
		// }

		ofstream FileOutput, PCB_FileOutput , ParticleMON_FileOutput;
		
		#if (FDMaxwell == true )
		ofstream ICP_FileOutput ;
		#endif
		
		map< int, CElectrical>::iterator Iter;

		double number_temp = 0.0 ;// This just for computation	
		
		if ( mpi_rank == MASTER_NODE ){
			PCB_FileOutput.open( "PCB.dat" , ios::out | ios::trunc ) ;
			PCB_FileOutput<<"VARIABLES=\"Time\", \"PowerAbs [W]\", \"Bias Voltage [V]\" , \"Residue\""<<endl ;
		
		#if (FDMaxwell == true )
			ICP_FileOutput.open( "ICP.dat" , ios::out | ios::trunc ) ;
			ICP_FileOutput<<"VARIABLES=\"Cycle\", \"Controlled power [W]\", \"EM_power [W]\" , \"Current [A]\" , \"ES_power [W]\", \"#Maxwell_equ_solved\" "<<endl ;
		#endif
		
			ParticleMON_FileOutput.open( "ParticleMON.dat" , ios::out | ios::trunc ) ;	
			ParticleMON_FileOutput<<"VARIABLES=\"Cycle\", " ;			
			for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ )	{	
			ParticleMON_FileOutput<<"\"Two_norm_diff["+to_string(jSpecies)+"]"+"\" " ;
			ParticleMON_FileOutput<<"\"Total_particle["+to_string(jSpecies)+"]"+"\" " ;			
			}
			ParticleMON_FileOutput<< endl ;		
		}
		
		double BiasVoltage=0.0 ;
		double Target_POWER = 100.0, Current_Voltage=0.0 ;
	//exit(1);

		#if (time_monitor == true )
		double time_monitor_start, time_monitor_end ;
		double time_cost = 0 ;
		#endif

		bool power_control = false ;

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
			
			/* TEST */
			/*--- Main step loop in eack cycle. ---*/
 			for ( int MainStep = 0 ; MainStep < Config->StepPerCycle ; MainStep++ ) {
 				
 				time_monitor_start = clock() ;

 				if( WRT_CYC_AVG and MainStep%(Config->StepPerCycle/Config->WRT_Insta_Freq) == 0  ) WRT_INS = true ;
 				else  WRT_INS 	= false ;

 				if( MON_CYC and MainStep%(Config->StepPerCycle/Config->MON_Insta_Freq) == 0  ) MON_INS = true ;
 				else  MON_INS 	= false ;

 				//cout<<"MON_CYC: "<<MON_CYC<<"\t"<<"MON_INS"<<MON_INS<<endl;
				if( mpi_rank == MASTER_NODE and MON_CYC and MON_INS ){
					cout<<"MainCycle: "<<MainCycle<<"\t"<<"MainStep: "<<MainStep<<"\t"<<"Voltage: "<<Var->Volt<<"  [V]"<<"\t"<<"PhysicalTime: "<<Var->PhysicalTime<<"  [s]"<<endl ; 
					#if (time_monitor == true )
					cout<<"Time cost per time step now is: "<< time_cost/CLOCKS_PER_SEC << " s" << endl  ;
					#endif
					cout<<"Poissn ksp iter : "<<poisson_solver->its<<endl ;
					for ( int iEqn = 0 ; iEqn < DriftDiffusionNum ; iEqn++ ) {
						cout<<"Continuity["<<iEqn<<"] ksp iter : "<<continuity_solver[ iEqn ]->its<<endl ;
						
					}
					#if (FDMaxwell == true )
					cout<<"FD_maxwell ksp iter : "<<FD_maxwell_solver->its<<endl ;
					cout<<"Power absorbed by electron from inductive field is " << Var->power_inductive << " W." << endl ; 
					cout<<"Power absorbed by electron from static field is " << Var->power_static << " W." << endl ; 
					cout<<"Coil current now is " << Var->Coil_Current << " A." << endl ; 
					cout<<"I have solved Maxwell " << Var->Maxwell_solver_count << " times" << endl ;
					#endif
					if(Config->PFM_Assumption == "LMEA") 
					cout<<"Energy ksp iter : "<<electron_energy_solver->its<<endl<<endl;
					
					number_temp = double(MainCycle) + double(MainStep)/double(Config->StepPerCycle) ;							
					ParticleMON_FileOutput<<number_temp<<"\t";

				} 
				
				//Monitor the change rate of all variable
				//if( MON_CYC and MON_INS ){	
				if( false ){ //close by KL	

					for ( int jSpecies = 0 ; jSpecies < Config->TotalSpeciesNum ; jSpecies++ ) {	

						Var->two_norm_diff = 0 ;
						Var->total_particle = 0 ;
						for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
							
							Cell *cell = plasma.get_cell(i) ;				

							if(Var->U0[ jSpecies ][ i ] > 0.1*Config->Species[jSpecies].InitialDensity) 							
							Var->two_norm_diff	+= ( 1 - Var->U0[ jSpecies ][ i ]/Var->PreU0[ jSpecies ][ i ] )*( 1 - Var->U0[ jSpecies ][ i ]/Var->PreU0[ jSpecies ][ i ] );

							Var->total_particle +=   Var->U0[ jSpecies ][ i ] * cell->volume;
						}

						Var->two_norm_diff =	 plasma.parallel_sum( &Var->two_norm_diff ) ;	
						Var->two_norm_diff = 	sqrt (Var->two_norm_diff)/plasma.Mesh.cell_number * Config->StepPerCycle ;
						
						Var->total_particle =	 plasma.parallel_sum( &Var->total_particle ) ;	
					
						if( mpi_rank == MASTER_NODE ){					
						ParticleMON_FileOutput<<Var->two_norm_diff<<"\t";
						ParticleMON_FileOutput<<Var->total_particle<<"\t";	
						}						
					} 	
					if( mpi_rank == MASTER_NODE ) {
							
						#if (FDMaxwell == true )

						ICP_FileOutput<<number_temp<<"\t"<<Var->Controlled_Coil_power<<"\t"<<Var->power_inductive<<"\t"<<Var->Coil_Current<<"\t"<<Var->power_static<<"\t"<<Var->Maxwell_solver_count<<endl;
						#endif	
						
						ParticleMON_FileOutput<<endl ;	
						
					}			
				} 				
			

 				/* Update solution: Copy solution (n+1) step -> (n) step */ 
 				Var->UpdateSolution( mesh ) ; 
 				
	 				#if ( Debug == true )
	 					PetscPrintf( PETSC_COMM_WORLD, "[%d] UpdateSolution done...\n", MainStep ) ;
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
				if ( poisson_test == true ) {
					poisson_solver->SOLVE_TEST( Config, Var ) ;
				} else {
					poisson_solver->SOLVE( Config, Var ) ;
				}
	 				#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, "poisson_solver done...\n" ) ;
	 				#endif
	 				
				/* Solve number density using D-D approximation for next time step (n+1). */
				for ( int iEqn = 0 ; iEqn < DriftDiffusionNum ; iEqn++ ) 
				{
					continuity_solver[ iEqn ]->Solve( mesh, Config, Var ) ;

					#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, "Continuity (DD) [%d] done.\n",iEqn ) ;
	 				#endif
				}

				/*--- Solve number density and flux for next time step (k+1). ---*/
				for ( int iEqn = 0 ; iEqn < FullEqnNum ; iEqn++ ) {
					if ( fluid_model_solver[ iEqn ]->SpeciesType != ION ) {
						fluid_model_solver[ iEqn ]->Solve_Continuity( mesh, Config, Var ) ;
					}

					fluid_model_solver[ iEqn ]->Solve_Momentum( mesh, Config, Var ) ;

	 				#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, "Full ion equations[%d] done.\n",iEqn ) ;
	 				#endif

				}
				
				/* Solve for maxwell. */
 				#if ( FDMaxwell == true )
 					if (Var->Controlled_Coil_power < Var->Coil_power) 
					Var->Controlled_Coil_power = Var->Controlled_Coil_power + Var->power_grows_rate / Config->StepPerCycle ;
 					
 					if ( abs(Var->power_inductive - Var->Controlled_Coil_power)/Var->Controlled_Coil_power > 0.001){ 
					 	FD_maxwell_solver->SOLVE( Config, Var ) ; 
					
						if( mpi_rank == MASTER_NODE ) {
						number_temp = double(MainCycle) + double(MainStep)/double(Config->StepPerCycle) ;							
						ICP_FileOutput<<number_temp<<"\t"<<Var->Controlled_Coil_power<<"\t"<<Var->power_inductive<<"\t"<<Var->Coil_Current<<"\t"<<Var->power_static<<"\t"<<Var->Maxwell_solver_count<<endl;				
						} 
					 
					}else if( MainStep%100 == 0 ){ FD_maxwell_solver->SOLVE( Config, Var ) ; 
						if( mpi_rank == MASTER_NODE ) {
						number_temp = double(MainCycle) + double(MainStep)/double(Config->StepPerCycle) ;							
						ICP_FileOutput<<number_temp<<"\t"<<Var->Controlled_Coil_power<<"\t"<<Var->power_inductive<<"\t"<<Var->Coil_Current<<"\t"<<Var->power_static<<"\t"<<Var->Maxwell_solver_count<<endl;				
						} 	
					}

 					FD_maxwell_solver->UltraMPPComputePowerAbsorptionFromMaxwell( Config, Var ) ;	
 					FD_maxwell_solver->UltraMPPComputeTotalPower( Config, Var ) ;
	 				#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, "FD_maxwell_solver done...\n" ) ;
	 				#endif 
	 			#endif

				Var->CalculateElectrodeCurrent( mesh, Config ) ;
				#if (Debug == true ) 
	 				PetscPrintf( PETSC_COMM_WORLD, "CalculateElectrodeCurrent done...\n" ) ;
	 			post->OutputFlow( mesh, Config, Var, 0, 0 ) ;	
	 			#endif 
				/* Solve energy density for next time step (n+1). */
				if ( Config->PFM_Assumption == "LMEA" and DriftDiffusionNum > 0 ) {
					electron_energy_solver->Solver( mesh, Config, Var ) ;
					#if (Debug == true ) 
	 					PetscPrintf( PETSC_COMM_WORLD, "Energy density (DD) done.\n" ) ;
	 				#endif
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
				if( WRT_CYC_AVG ) {
					Var->AddAverage( mesh, Config ) ;
				}
				Var->AddAverage_PowerAbs( mesh, Config ) ;
				Var->AddAverage_Electrode( mesh, Config ) ;
 				Var->PhysicalTime += Var->Dt ;

 				time_monitor_end = clock() ;
 				time_cost = time_monitor_end - time_monitor_start ;
			}//End Main Step


			if ( MainCycle % 10 == 0 and power_control ){
				/* Self-Bias */
				for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
    				if ( Iter->first == POWER and fabs(Var->I_AvgPowerElectrode) > 5.E-5 ) {
	    				//Iter->second.BiasVoltage += (Var->I_AvgPowerElectrode*Var->Dt*Config->StepPerCycle)/(500*1.0E-12)*0.5 ;
    					//BiasVoltage += (Var->I_AvgPowerElectrode*Var->Dt*Config->StepPerCycle)/(500*1.0E-12)*0.5 ;
    				}
				}
				/* Power Control */
				if (       Var->AvgPowerAbs < (Target_POWER+0.05*Target_POWER) ){
					for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
    				if ( Iter->first == POWER ) {
    					Current_Voltage = Iter->second.Voltage_p2p ;
	    				Iter->second.Voltage_p2p = Current_Voltage + 0.5*fabs(Target_POWER - Var->AvgPowerAbs );
    				}
					}
				} else if( Var->AvgPowerAbs < (Target_POWER-0.05*Target_POWER) ){
					for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
    				if ( Iter->first == POWER ) {
    					Current_Voltage = Iter->second.Voltage_p2p ;
	    				Iter->second.Voltage_p2p = Current_Voltage - 0.5*fabs(Target_POWER - Var->AvgPowerAbs );
    				}
					}
				} else {
					for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
    				if ( Iter->first == POWER ) {
    					Current_Voltage = Iter->second.Voltage_p2p ;
    				}
					}
				}
			}


			if( mpi_rank == MASTER_NODE and MON_CYC ){
				cout<<"BiasVoltage: "<< BiasVoltage <<endl ;
				cout<<"Power Current[A] : "<< Var->I_AvgPowerElectrode <<endl ;
				cout<<"POWER [W]: "<<Var->AvgPowerAbs <<endl;
				for ( Iter=Config->ElectricalMap.begin() ; Iter!=Config->ElectricalMap.end() ; ++Iter ) {
    			if ( Iter->first == POWER ) {
						cout<<"Voltage [V]:"<<Iter->second.Voltage_p2p<<endl;
    			}
				}
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






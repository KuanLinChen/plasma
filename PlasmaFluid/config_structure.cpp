#include "config_structure.hpp"

CConfig::CConfig()
{
}
void CConfig::Init( string Path )
{
	MaxFrequency 		= 0.0 ;
	StepPerCycle 		= 400 ;
	SimulationCycles 	= 100 ;
	DC_TimeStep 		= -999.999 ;
	BGMass_Kg = 0.0 ;
	BGTemperature_K = 0.0 ;
	Cycle 				= 0 ;
	CasePath = Path ;
	eMeanEnergyFile = "NULL";
	eEnergyLossFile = "";
	SecondaryElectronEmissionEnergy = 0.0 ;
	SecondaryElectronEmissionCoeff 	= 0.0 ;
	ReflectionCoefficient 			= 0.0 ;
	ReadSolverFile (  CasePath+"Solver.inp" ) ;
	ReadSpeciesFile(  CasePath+"Species.inp" ) ;
}
void CConfig::ReadSpeciesFile( string FileName )
{
	// Qe = 1.602e-19 ;
 	// PI = 4.0*atan(1.0) ;
	// Me = 9.1095E-31 ;//[ kg ]
	// Kb = 1.38065E-23 ;//[ m2 kg s-2 K-1]
	// K2eV = Qe/Kb ;
	double K2eV = 1.602e-19/1.38065E-23 ;
	json config ;
	ifstream jsonFile ;

	jsonFile.open( FileName.c_str(), ifstream::in ) ; 

	if( !jsonFile ){
		cout<<"Please checked the \"Species.inp\" file."<<endl ;
		exit(1) ;
	} else {
		jsonFile >> config ;
	}

	CSpecies species_initial ;
	species_initial.Type                      = -99999 ;
	species_initial.Name                      = "NULL" ;
	species_initial.Index                     = -99999 ;
	species_initial.MobilityType              = -99999 ;
	species_initial.DiffusivityType           = -99999 ;
	species_initial.Charge                    = -99999 ;
	species_initial.Gamma                     = 1.666666667 ;
	species_initial.Mass_Kg                   = -9999. ;
	species_initial.InitialDensity            = -9999. ;
	species_initial.InitialPressure           = -9999. ;
	species_initial.InitialTemperature        = -9999. ;
	species_initial.MobilityValue             = -9999. ;
	species_initial.DiffusivityValue          = -9999. ;
	species_initial.MobilityFile              = "NULL" ;
	species_initial.DiffusivityFile           = "NULL" ;
	species_initial.CollisionFreqFile         = "NULL" ;
	species_initial.ConstantMobilityUpdate    = false ;
	species_initial.ConstantDiffusivityUpdate = false ;
	species_initial.LJ_Potential              = false ;
	species_initial.LJ_Parameter              = false ;
	species_initial.Activate                  = false ;

	if ( mpi_rank == 0 ) cout<<"/*-----  Electron group informations -----*/"<<endl ;

	/*--- Electron information ---*/
		json ElectronGroup, SubElectron ;
		if ( config.find ("ElectronGroup") != config.end() ) 
		{
			ElectronGroup = config[ "ElectronGroup" ] ; 

			ElelectronNum = 1 ;

			Species.push_back( species_initial ) ;

			SubElectron = ElectronGroup["SubGroup"][ 0 ] ;

			Species[ 0 ].Type  = ELECTRON ;
			Species[ 0 ].Index = 0 ;

			if ( SubElectron.find ("Name") 				!= SubElectron.end() ) Species[ 0 ].Name 				=  SubElectron[ "Name" ] ;

			if ( SubElectron.find ("MobilityType"	) 	!= SubElectron.end() ) Species[ 0 ].MobilityType 		=  SubElectron[ "MobilityType" 	] ;
			if ( SubElectron.find ("MobilityValue"	) 	!= SubElectron.end() ) Species[ 0 ].MobilityValue 		=  SubElectron[ "MobilityValue" ] ;
			if ( SubElectron.find ("MobilityFile"	) 	!= SubElectron.end() ) Species[ 0 ].MobilityFile 		=  SubElectron[ "MobilityFile" 	] ;


			if ( SubElectron.find ("DiffusivityType"	) 	!= SubElectron.end() ) Species[ 0 ].DiffusivityType 	=  SubElectron[ "DiffusivityType" 	] ;
			if ( SubElectron.find ("DiffusivityFile"	) 	!= SubElectron.end() ) Species[ 0 ].DiffusivityFile 	=  SubElectron[ "DiffusivityFile" 	] ;
			if ( SubElectron.find ("DiffusivityValue"	) 	!= SubElectron.end() ) Species[ 0 ].DiffusivityValue 	=  SubElectron[ "DiffusivityValue" 	] ;
			
			if ( SubElectron.find ("CollisionFreqFile") 	!= SubElectron.end() ) Species[ 0 ].CollisionFreqFile 	=  SubElectron[ "CollisionFreqFile" ] ;

			if ( SubElectron.find ("Charge") 			!= SubElectron.end() ) {
				Species[ 0 ].Charge 				=  SubElectron[ "Charge" ] ;
				Species[ 0 ].sgn 					=  -1.0 ;
			}
			if ( SubElectron.find ("Gamma") 			!= SubElectron.end() ) Species[ 0 ].Gamma 				=  SubElectron[ "Gamma"  ] ;
			if ( SubElectron.find ("Mass_Kg") 			!= SubElectron.end() ) Species[ 0 ].Mass_Kg 			=  SubElectron[ "Mass_Kg" ] ;
			if ( SubElectron.find ("InitialDensity") 	!= SubElectron.end() ) Species[ 0 ].InitialDensity 		=  SubElectron[ "InitialDensity" ] ;
			if ( SubElectron.find ("InitialTemperature")!= SubElectron.end() ) Species[ 0 ].InitialTemperature 	=  SubElectron[ "InitialTemperature" ] ;
			if ( SubElectron.find ("InitialTemperature_K")!= SubElectron.end() ){
					double tmp_Temperature_K  =  SubElectron[ 0 ][ "InitialTemperature_K" ] ;
					Species[ 0 ].InitialTemperature =  tmp_Temperature_K/K2eV ;//To eV
			} 

			if ( SubElectron.find ("Activate")!= SubElectron.end() ) {
				if ( SubElectron[ "Activate" ] == 1 )  Species[ 0 ].Activate = true ;
			} else {
				if ( mpi_rank == 0 ) cout<<"U should specify Activate or not "<<endl ;
			}
			/*--- For debug ---*/
			if ( mpi_rank == 0 ){
				// cout<<"Name              : "<<Species[ 0 ].Name<<endl;
				// cout<<"Index             : "<<Species[ 0 ].Index<<endl;
				// cout<<"Charge            : "<<Species[ 0 ].Charge<<endl;
				// cout<<"Gamma             : "<<Species[ 0 ].Gamma<<endl;
				// cout<<"Mass_Kg           : "<<Species[ 0 ].Mass_Kg<<endl;
				// cout<<"MobilityType      : "<<Species[ 0 ].MobilityType<<endl;
				// cout<<"MobilityFile      : "<<Species[ 0 ].MobilityFile<<endl;
				// cout<<"DiffusivityType   : "<<Species[ 0 ].DiffusivityType<<endl;
				// cout<<"DiffusivityFile   : "<<Species[ 0 ].DiffusivityFile<<endl;
				// cout<<"InitialDensity    : "<<Species[ 0 ].InitialDensity<<endl;
				// cout<<"InitialTemperature: "<<Species[ 0 ].InitialTemperature<<endl<<endl;
				if ( Species[ 0 ].MobilityType == 0 and Species[ 0 ].MobilityValue < 0.0 ) {
					cout<<"The electron Mobility is choose be constant, please assign the value [ e.g. \"MobilityValue\": 30.0 ]"<<endl;
					exit(1) ;
				} else if( Species[ 0 ].DiffusivityType == 0 and Species[ 0 ].DiffusivityValue < 0.0 ){
					cout<<"The electron Diffusivity is choose be constant, please assign the value [ e.g. \"DiffusivityValue\": 120.0 ]"<<endl;
					exit(1) ;
				}
			}

		} else {

			cout<<"Can't not find \"Electron\" Group in "<<FileName<<endl ; 
			exit(1) ; 

		} 


	/*--- Ion species information ---*/
		if ( mpi_rank == 0 ) cout<<"/*----- Ion group informations -----*/"<<endl;
		json IonGroup, *SubIon ;
		int ii = ElelectronNum ;
		if ( config.find ("IonGroup") != config.end() ) { 

			IonGroup = config[ "IonGroup" ] ; 

			IonNum = IonGroup["SubGroup"].size() ; 
			if( mpi_rank == 0 ) cout<<"Number of ion Specie: "<<IonNum<<endl;

			SubIon = new json [ IonNum ] ;

			for( int i = 0 ; i < IonNum ; i++ ) Species.push_back( species_initial ) ;

			for( int i = 0 ; i < IonNum ; i++ ) SubIon[ i ] = IonGroup["SubGroup"][ i ] ;

			for( int i = 0 ; i < IonNum ; i++ ){
				
				Species[ i+ii ].Type  =  ION ;
				Species[ i+ii ].Index =  i+ii ;

				if ( SubIon[ i ].find ("Name") 				!= SubIon[ i ].end() ) Species[ i+ii ].Name 			 =  SubIon[ i ][ "Name" ] ;

				if ( SubIon[ i ].find ("MobilityType"	) 	!= SubIon[ i ].end() ) Species[ i+ii ].MobilityType 	 =  SubIon[ i ][ "MobilityType" ] ;
				if ( SubIon[ i ].find ("MobilityFile"	) 	!= SubIon[ i ].end() ) Species[ i+ii ].MobilityFile 	 =  SubIon[ i ][ "MobilityFile" ] ;
				if ( SubIon[ i ].find ("MobilityValue"	) 	!= SubIon[ i ].end() ) Species[ i+ii ].MobilityValue 	 =	SubIon[ i ][ "MobilityValue"] ;

				if ( SubIon[ i ].find ("DiffusivityType"	) 	!= SubIon[ i ].end() ) Species[ i+ii ].DiffusivityType 	 =  SubIon[ i ][ "DiffusivityType" ] ;
				if ( SubIon[ i ].find ("DiffusivityFile"	) 	!= SubIon[ i ].end() ) Species[ i+ii ].DiffusivityFile 	 =  SubIon[ i ][ "DiffusivityFile" ] ;
				if ( SubIon[ i ].find ("DiffusivityValue"	) 	!= SubIon[ i ].end() ) Species[ i+ii ].DiffusivityValue  =  SubIon[ i ][ "DiffusivityValue"] ;

				if ( SubIon[ i ].find ("Charge") 			!= SubIon[ i ].end() ) {
					Species[ i+ii ].Charge 			 =  SubIon[ i ][ "Charge" ] ;
					if ( Species[ i+ii ].Charge < 0 ) 	Species[ i+ii ].sgn = -1.0 ;
					else 								Species[ i+ii ].sgn =  1.0 ;
				}
				if ( SubIon[ i ].find ("Gamma") 			!= SubIon[ i ].end() ) Species[ i+ii ].Gamma 			 =  SubIon[ i ][ "Gamma" ] ;
				if ( SubIon[ i ].find ("Mass_Kg") 			!= SubIon[ i ].end() ) Species[ i+ii ].Mass_Kg 			 =  SubIon[ i ][ "Mass_Kg" ] ;
				if ( SubIon[ i ].find ("Mass_AMU") 			!= SubIon[ i ].end() ){
					double tmp_mass = 0.0 ;
					tmp_mass			 	=  SubIon[ i ][ "Mass_AMU" ] ;
					Species[ i+ii ].Mass_Kg =  tmp_mass*( 1.66058e-27 ) ;
				}
				if ( SubIon[ i ].find ("InitialDensity") 	!= SubIon[ i ].end() ) Species[ i+ii ].InitialDensity 	  = SubIon[ i ][ "InitialDensity" ] ;
				if ( SubIon[ i ].find ("InitialTemperature")!= SubIon[ i ].end() ) Species[ i+ii ].InitialTemperature = SubIon[ i ][ "InitialTemperature" ] ;
				if ( SubIon[ i ].find ("InitialTemperature_K")!= SubIon[ i ].end() ){
					double tmp_Temperature_K 	 	 			 =  SubIon[ i ][ "InitialTemperature_K" ] ;
					Species[ i+ii ].InitialTemperature 	 =  tmp_Temperature_K/K2eV ;//To eV
				} 

				if ( SubIon[ i ].find ("Activate")!= SubIon[ i ].end() ) {
					if ( SubIon[ i ][ "Activate" ] == 1 )  Species[ i+ii ].Activate = true ;
				} else {
					if ( mpi_rank == 0 ) cout<<"U should specify Activate or not for ion species "<<endl ;
				}	


				if ( mpi_rank == 0 ){
					// cout<<"Name              : "<<Species[ i+ii ].Name<<endl;
					// cout<<"Index             : "<<Species[ i+ii ].Index<<endl;
					// cout<<"Charge            : "<<Species[ i+ii ].Charge<<endl;
					// cout<<"Gamma             : "<<Species[ i+ii ].Gamma<<endl;
					// cout<<"Mass_Kg           : "<<Species[ i+ii ].Mass_Kg<<endl;
					// cout<<"MobilityType      : "<<Species[ i+ii ].MobilityType<<endl;
					// cout<<"MobilityFile      : "<<Species[ i+ii ].MobilityFile<<endl;
					// cout<<"DiffusivityType   : "<<Species[ i+ii ].DiffusivityType<<endl;
					// cout<<"DiffusivityFile   : "<<Species[ i+ii ].DiffusivityFile<<endl;
					// cout<<"InitialDensity    : "<<Species[ i+ii ].InitialDensity<<endl;
					// cout<<"InitialTemperature: "<<Species[ i+ii ].InitialTemperature<<endl<<endl;
					if ( Species[ i+ii ].MobilityType == 0 and Species[ i+ii ].MobilityValue < 0.0 ) {
						cout<<"The "<<Species[ i+ii ].Name<<" mobility is choose be constant, please assign the value [ e.g. \"MobilityValue\": 0.14 ]"<<endl;
						exit(1) ;
					} else if( Species[ i+ii ].DiffusivityType == 0 and Species[ i+ii ].DiffusivityValue < 0.0 ){
						cout<<"The "<<Species[ i+ii ].Name<<"Diffusivity is choose be constant, please assign the value [ e.g. \"DiffusivityValue\": 4.0E-3 ]"<<endl;
						exit(1) ;
					}
				}

			}
		} else {

			cout<<"Can't not find \"Ion\" Group in "<<FileName<<endl ; 
			exit(1) ; 

		} 
		ChargeSpeciesNum = ElelectronNum + IonNum ;

	/*--- Nertral species information ---*/
		if ( mpi_rank == 0 ) cout<<"/*----- Nertral group informations -----*/"<<endl;
		json NeutralGroup, *SubNeutral ;	
		ii = ElelectronNum+IonNum ;

		if ( config.find ("NeutralGroup") != config.end() ) { 

			NeutralGroup = config[ "NeutralGroup" ] ; 

			NeutralNum = NeutralGroup["SubGroup"].size() ; if( mpi_rank == 0 ) cout<<"Neutral Specie Number: "<<NeutralNum<<endl;

			SubNeutral = new json [ NeutralNum ] ;

			for( int i = 0 ; i < NeutralNum ; i++ ) Species.push_back( species_initial ) ;

			for( int i = 0 ; i < NeutralNum ; i++ ) SubNeutral[ i ] = NeutralGroup["SubGroup"][ i ] ;

			for( int i = 0 ; i < NeutralNum ; i++ ){

				Species[ i+ii ].Type  =  NEUTRAL ;
				Species[ i+ii ].Index =  i+ii ;

				if ( SubNeutral[ i ].find ("Name") != SubNeutral[ i ].end() ) Species[ i+ii ].Name =  SubNeutral[ i ][ "Name" ] ;

				Species[ i+ii ].MobilityType  =  0 ;
				Species[ i+ii ].MobilityFile  =  "N/A" ;
				Species[ i+ii ].MobilityValue =  0.0 ;

				if ( SubNeutral[ i ].find ("DiffusivityType") 	!= SubNeutral[ i ].end() ) Species[ i+ii ].DiffusivityType 	 =  SubNeutral[ i ][ "DiffusivityType" ] ;
				if ( SubNeutral[ i ].find ("DiffusivityValue") 	!= SubNeutral[ i ].end() ) Species[ i+ii ].DiffusivityValue  =  SubNeutral[ i ][ "DiffusivityValue" ] ;
				if ( SubNeutral[ i ].find ("DiffusivityFile") 	!= SubNeutral[ i ].end() ) Species[ i+ii ].DiffusivityFile 	 =  SubNeutral[ i ][ "DiffusivityFile" ] ;

				if ( SubNeutral[ i ].find ("Charge") 			!= SubNeutral[ i ].end() ) {
					Species[ i+ii ].Charge 			 =  SubNeutral[ i ][ "Charge" ] ;
					if ( Species[ i+ii ].Charge == 0.0 ) {
						 Species[ i+ii ].sgn = 0.0 ;
					}else{
						 cout<<"Error in neutral charge"<<endl;exit(1);
					}
				}
				if ( SubNeutral[ i ].find ("Gamma") 			!= SubNeutral[ i ].end() ) Species[ i+ii ].Gamma 			 =  SubNeutral[ i ][ "Gamma" ] ;
				if ( SubNeutral[ i ].find ("Mass_Kg") 			!= SubNeutral[ i ].end() ) Species[ i+ii ].Mass_Kg 			 =  SubNeutral[ i ][ "Mass_Kg" ] ;
				if ( SubNeutral[ i ].find ("Mass_AMU") 			!= SubNeutral[ i ].end() ){
					double tmp_mass = 0.0 ;
					tmp_mass			 	=  SubNeutral[ i ][ "Mass_AMU" ] ;
					Species[ i+ii ].Mass_Kg =  tmp_mass*( 1.66058e-27 ) ;
				}
				if ( SubNeutral[ i ].find ("InitialDensity") 	!= SubNeutral[ i ].end() ) Species[ i+ii ].InitialDensity 	  = SubNeutral[ i ][ "InitialDensity" ] ;
				if ( SubNeutral[ i ].find ("InitialTemperature")!= SubNeutral[ i ].end() ) Species[ i+ii ].InitialTemperature = SubNeutral[ i ][ "InitialTemperature" ] ;

				if ( SubNeutral[ i ].find ("InitialTemperature_K")!= SubNeutral[ i ].end() ){
					double tmp_Temperature_K 	 	 			 =  SubNeutral[ i ][ "InitialTemperature_K" ] ;
					Species[ i+ii ].InitialTemperature 	 =  tmp_Temperature_K/K2eV ;//To eV
				} 

				if ( SubNeutral[ i ].find ("Activate")!= SubNeutral[ i ].end() ) {
					if ( SubNeutral[ i ][ "Activate" ] == 1 )  Species[ i+ii ].Activate = true ;
				} else {
					if ( mpi_rank == 0 ) cout<<"U should specify Activate or not for ion species "<<endl ;
				}	

				if ( mpi_rank == 0 ){
					// cout<<"Name              : "<<Species[ i+ii ].Name<<endl;
					// cout<<"Index             : "<<Species[ i+ii ].Index<<endl;
					// cout<<"Charge            : "<<Species[ i+ii ].Charge<<endl;
					// cout<<"Gamma             : "<<Species[ i+ii ].Gamma<<endl;
					// cout<<"Mass_Kg           : "<<Species[ i+ii ].Mass_Kg<<endl;
					// cout<<"MobilityType      : "<<Species[ i+ii ].MobilityType<<endl;
					// cout<<"MobilityFile      : "<<Species[ i+ii ].MobilityFile<<endl;
					// cout<<"DiffusivityType   : "<<Species[ i+ii ].DiffusivityType<<endl;
					// cout<<"DiffusivityFile   : "<<Species[ i+ii ].DiffusivityFile<<endl;
					// cout<<"InitialDensity    : "<<Species[ i+ii ].InitialDensity<<endl;
					// cout<<"InitialTemperature: "<<Species[ i+ii ].InitialTemperature<<endl<<endl;
					if ( Species[ i+ii ].MobilityType == 0 and Species[ i+ii ].MobilityValue < 0.0 ) {
						cout<<"The "<<Species[ i+ii ].Name<<" mobility should be zero, please assign the value [ e.g. \"MobilityValue\": 0.0 ]"<<endl;
						exit(1) ;
					} else if( Species[ i+ii ].DiffusivityType == 0 and Species[ i+ii ].DiffusivityValue < 0.0 ){
						cout<<"The "<<Species[ i+ii ].Name<<" Diffusivity is constant, please assign the value [ e.g. \"DiffusivityValue\": 4.0E-3 ]"<<endl;
						exit(1) ;
					}
				}
			}
		} else { cout<<"Can't not find \"Neutral\" Group in "<<FileName<<endl ; exit(1) ; }


	/*--- Background species information ---*/
		if ( mpi_rank == 0 ) cout<<"/*----- Background group informations -----*/"<<endl ;
		json BackGroundGroup, *SubBackGround ;
		ii = ElelectronNum+IonNum+NeutralNum ;
		if ( config.find ("BackGroundGroup") != config.end() ){ 

			BackGroundGroup = config[ "BackGroundGroup" ] ; 

			BackGroundNum = BackGroundGroup["SubGroup"].size() ; 

			if( mpi_rank == 0 ) cout<<"BackGround Specie Number: "<<BackGroundNum<<endl;

			SubBackGround = new json [ BackGroundNum ] ;

			for( int i = 0 ; i < BackGroundNum ; i++ ) Species.push_back( species_initial ) ;

			for( int i = 0 ; i < BackGroundNum ; i++ ) SubBackGround[ i ] = BackGroundGroup["SubGroup"][ i ] ;


			P_back = 0.0 ;
			T_back = 0.0 ;
			for( int i = 0 ; i < BackGroundNum ; i++ ){

				Species[ i+ii ].Type  =  BACKGROUND ;
				Species[ i+ii ].Index =  i+ii ;

				if ( SubBackGround[ i ].find ("Name") 				 != SubBackGround[ i ].end() ) Species[ i+ii ].Name 			 	 =  SubBackGround[ i ][ "Name" ] ;

				Species[ i+ii ].MobilityType  =  0 ;
				Species[ i+ii ].MobilityFile  =  "N/A" ;
				Species[ i+ii ].MobilityValue =  0.0 ;
				
				if ( SubBackGround[ i ].find ("DiffusivityType") 	!= SubBackGround[ i ].end() ) Species[ i+ii ].DiffusivityType 	 =  SubBackGround[ i ][ "DiffusivityType" ] ;
				if ( SubBackGround[ i ].find ("DiffusivityValue") 	!= SubBackGround[ i ].end() ) Species[ i+ii ].DiffusivityValue   =  SubBackGround[ i ][ "DiffusivityValue" ] ;
				if ( SubBackGround[ i ].find ("DiffusivityFile") 	!= SubBackGround[ i ].end() ) Species[ i+ii ].DiffusivityFile 	 =  SubBackGround[ i ][ "DiffusivityFile" ] ;

				if ( SubBackGround[ i ].find ("InitialPressure") 	 != SubBackGround[ i ].end() ) Species[ i+ii ].InitialPressure 	 	 =  SubBackGround[ i ][ "InitialPressure" ] ;
				if ( SubBackGround[ i ].find ("InitialTemperature")  != SubBackGround[ i ].end() ) Species[ i+ii ].InitialTemperature	 =  SubBackGround[ i ][ "InitialTemperature" ] ;

				if ( SubBackGround[ i ].find ("InitialTemperature_K")!= SubBackGround[ i ].end() ){
					double tmp_Temperature_K 	 	 			 =  SubBackGround[ i ][ "InitialTemperature_K" ] ;
					T_back = tmp_Temperature_K ;
					Species[ i+ii ].InitialTemperature 	 =  tmp_Temperature_K/K2eV ;//To eV
					double Pressure = Species[ i+ii ].InitialPressure ;
					P_back += Pressure ;
					double BackGroundNumberDensity = 133.2894 * Pressure / 1.38065E-23 / tmp_Temperature_K ;
					Species[ i+ii ].InitialDensity = BackGroundNumberDensity ;
					if( tmp_Temperature_K > BGTemperature_K ) BGTemperature_K = tmp_Temperature_K ;
				} 
					
					if ( SubBackGround[ i ].find ("Charge")  != SubBackGround[ i ].end() ) Species[ i+ii ].Charge = SubBackGround[ i ][ "Charge" ] ;
					if ( SubBackGround[ i ].find ("Gamma")   != SubBackGround[ i ].end() ) Species[ i+ii ].Gamma  = SubBackGround[ i ][ "Gamma"  ] ;
					if ( SubBackGround[ i ].find ("Mass_Kg") != SubBackGround[ i ].end() ) Species[ i+ii ].Mass_Kg= SubBackGround[ i ][ "Mass_Kg"] ;
					if ( SubBackGround[ i ].find ("Mass_AMU")!= SubBackGround[ i ].end() ) {
					double tmp_mass                                     = 0.0 ;
					tmp_mass                                            =  SubBackGround[ i ][ "Mass_AMU" ] ;
					Species[ i+ii ].Mass_Kg                             =  tmp_mass*( 1.66058e-27 ) ;
					if( Species[ i+ii ].Mass_Kg > BGMass_Kg ) BGMass_Kg = Species[ i+ii ].Mass_Kg ;
					}
					if ( SubBackGround[ i ].find ("Activate")           != SubBackGround[ i ].end() ) {
					if ( SubBackGround[ i ][ "Activate" ]               == 1 )  Species[ i+ii ].Activate = true ;
					} else {
					if ( mpi_rank                                       == 0 ) cout<<"U should specify Activate or not for ion species "<<endl ;
					}	



				if ( SubBackGround[ i ].find ("InitialDensity")   != SubBackGround[ i ].end() ) Species[ i+ii ].InitialDensity 	 =  SubBackGround[ i ][ "InitialDensity" ] ;
				if ( mpi_rank == 0 ){
					cout<<"Name              : "<<Species[ i+ii ].Name<<endl;
					cout<<"Index             : "<<Species[ i+ii ].Index<<endl;
					cout<<"Charge            : "<<Species[ i+ii ].Charge<<endl;
					cout<<"Gamma             : "<<Species[ i+ii ].Gamma<<endl;
					cout<<"Mass_Kg           : "<<Species[ i+ii ].Mass_Kg<<endl;
					cout<<"MobilityType      : "<<Species[ i+ii ].MobilityType<<endl;
					cout<<"MobilityFile      : "<<Species[ i+ii ].MobilityFile<<endl;
					cout<<"DiffusivityType   : "<<Species[ i+ii ].DiffusivityType<<endl;
					cout<<"DiffusivityFile   : "<<Species[ i+ii ].DiffusivityFile<<endl;
					cout<<"InitialDensity    : "<<Species[ i+ii ].InitialDensity<<endl;
					cout<<"InitialPressure   : "<<Species[ i+ii ].InitialPressure<<endl;
					cout<<"InitialTemperature: "<<Species[ i+ii ].InitialTemperature<<endl<<endl;
					if ( Species[ i+ii ].MobilityType == 0 and Species[ i+ii ].MobilityValue < 0.0 ) {
						cout<<"The "<<Species[ i+ii ].Name<<" mobility should be zero, please assign the value [ e.g. \"MobilityValue\": 0.0 ]"<<endl;
						exit(1) ;
					} else if( Species[ i+ii ].DiffusivityType == 0 and Species[ i+ii ].DiffusivityValue < 0.0 ){
						cout<<"The "<<Species[ i+ii ].Name<<" Diffusivity is constant, please assign the value [ e.g. \"DiffusivityValue\": 4.0E-3 ]"<<endl;
						exit(1) ;
					}
				}
			}

		} else {

			cout<<"Can't not find \"BackGround\" Group in "<<FileName<<endl ; 
			exit(1) ; 

		}
	TotalSpeciesNum = ElelectronNum + IonNum + NeutralNum + BackGroundNum ;
	//exit(1);
	if ( mpi_rank==0 ) cout<<"TotalSpeciesNum: "<<TotalSpeciesNum<<endl ;
}
void CConfig::ReadSolverFile( string FileName )
{
	json config ;
	ifstream jsonFile ;

	jsonFile.open( FileName.c_str(), ifstream::in ) ; 

	if( !jsonFile ){
		cout<<"Can not find the \"Solver.inp\" file."<<endl ;
		exit(1) ;
	} else {
		jsonFile >> config ;
	}

	/*--- Model Assumptions ---*/
	json ModelAssumptions ;
	if ( config.find ("ModelAssumptions") != config.end() ) { ModelAssumptions = config[ "ModelAssumptions" ] ; }


	string Normalize_Read="" ;
	if ( config.find ("Normalize") != config.end() ) { 

		Normalize_Read = config[ "Normalize" ] ; 
		if ( Normalize_Read == "ON" ) {
			Normalize = true ;
		} else {
			Normalize = false ;
		}
	}

	/*------*/
	if ( ModelAssumptions.find ("PFM_Assumption") != ModelAssumptions.end() ) PFM_Assumption = ModelAssumptions["PFM_Assumption"] ;

	/*------*/
	if ( ModelAssumptions.find ("PFM_SubModel"  ) != ModelAssumptions.end() ) PFM_SubModel = ModelAssumptions["PFM_SubModel"] ;

	/*------*/
	if ( ModelAssumptions.find ("StepPerCycle") != ModelAssumptions.end() ) StepPerCycle = ModelAssumptions["StepPerCycle"] ;


	/*------*/
	if ( ModelAssumptions.find ("TimeStepSize") != ModelAssumptions.end() ) DC_TimeStep = ModelAssumptions["TimeStepSize"] ;

	if ( ModelAssumptions.find ("ExitCycle") != ModelAssumptions.end() ) ExitCycle = ModelAssumptions["ExitCycle"] ;

	
	/*------*/
	if ( ModelAssumptions.find ( "SecEleEmissCoeff"	) != ModelAssumptions.end() ) SecondaryElectronEmissionCoeff  	= ModelAssumptions[ "SecEleEmissCoeff"	] ;
	if ( ModelAssumptions.find ( "SecEleEnergy"		) != ModelAssumptions.end() ) SecondaryElectronEmissionEnergy 	= ModelAssumptions[ "SecEleEnergy"		] ;
	if ( ModelAssumptions.find ( "ReflectionCoeff"	) != ModelAssumptions.end() ) ReflectionCoefficient 			= ModelAssumptions[ "ReflectionCoeff"		] ;

	/*------*/
	if ( PFM_Assumption == "LFA" ){
		if ( ModelAssumptions.find ("eMeanEnergyFile") != ModelAssumptions.end() ){
			eMeanEnergyFile = ModelAssumptions["eMeanEnergyFile"] ;
			if ( mpi_rank == MASTER_NODE ) cout<<"PFM_Assumption: "<<PFM_Assumption<<", eMeanEnergyFile: "<<eMeanEnergyFile<<endl;
		} else {
			if ( mpi_rank == MASTER_NODE ) cout<<"If choose LFA, you need to specify the reduce electric field vs Mean energy\n"<<endl; exit(1) ;
		}
	}

	if ( ModelAssumptions.find ("eEnergyLossFile") != ModelAssumptions.end() ) {

		eEnergyLossFile = ModelAssumptions["eEnergyLossFile"] ;

		if ( mpi_rank == MASTER_NODE ) cout<<"eEnergyLossFile: " << eEnergyLossFile<<endl ;
	} else {

		eEnergyLossFile = "" ;

	}


	/*--- Electrical Control ---*/
	json ElectricalControl, *SubElectricalControl ;
	if ( config.find ("ElectricalControl") != config.end() ) { ElectricalControl = config[ "ElectricalControl" ] ; }
	ElectricalComponent = ElectricalControl["SubGroup"].size() ;

	SubElectricalControl = new json [ ElectricalComponent ] ;
	for( int i = 0 ; i < ElectricalComponent ; i++ ) SubElectricalControl[ i ] = ElectricalControl["SubGroup"][ i ] ;
	CElectrical electrical_initial ;
	electrical_initial.Name 			= "NULL" ;
	electrical_initial.Type 			= -99999 ;
	electrical_initial.Frequency  		= -99999 ;
	electrical_initial.Voltage_p2p 		= -99999 ;
	electrical_initial.BiasVoltage 		= -99999 ;
	electrical_initial.VoltageFileName  = "NULL" ;

	for( int i = 0 ; i < ElectricalComponent ; i++ ) {

		Electrical.push_back( electrical_initial ) ;
		
	}

	//map< int, CElectrical> ElectricalMap ;
	//mapStudent.insert(pair<string, string>("r000", "student_zero"));	
	//CElectrical temp ;
	//ElectricalMap.insert( pair< int, CElectrical>( POWER, electrical_initial ) ) ;
	//iter = mapStudent.find("r123");
	



	for( int i = 0 ; i < ElectricalComponent ; i++ ){
		if ( SubElectricalControl[ i ].find ("Name") 			!= SubElectricalControl[ i ].end() ) Electrical[ i ].Name 			=  SubElectricalControl[ i ][ "Name" 			] ;
		if ( SubElectricalControl[ i ].find ("Voltage_p2p") 	!= SubElectricalControl[ i ].end() ) Electrical[ i ].Voltage_p2p 	=  SubElectricalControl[ i ][ "Voltage_p2p" 	] ;
		if ( SubElectricalControl[ i ].find ("Frequency") 		!= SubElectricalControl[ i ].end() ) Electrical[ i ].Frequency 	 	=  SubElectricalControl[ i ][ "Frequency" 		] ;
		if ( SubElectricalControl[ i ].find ("BiasVoltage")  	!= SubElectricalControl[ i ].end() ) Electrical[ i ].BiasVoltage	=  SubElectricalControl[ i ][ "BiasVoltage" 	] ;
		//	double RampFactor,T_min, T_max, RampTime, InitialRampFactot ;
		// if ( SubElectricalControl[ i ].find ("T_max") != SubElectricalControl[ i ].end() ) Electrical[ i ].T_max=  SubElectricalControl[ i ][ "T_max" ] ;
		// if ( SubElectricalControl[ i ].find ("T_min") != SubElectricalControl[ i ].end() ) Electrical[ i ].T_min=  SubElectricalControl[ i ][ "T_min" ] ;
		// if ( SubElectricalControl[ i ].find ("RampTime") != SubElectricalControl[ i ].end() ) Electrical[ i ].RampTime=  SubElectricalControl[ i ][ "RampTime" ] ;
		// if ( SubElectricalControl[ i ].find ("InitialRampFactot") != SubElectricalControl[ i ].end() ) Electrical[ i ].InitialRampFactot=  SubElectricalControl[ i ][ "InitialRampFactot" ] ;
		
		if ( SubElectricalControl[ i ].find ("VoltageFileName") != SubElectricalControl[ i ].end() ) Electrical[ i ].VoltageFileName=  SubElectricalControl[ i ][ "VoltageFileName" ] ;

		if 		( Electrical[ i ].Name == "POWER" 	) Electrical[ i ].Type = POWER ;
		// else if ( Electrical[ i ].Name == "POWER_1" ) Electrical[ i ].Type = POWER_1 ;
		// else if ( Electrical[ i ].Name == "POWER_2" ) Electrical[ i ].Type = POWER_2 ;
		// else if ( Electrical[ i ].Name == "POWER_3" ) Electrical[ i ].Type = POWER_3 ;
		// else if ( Electrical[ i ].Name == "POWER_4" ) Electrical[ i ].Type = POWER_4 ;
		// else if ( Electrical[ i ].Name == "POWER_5" ) Electrical[ i ].Type = POWER_5 ;
		else if ( Electrical[ i ].Name == "GROUND" 	) Electrical[ i ].Type = GROUND ;
		// else if ( Electrical[ i ].Name == "GROUND_1" ) Electrical[ i ].Type = GROUND_1 ;
		// else if ( Electrical[ i ].Name == "GROUND_2" ) Electrical[ i ].Type = GROUND_2 ;
		// else if ( Electrical[ i ].Name == "GROUND_3" ) Electrical[ i ].Type = GROUND_3 ;
		// else if ( Electrical[ i ].Name == "GROUND_4" ) Electrical[ i ].Type = GROUND_4 ;
		// else if ( Electrical[ i ].Name == "GROUND_5" ) Electrical[ i ].Type = GROUND_5 ;
		else if ( mpi_rank==0 ) cout<<"Too much power electrode pls contex K. L. Chen"<<endl;
	}

	CElectrical container ;
	map< int, CElectrical>::iterator ElectricalIter;
	for( int i = 0 ; i < ElectricalComponent ; i++ ){

		if ( SubElectricalControl[ i ].find ("Name") 			!= SubElectricalControl[ i ].end() ) container.Name 		=  SubElectricalControl[ i ][ "Name" 			] ;
		if ( SubElectricalControl[ i ].find ("Voltage_p2p") 	!= SubElectricalControl[ i ].end() ) container.Voltage_p2p 	=  SubElectricalControl[ i ][ "Voltage_p2p" 	] ;

		if ( SubElectricalControl[ i ].find ("Frequency") 		!= SubElectricalControl[ i ].end() ){
			container.Frequency 	=  SubElectricalControl[ i ][ "Frequency" ] ;

			if ( container.Frequency > MaxFrequency ) MaxFrequency = container.Frequency ; 

			if ( container.Frequency <= 0.0 ){

				container.Period = 1.E+20 ;
				
			} else {
				container.Period 		= 1.0/container.Frequency ;
			}

		} 

		if ( SubElectricalControl[ i ].find ("BiasVoltage")  	!= SubElectricalControl[ i ].end() ) container.BiasVoltage	=  SubElectricalControl[ i ][ "BiasVoltage" 	] ;
		if ( SubElectricalControl[ i ].find ("VoltageFileName") != SubElectricalControl[ i ].end() ) container.VoltageFileName=  SubElectricalControl[ i ][ "VoltageFileName" ] ;

		if 		( container.Name == "POWER"   ) container.Type = POWER ;
		// else if ( container.Name == "POWER_1" ) container.Type = POWER_1 ;
		// else if ( container.Name == "POWER_2" ) container.Type = POWER_2 ;
		// else if ( container.Name == "POWER_3" ) container.Type = POWER_3 ;
		// else if ( container.Name == "POWER_4" ) container.Type = POWER_4 ;
		// else if ( container.Name == "POWER_5" ) container.Type = POWER_5 ;
		else if ( container.Name == "GROUND"  ) container.Type = GROUND ;
		// else if ( container.Name == "GROUND_1" ) container.Type = GROUND_1 ;
		// else if ( container.Name == "GROUND_2" ) container.Type = GROUND_2 ;
		// else if ( container.Name == "GROUND_3" ) container.Type = GROUND_3 ;
		// else if ( container.Name == "GROUND_4" ) container.Type = GROUND_4 ;
		// else if ( container.Name == "GROUND_5" ) container.Type = GROUND_5 ;
		else if ( mpi_rank==0 ) cout<<"Too much power electrode pls contex K. L. Chen"<<endl;
		ElectricalMap.insert( pair< int, CElectrical>( container.Type, container ) ) ;
		ElectricalIter = ElectricalMap.find( container.Type ) ;
		//cout<<"iter: " <<ElectricalIter->second.Frequency<<endl;
	}


	json PoissonEqnControl, ElectronEqnControl, IonEqnControl, NeutralEqnControl, BackgroundEqnControl, OutputControl ;
	for ( int iEqn = 0 ; iEqn < 5 ; iEqn ++ ) Equation[ iEqn ].Init() ;

	/*--- Poisson equation control ---*/
	if ( config.find ("PoissonEqnControl") != config.end() ) 
	{ 
		PoissonEqnControl = config[ "PoissonEqnControl" ] ; 

		if ( PoissonEqnControl.find ("Equation") 	!= PoissonEqnControl.end() ) 
			Equation[ POISSON ].Equation 	= PoissonEqnControl["Equation"] ;
	}

	/*--- Electron equation control ---*/
	if ( config.find ("ElectronEqnControl") != config.end() ) 
	{ 
		ElectronEqnControl = config[ "ElectronEqnControl" ] ; 

		if ( ElectronEqnControl.find ("Equation") 	!= ElectronEqnControl.end() ) 
			Equation[ ELECTRON ].Equation 	= ElectronEqnControl["Equation"] ;

		if ( ElectronEqnControl.find ("WallBoundaryType") != ElectronEqnControl.end() ) 
			Equation[ ELECTRON ].WallBoundaryType 	= ElectronEqnControl["WallBoundaryType"] ;
	}

	/*--- Ion equation control ---*/
	if ( config.find ("IonEqnControl") != config.end() ) 
	{
		IonEqnControl = config[ "IonEqnControl" ] ; 

		if ( IonEqnControl.find ("Equation") 	!= IonEqnControl.end() ) 
			Equation[ ION ].Equation 	= IonEqnControl["Equation"] ;

		if ( IonEqnControl.find ("WallBoundaryType")  != IonEqnControl.end() ) 
			Equation[ ION ].WallBoundaryType 	= IonEqnControl["WallBoundaryType"] ;
	}

	/*--- Neutral equation control ---*/
	if ( config.find ("NeutralEqnControl") != config.end() ) 
	{ 
		NeutralEqnControl = config[ "NeutralEqnControl" ] ; 

		if ( NeutralEqnControl.find ("Equation") 	!= NeutralEqnControl.end() ) 
			Equation[ NEUTRAL ].Equation 	= NeutralEqnControl["Equation"] ;
		
		if ( NeutralEqnControl.find ("WallBoundaryType")  != NeutralEqnControl.end() ) 
			Equation[ NEUTRAL ].WallBoundaryType 	= NeutralEqnControl["WallBoundaryType"] ;	

	}
	/*--- Background equation control ---*/
	if ( config.find ("BackgroundEqnControl") != config.end() ) 
	{ 
		BackgroundEqnControl = config[ "BackgroundEqnControl" ] ; 

		if ( BackgroundEqnControl.find ("Equation") 	!= BackgroundEqnControl.end() ) 
			Equation[ BACKGROUND ].Equation 	= BackgroundEqnControl["Equation"] ;	
	}


	/*--- Output Control ---*/
	if ( config.find ("OutputControl") != config.end() ) OutputControl = config[ "OutputControl" ] ; 

	if ( OutputControl.find ("WRT_Insta_Freq")  != OutputControl.end() ) WRT_Insta_Freq 	= OutputControl["WRT_Insta_Freq"] ;
	if ( OutputControl.find ("WRT_Cycle_Freq") 	!= OutputControl.end() ) WRT_Cycle_Freq 	= OutputControl["WRT_Cycle_Freq"] ;
	if ( OutputControl.find ("MON_Insta_Freq")  != OutputControl.end() ) MON_Insta_Freq 	= OutputControl["MON_Insta_Freq"] ;
	if ( OutputControl.find ("MON_Cycle_Freq") 	!= OutputControl.end() ) MON_Cycle_Freq 	= OutputControl["MON_Cycle_Freq"] ;

	if ( OutputControl.find ("AVERAGE_SW") 		!= OutputControl.end() ){
		if ( OutputControl["AVERAGE_SW"] == 1 ){
			Average_Switch 	= true ;
		} else {
			Average_Switch = false ;
		}
	} 

	if ( OutputControl.find ("OutputFormat") 		!= OutputControl.end() ){
		if ( OutputControl["OutputFormat"] == "TECPLOT" ){
			OutputFormat 	= 0 ;
		} else if ( OutputControl["OutputFormat"] == "TECPLOT1D" ){
			OutputFormat 	= 2 ;
		}
	} else {
		OutputFormat 	= 0 ;
	}

}

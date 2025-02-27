#include "variable_structure.hpp"
using namespace std ;
CVariable::CVariable()
{
	bool Normalize = false ;

	if ( Normalize ) {

		Ref_L 		= 1.0 ;
		Ref_N 		= 1.E+15 ;
		Ref_Mass 	= 9.1095E-31 ;
		Ref_Qe 		= 1.602e-19 ;
		Ref_Kb 		= 1.38065E-23 ;
		Ref_Phi 	= 1.0 ;

	} else {

		Ref_L 		= 1.0 ;
		Ref_N 		= 1.0 ;
		Ref_Mass 	= 1.0 ;
		Ref_Qe 		= 1.0 ;
		Ref_Kb 		= 1.0 ;
		Ref_Phi 	= 1.0 ;

	}
	double coeff = Ref_Qe/Ref_Kb ;
	Ref_Te 		= Ref_Phi ; // [ev]
	//RefVelocity = sqrt( Refkb*(RefTe*coeff)/RefMass ) ;// [m/s]
	Ref_V 		= sqrt( Ref_Kb*(Ref_Te*coeff)/Ref_Mass ) ;
	Ref_Flux 	= Ref_V*Ref_N ;
	Ref_t 		= Ref_L/Ref_V ;
	//RefDiff = RefLength*RefVelocity ; //[m^2/s]
	Ref_Diff 	= Ref_L*Ref_V ; 
	//RefMobi = RefLength* RefQe / RefMass / RefVelocity ;
	Ref_Mu 		= Ref_L*Ref_Qe/Ref_Mass/Ref_V ;
	//RefPermittivity 	= (RefQe*RefQe)* (RefLength*RefLength) * RefDensity / (RefMass*RefVelocity*RefVelocity );
	Ref_Eps 	= (Ref_Qe*Ref_Qe)*(Ref_L*Ref_L)*Ref_N/ (Ref_Mass*Ref_V*Ref_V ) ;
	//RefEnergyDensity  	= RefTe * RefDensity ;
	Ref_EN 		= Ref_Te * Ref_N ;
	Ref_Rho 	= Ref_Qe*Ref_N ;
	Ref_SQ 		= Ref_Phi*Ref_Eps/Ref_L ;
	//RefSurfaceCharge  	= Ref_Phi*Ref_Eps/Ref_L ;
	Ref_JD 		= Ref_Eps*(Ref_Phi/Ref_L)/Ref_t ;
	Ref_ES 		= Ref_EN/Ref_t ;
  Ref_SS  	= Ref_N * Ref_V/Ref_L ;
  Ref_EField	= Ref_Phi/Ref_L ;
}
void CVariable::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	if( mpi_rank == MASTER_NODE ) cout<<"Start buliding variable module ... "<<endl;
	/*--- Read Table: for calculate electron transport coefficient ---*/
	if ( config->Species[ 0 ].MobilityType == 1 ){
		CollTable.Init( config->CasePath+config->Species[ 0 ].CollisionFreqFile ) ;
	}

	/*--- Read Table: when using LFA, mean energy is from table ---*/
	if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "NONE" ) {

		MeanEnergyTable.Init( config->CasePath+"4MeanEnergy.inp" ) ;

	}else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "2Fluid" ){

		AlphaTable.Init( config->CasePath+"6Alpha.inp" ) ;
		MeanEnergyTable.Init( config->CasePath+"4MeanEnergy.inp" ) ;

	}else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "3Fluid" ){

		AlphaTable.Init( config->CasePath+"6Alpha.inp" ) ;
		EtaTable.Init( config->CasePath+"7Eta.inp"   ) ;
		MeanEnergyTable.Init( config->CasePath+"4MeanEnergy.inp" ) ;

	}else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "Streamer" ){
		AlphaTable.Init( config->CasePath+"6Alpha.inp" ) ;
	}


	/*--- If electron energy loss term is calculate from file ---*/
	if ( config->eEnergyLossFile.size() > 0 ) {

//		eEnergyLossTable.Init( config->CasePath+"eEnergyLoss.inp" ) ;

	}

	/*------*/
	Qe = 1.602e-19/Ref_Qe ;
 	PI = 4.0*atan(1.0) ;
	Me = 9.1095E-31/Ref_Mass ;//[ kg ]
	Kb = 1.38065E-23/Ref_Kb ;//[ m2 kg s-2 K-1]
	K2eV = Qe/Kb ;
	PhysicalTime = 0.0 ;

	if ( config->DC_TimeStep > 0.0 ) {

	 	Dt = config->DC_TimeStep ;
	 	//config->ExitCycle = 2 ;
	} else {

	 	double AC_Period = 1.0 / config->MaxFrequency ;
	 	Dt =  AC_Period / config->StepPerCycle ;
	 	Dt /= Ref_t ;

	}

	/*--- Could be delete ---*/
		Momentum_Term[ 0 ].initial ( "M_Convection" ) ;
		Momentum_Term[ 1 ].initial ( "M_Force" ) ;
		Momentum_Term[ 2 ].initial ( "M_Collision" ) ;

		Energy_Term[ 0 ].initial ( "Convection" ) ;
		Energy_Term[ 1 ].initial ( "Conduction" ) ;
		Energy_Term[ 2 ].initial ( "Joule" ) ;
		Energy_Term[ 3 ].initial ( "Collision" ) ;
		Energy_Term[ 4 ].initial ( "Pressure" ) ;
		Energy_Term[ 5 ].initial ( "velocity" ) ;
	/*--- Could be delete ---*/


	if ( mpi_rank == MASTER_NODE ) {
		cout<<"Max. frequency: "<< config->MaxFrequency <<" Hz."<<endl;
		cout<<"Time setp size: "<< Dt <<" sec."<<endl;
	} 

	/*--- Chemistry Module ---*/
	Chemistry.init( config->CasePath + "chemistry_parameter.txt", plasma.Mesh.cell_number ) ;
	//exit(1) ;
	if ( (Chemistry.species_size-config->TotalSpeciesNum) != 0 ) {
		cout<<"config->TotalSpeciesNum: "<<config->TotalSpeciesNum<<endl;
		cout<<"Chemistry.species_size: "<<Chemistry.species_size<<endl;
		cout<<"Species Number doesn't not match pls check Chemistry and Configure"<<endl;exit(1) ;
	}
	//exit(1) ;
	//if( mpi_rank == 0 ) cout<< (Chemistry.species_size-config->TotalSpeciesNum)<<endl;
	//exit(1) ;
	InputNumberDensity = new double[ plasma.Mesh.cell_number*config->TotalSpeciesNum ] ;
	InputTemperature   = new double[ plasma.Mesh.cell_number*config->TotalSpeciesNum ] ;
	InputEField        = new double[ plasma.Mesh.cell_number*config->TotalSpeciesNum ] ;
	coll_freq          = new double[ plasma.Mesh.cell_number ] ;

	//MobilityPoint		= Chemistry.ptr_mobility();
	//DiffusivityPoint	= Chemistry.ptr_diffusion();
	//chemistry module - calculate transport coefficients, chemical source and electron energy loss
	ReactionRatePoint	= Chemistry.ptr_source_sink();
	EnergySourcePoint	= Chemistry.ptr_energy_loss();
	if ( config->eEnergyLossFile.size() > 0 ) {
		//Lookup table
	} else {

	}


		//ElectricalMap.insert( pair< int, CElectrical>( container.Type, container ) ) ;
		//ElectricalIter = ElectricalMap.find( container.Type ) ;
		//cout<<"iter: " <<ElectricalIter->second.Frequency<<endl;
		CTable container ;

		map< int, CTable>::iterator MobilityIter, DiffusivityIter ;

    	for ( int iSpecies = 0 ; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) {

    		if ( config->Species[ iSpecies ].MobilityType == 2 
    		or 	 config->Species[ iSpecies ].MobilityType == 3 
    		or   config->Species[ iSpecies ].MobilityType == 4 
    		or   config->Species[ iSpecies ].MobilityType == 6){

   		 		MobilityMap.insert( pair< int, CTable>( iSpecies, container ) ) ;
    			MobilityIter = MobilityMap.find( iSpecies ) ;
				MobilityIter->second.Init( config->CasePath + config->Species[ iSpecies ].MobilityFile) ;
				//cout<<"Mu: "<<MobilityIter->second.GetValueLog( 100. )<<endl;
    		}

    		if ( config->Species[ iSpecies ].DiffusivityType == 2 
    		or   config->Species[ iSpecies ].DiffusivityType == 3
    		or   config->Species[ iSpecies ].DiffusivityType == 4
    		or   config->Species[ iSpecies ].DiffusivityType == 6 ){
    		
   		 		DiffusivityMap.insert( pair< int, CTable>( iSpecies, container ) ) ;
    			DiffusivityIter = DiffusivityMap.find( iSpecies ) ;
				DiffusivityIter->second.Init( config->CasePath + config->Species[ iSpecies ].DiffusivityFile) ;
    		}
    	}
		//cout<<"Variable: Mobi. & Diff. mpa initialzation"<<endl ; exit(1) ;

		MPI_ID.initial ( "MPI_ID" ) ;
		DebyeLength.initial ( "λ<sub>D</sub>/L" ) ;
		AvgDebyeLength.initial ( "λ<sub>D</sub>/L" ) ;

		CFL.initial ( "CFL" ) ;
		AvgCFL.initial ( "CFL" ) ;

		

	/*--- Electrical Field Variables ---*/
		//Phi.initial ( "Φ [V]" ) ;


		UltraMPPVarInit()	;

		UltraMPPAvgVarInit() ;

	//	for ( auto mpp = VarTag.cbegin(); mpp != VarTag.cend(); ++mpp ) {
		//	cout<<mpp->first<<"\t"<<mpp->second<<endl;
		//}
		//exit(0);


	/* End UltraMPP Variables */



		 //Scalar eEnergyLoss, eAvgEnergyLoss ;
		eEnergyLoss.initial ( "ε_Loss" ) ;
		eAvgEnergyLoss.initial ( "ε_Loss" ) ;

		TotalNumberDensity.initial ( "Total Number Density" ) ;

		for ( int nDim = 0 ; nDim < 3 ; nDim++ ) {
			  // EField[ nDim ].initial( "E_"+to_string(nDim) ) ;
			//PreEField[ nDim ].initial( "PreE_"+to_string(nDim) ) ;
		}

	/*--- Transport Coefficients ---*/
		Mobi = new CScalar [ config->TotalSpeciesNum ] ;
		for ( int i = 0 ; i < config->TotalSpeciesNum ; i++ ) Mobi[ i ].initial( "μ_"+config->Species[ i ].Name ) ;

		Diff = new CScalar [ config->TotalSpeciesNum ] ;
		for ( int i = 0 ; i < config->TotalSpeciesNum ; i++ ) Diff[ i ].initial( "D_"+config->Species[ i ].Name  ) ;

		Debug	= new CScalar [ config->TotalSpeciesNum ] ;
		for ( int i = 0 ; i < config->TotalSpeciesNum ; i++ ) Debug[ i ].initial( "Debug"  ) ;
			
		ProductionRate = new CScalar [ config->TotalSpeciesNum ] ;
		for ( int i = 0 ; i < config->TotalSpeciesNum ; i++ ) ProductionRate[ i ].initial( "Rdot-"+config->Species[ i ].Name  ) ;
			//exit(1) ;
	/*--- Solution Variables ---*/
		DD_Convection 		= new CScalar [ config->TotalSpeciesNum ] ;

		T 		= new CScalar [ config->TotalSpeciesNum ] ;
		AvgT	= new CScalar [ config->TotalSpeciesNum ] ;
		PreT 	= new CScalar [ config->TotalSpeciesNum ] ;

		U0 		= new CScalar [ config->TotalSpeciesNum ] ;
		PreU0 = new CScalar [ config->TotalSpeciesNum ] ;
		AvgU0 = new CScalar [ config->TotalSpeciesNum ] ;

		U1 		= new CScalar [ config->TotalSpeciesNum ] ;
		PreU1 = new CScalar [ config->TotalSpeciesNum ] ;
		AvgU1 = new CScalar [ config->TotalSpeciesNum ] ;
		
		U2 		= new CScalar [ config->TotalSpeciesNum ] ;
		PreU2 = new CScalar [ config->TotalSpeciesNum ] ;
		AvgU2 = new CScalar [ config->TotalSpeciesNum ] ;

		U3 		= new CScalar [ config->TotalSpeciesNum ] ;
		PreU3 = new CScalar [ config->TotalSpeciesNum ] ;
		AvgU3 = new CScalar [ config->TotalSpeciesNum ] ;
		
		U4 		= new CScalar [ config->TotalSpeciesNum ] ;
		PreU4 = new CScalar [ config->TotalSpeciesNum ] ;
		AvgU4 = new CScalar [ config->TotalSpeciesNum ] ;

		JouleHeating 		= new CScalar [ config->TotalSpeciesNum ] ;
		AvgJouleHeating = new CScalar [ config->TotalSpeciesNum ] ;

		
		LFASourceSink 	= new CScalar [ 3 ] ;
		LFASourceSink[ 0 ].initial( "LFASourceSink_"+config->Species[ 0 ].Name  ) ;
		LFASourceSink[ 1 ].initial( "LFASourceSink_"+config->Species[ 1 ].Name  ) ;
		LFASourceSink[ 2 ].initial( "LFASourceSink_"+config->Species[ 2 ].Name  ) ;

		for ( int i = 0 ; i < config->TotalSpeciesNum ; i++ ) {

				DD_Convection[ i ].initial( "DD_Convection_"+config->Species[ i ].Name  ) ;
			    T[ i ].initial( "T_"+config->Species[ i ].Name  ) ;
			   U0[ i ].initial( "N_"+config->Species[ i ].Name ) ;
			   U1[ i ].initial( "U1_"+config->Species[ i ].Name ) ;
			   U2[ i ].initial( "U2_"+config->Species[ i ].Name ) ;
			   U3[ i ].initial( "U3_"+config->Species[ i ].Name ) ;
			   U4[ i ].initial( "U4_"+config->Species[ i ].Name ) ;
			JouleHeating[ i ].initial( "JdotE_"+config->Species[ i ].Name ) ;

			 PreT[ i ].initial( "TO"+config->Species[ i ].Name ) ;
			PreU0[ i ].initial( "NO_"+config->Species[ i ].Name  ) ;
			PreU1[ i ].initial( "U1O_"+config->Species[ i ].Name  ) ;
			PreU2[ i ].initial( "U2O_"+config->Species[ i ].Name  ) ;
			PreU3[ i ].initial( "U3O_"+config->Species[ i ].Name  ) ;
			PreU4[ i ].initial( "N4O_"+config->Species[ i ].Name  ) ;

			 AvgT[ i ].initial( "Tavg_"+config->Species[ i ].Name  ) ;
			AvgU0[ i ].initial( "AvgN_"+config->Species[ i ].Name  ) ;
			AvgU1[ i ].initial( "AvgUx_"+config->Species[ i ].Name  ) ;
			AvgU2[ i ].initial( "AvgUy_"+config->Species[ i ].Name  ) ;
			AvgU3[ i ].initial( "AvgUz_"+config->Species[ i ].Name  ) ;
			AvgU4[ i ].initial( "AvgNε_"+config->Species[ i ].Name  ) ;
			AvgJouleHeating[ i ].initial( "AvgJdotE_"+config->Species[ i ].Name ) ;
		}

		GradU0 	= new CScalar *[ config->TotalSpeciesNum ] ;
		GradU4 	= new CScalar *[ config->TotalSpeciesNum ] ;
		for ( int iSpecies = 0 ; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) {
			GradU0[ iSpecies ] = new CScalar [ 3 ] ;
			GradU4[ iSpecies ] = new CScalar [ 3 ] ;
			for ( int i = 0 ; i < 3 ; i++ ){
				GradU0[iSpecies][ i ].initial( "GradientN["+to_string(i)+"]_"+config->Species[ iSpecies ].Name ) ;			
				GradU4[iSpecies][ i ].initial( "GradientN<sub><greek>e</greek></sub>["+to_string(i)+"]_"+config->Species[ iSpecies ].Name ) ;			
			}
		} 

		GradT 	= new CScalar *[ config->TotalSpeciesNum ] ;
		for ( int iSpecies = 0 ; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) {
			GradT[ iSpecies ] = new CScalar [ 3 ] ;
			for ( int i = 0 ; i < 3 ; i++ ){
				GradT[iSpecies][ i ].initial( "GradT["+to_string(i)+"]_"+config->Species[ iSpecies ].Name ) ;			
			}
		} 
		Beta.initial( "β" ) ;

		/*--- Current density ---*/	
		CondJD 	= new CScalar *[ config->ChargeSpeciesNum ] ;

		CondI_PowerElectrode_local 	= new double [ config->ChargeSpeciesNum ] ;
		CondI_GroundElectrode_local	= new double [ config->ChargeSpeciesNum ] ;

		CondI_PowerElectrode_global_sum	= new double [ config->ChargeSpeciesNum ] ;
		CondI_GroundElectrode_global_sum= new double [ config->ChargeSpeciesNum ] ;

		for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ) {
			CondJD[ iSpecies ] = new CScalar [ 3 ] ;
			for ( int nDim = 0 ; nDim < 3 ; nDim++ ) {
				CondJD[iSpecies][ nDim ].initial( "CondJD"+to_string(nDim)+"_"+config->Species[ iSpecies ].Name ) ;
			}
		}

		for ( int nDim = 0 ; nDim < 3 ; nDim++ ) {
			 DispJD[ nDim ].initial( "DispJD_"+to_string(nDim) ) ;
			TotalJD[ nDim ].initial( "TotalJD_"+to_string(nDim) ) ;
		}

		TotalGasPressure.initial( "P<sub>gas</sub>" ) ;
		InitialConditions( m, config ) ;

		CalTotalPressure( m, config ) ;
		ChemistryUpdate( m, config ) ;



	/*--- Material Properties ---*/
		Eps .initial( "Effective ε" ) ;
		Eps0.initial( "ε" ) ;
		NetQ.initial( "NetQ" ) ;
		CellProperties( m ) ; Eps0 =Eps0 ;
		if( mpi_rank == MASTER_NODE ) cout<<"End buliding variable module ... "<<endl;

}
void CVariable::UltraMPPVarInit()
{
		plasma.set_parallel_variable( &Potential, "Potential" );
		for ( int i = 0; i < plasma.Mesh.face_number; i++ ) {
			Potential.face[i] = 0.0 ;
		}
		VarTag["permittivity"          ] = plasma.set_parallel_cell_data(&eps      , "permittivity"           ) ;
		VarTag["effective_permittivity"] = plasma.set_parallel_cell_data(&eps_eff  , "effective_permittivity" ) ;
		VarTag["effective_permittivity"] = plasma.set_parallel_cell_data(&eps_eff_A  , "effective_permittivity_A" ) ;
		VarTag["ChargeDen"             ] = plasma.set_parallel_cell_data(&ChargeDen, "ChargeDen"              ) ;
		VarTag["ChargeDen    [m^-3]"   ] = plasma.set_parallel_cell_data(&RealChargeDen, "ChargeDen    [m^-3]"              ) ;
		
		VarTag["Ex"   ] = plasma.set_parallel_cell_data(    &Ex, "Ex  [V/m]" ) ;
		VarTag["Ey"   ] = plasma.set_parallel_cell_data(    &Ey, "Ey  [V/m]" ) ;
		VarTag["Ez"   ] = plasma.set_parallel_cell_data(    &Ez, "Ez  [V/m]" ) ;

		VarTag["PreEx"] = plasma.set_parallel_cell_data( &PreEx, "Ex [V/m]" ) ;
		VarTag["PreEy"] = plasma.set_parallel_cell_data( &PreEy, "Ey [V/m]" ) ;
		VarTag["PreEz"] = plasma.set_parallel_cell_data( &PreEz, "Ez [V/m]" ) ;

		VarTag["Etd" ] = plasma.set_parallel_cell_data(  &Etd,  "Etd" ) ;
		VarTag["Emag"] = plasma.set_parallel_cell_data( &Emag, "Emag" ) ;
		VarTag["Tx"] = plasma.set_parallel_cell_data( &Tx, "Tx" ) ;
		VarTag["Ty"] = plasma.set_parallel_cell_data( &Ty, "Ty" ) ;
		VarTag["Tz"] = plasma.set_parallel_cell_data( &Tz, "Tz" ) ;
		//GradTe
		VarTag["Kappa"] = plasma.set_parallel_cell_data( &Kappa, "kappa" ) ;

		VarTag["plot_var"] = plasma.set_parallel_cell_data( &plot_var, "plot_var" ) ;

		VarTag["surface_charge"] = plasma.set_face_data( &surface_charge, "surface_charge" ) ;


		//---------Variable definition of ICP simulation ---------------------------
		#if ( FDMaxwell == true ) 
		FDMaxwell_coupled_eqs.set_parallel_variable(&E_phi_Re,"E_phi_Re") ;
    	FDMaxwell_coupled_eqs.set_parallel_variable(&E_phi_Im,"E_phi_Im") ;
    	FDMaxwell_Re.set_multi_variable( 0, E_phi_Re.tag_current );
    	FDMaxwell_Im.set_multi_variable( 0, E_phi_Im.tag_current );
	    FDMaxwell_coupled_eqs.set_multi_variable( 0, E_phi_Re.tag_current );
	    FDMaxwell_coupled_eqs.set_multi_variable( 1, E_phi_Im.tag_current );
		VarTag["sigma_p_Re_plasma" ] = plasma.set_parallel_cell_data		  (    &sigma_p_Re_plasma 		, "sigma_p_Re_plasma" 			) ;
		VarTag["sigma_p_Im_plasma" ] = plasma.set_parallel_cell_data		  (    &sigma_p_Im_plasma 		, "sigma_p_Im_plasma" 			) ;
		VarTag["sigma_p_Re_FVFD"   ] = FDMaxwell_Re.set_parallel_cell_data(    &sigma_p_Re_FVFD			, "sigma_p_Re_FVFD" 			) ;
		VarTag["sigma_p_Im_FVFD"   ] = FDMaxwell_Re.set_parallel_cell_data(    &sigma_p_Im_FVFD			, "sigma_p_Im_FVFD" 			) ;
		VarTag["Re_eq_source"      ] = FDMaxwell_Re.set_parallel_cell_data(    &Re_eq_source				, "Re_eq_source" 				) ;
		VarTag["Im_eq_source"      ] = FDMaxwell_Re.set_parallel_cell_data(    &Im_eq_source				, "Im_eq_source" 				) ;
		VarTag["permittivity_FVFD"       ] = FDMaxwell_Re.set_parallel_cell_data(	 &eps_FVFD      			, "permittivity_FVFD"           ) ;
		VarTag["k_square_Re"  			 ] = FDMaxwell_Re.set_parallel_cell_data(    &k_square_Re			 	, "k_square_Re" 				) ;
		VarTag["k_square_Im"  			 ] = FDMaxwell_Re.set_parallel_cell_data(    &k_square_Im			 	, "k_square_Im" 				) ;
		VarTag["CurrentDen"  		 	 ] = FDMaxwell_Re.set_parallel_cell_data(    &CurrentDen			 	, "CurrentDen" 					) ;
		VarTag["Power_Absorption_FVFD"   ] = FDMaxwell_Re.set_parallel_cell_data(    &Power_Absorption_FVFD	 	, "Power_Absorption_FVFD" 		) ;
		VarTag["Power_Absorption_plasma" ] = plasma.set_parallel_cell_data		(    &Power_Absorption_plasma	,"Power_Absorption_plasma" 		) ;	
		#endif	
		

}
void CVariable::UltraMPPAvgVarInit()
{
	VarTag["AvgPotential"]  = plasma.set_parallel_cell_data(   &AvgPotential, "Φ [V]" ) ;
	VarTag["AvgEx"  ] = plasma.set_parallel_cell_data(   &AvgEx, "Ex [V/m]" ) ;
	VarTag["AvgEy"  ] = plasma.set_parallel_cell_data(   &AvgEy, "Ey [V/m]" ) ;
	VarTag["AvgEz"  ] = plasma.set_parallel_cell_data(   &AvgEz, "Ez [V/m]" ) ;
	VarTag["AvgE_Td"] = plasma.set_parallel_cell_data( &AvgE_Td, "reduce E (Td)" ) ;
	VarTag["avg_plot_var"] = plasma.set_parallel_cell_data( &avg_plot_var, "avg_plot_var" ) ;
}
void CVariable::UltraMPPInitialCellParameter()
{
  for ( int cth=0 ; cth<plasma.Mesh.cell_number; cth++ ) {
    Cell *cell = plasma.get_cell( cth ) ;
    eps[ cth ] = cell_parameter[ plasma.get_cell_typename( cell->data_id ) ] * vacuum_permittivity ;
  }
  plasma.syn_parallel_cell_data( VarTag["permittivity"] );
}
void CVariable::CellProperties( boost::shared_ptr<CDomain> &m )
{

	// Cell *Cell_i ;

	// for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	Cell_i  = plasma.get_cell( i ) ;

	// 	if( Cell_i->type == MPP_cell_tag["PLASMA"] ) {

	// 		Eps0[ i ] = 1.0*vacuum_permittivity/Ref_Eps ;

	// 	}	else if ( Cell_i->type == MPP_cell_tag["POWER"] ) {

	// 		Eps0[ i ] = 1.0E+10 /Ref_Eps ;

	// 	} else if ( Cell_i->type == MPP_cell_tag["GROUND"] ) {

	// 		Eps0[ i ] = 1.0E+10 /Ref_Eps ;

	// 	}	else if ( Cell_i->type == MPP_cell_tag["DIELECTRIC"] ) {

	// 		Eps0[ i ] = 4.0*vacuum_permittivity/Ref_Eps ;
	// 	}
	// 	Eps[ i ] = Eps0[ i ] ;

	// }
	// Eps0 = Eps0 ;
	// Eps  = Eps  ;
}
void CVariable::InitialConditions( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	UltraMPPInitialCellParameter() ;
	
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		MPI_ID[ i ] = mpi_rank ;
	}//End cell loop
	double N0 = 5.0E+18/Ref_N ;
	double sigma=0.4/1000.0/Ref_L;
	double r0=0.0, z0=1.0/100.0/Ref_L, r=0.0, z=0.0 ;
	Cell *Cell_i ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		TotalNumberDensity[ i ] = 0.0 ;

		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			for ( int jSpecies = 0; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

				U0[ jSpecies ][ i ] = config->Species[jSpecies].InitialDensity/Ref_N ; 
				U1[ jSpecies ][ i ] = config->Species[jSpecies].InitialDensity*0.E-15/Ref_Flux ;
				U2[ jSpecies ][ i ] = config->Species[jSpecies].InitialDensity*0.E-15/Ref_Flux ;
				U3[ jSpecies ][ i ] = config->Species[jSpecies].InitialDensity*0.E-15/Ref_Flux ;//config->Species[jSpecies].InitialDensity ;
			}

			if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "Streamer" ) {
				r = Cell_i->r[0] ;
				z = Cell_i->r[1] ;
				U0[ 1 ][ i ] += N0*exp( -(r*r+pow(z-z0, 2.0)) / pow(sigma, 2.0) ) ;
			}

		}//End plasma cells
	}//End cell loop
	//double *tmp ;
	//int tag_tmp   =	plasma.set_parallel_cell_data( &tmp, "N_Ar" ) ;

	/*read the initial conditions from file */
	for ( int i = 0; i < config->TotalSpeciesNum ; i++ ) {
		if( !config->Species[i].initial_data_file_name.empty() )
		{
			plasma.set_initial_variable( U0[ i ].data_id, config->Species[i].initial_data_variable_name  ) ;
			//plasma.set_initial_variable( tag_tmp ) ;
			plasma.get_initial_value( config->Species[i].initial_data_file_name ) ;
			//cout<<config->Species[i].initial_data_file_name<<endl;
			//cout<<U0[ i ].data_name<<endl;
			// plasma.set_output( "tmp" ) ;
			// plasma.set_output( tag_tmp ) ;
			// plasma.write_output( "test-restart" ) ;
			//cout<<"not support yet"<<endl; 
			//PetscEnd() ;
		}
	}




	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		Cell_i  = plasma.get_cell( i ) ;
		TotalNumberDensity[ i ] = 0.0 ;
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){
			for ( int jSpecies = 0; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
				TotalNumberDensity[ i ] += U0[ jSpecies ][ i ]*Ref_N ;
			}
		}//End plasma cells
	}//End cell loop
	TotalNumberDensity = TotalNumberDensity ;	






	/*---  Temperature ---*/
	double Internal_E=0.0, Kinetic_E=0.0 ;
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
		Cell_i  = plasma.get_cell( i ) ;
		if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

			for ( int jSpecies = 0; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

				T[ jSpecies ][ i ] = config->Species[jSpecies].InitialTemperature/Ref_Te ;
				
				//cout<< config->Species[jSpecies].InitialTemperature<<endl;

				//cout<<Ref_Te<<endl;

				if ( config->PFM_Assumption != "LFA" ) {
					U4[ jSpecies ][ i ] = (3.0/2.0)*U0[ jSpecies ][ i ]*T[jSpecies][ i ] ;
					// Internal_E = (U0[ jSpecies ][ i ]*T[jSpecies][ i ]*Qe)/(config->Species[ jSpecies ].Gamma-1.0) ;
					//  Kinetic_E = 0.5*(config->Species[ jSpecies ].Mass_Kg/Ref_Mass)*( U1[jSpecies][ i ]*U1[jSpecies][ i ]+U2[jSpecies][ i ]*U2[jSpecies][ i ]+U3[jSpecies][ i ]*U3[jSpecies][ i ])/U0[jSpecies][ i ] ;
					// U4[jSpecies][ i ] =  Internal_E + Kinetic_E ;
				}

			} 
			Beta[ i ] = 1.0 ;
		}	
	}


	/*--- Update overlap cell value ---*/
	for ( int jSpecies = 0; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
		U0[ jSpecies ] = U0[ jSpecies ] ;
		U1[ jSpecies ] = U1[ jSpecies ] ;
		U2[ jSpecies ] = U2[ jSpecies ] ;
		U3[ jSpecies ] = U3[ jSpecies ] ;
		U4[ jSpecies ] = U4[ jSpecies ] ;
		 T[ jSpecies ] =  T[ jSpecies ] ;
	}
	Beta = Beta ;
	//cout<<"UpdateSolution"<<endl;
	UpdateSolution( m ) ;
}
void CVariable::CalTotalPressure( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	double T_kelvin=0.0, P_pascal=0.0 ;
	for ( int jSpecies = 0; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
		if ( config->Species[ jSpecies ].Type == BACKGROUND )
		for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
			T_kelvin = T[ jSpecies ][ i ] * 11604.52500617 ;
			P_pascal = (Kb*Ref_Kb)*(U0[ jSpecies ][ i ]*Ref_N)*T_kelvin ; //[Pa]=[kg/m/s^2]
			TotalGasPressure[ i ] = P_pascal*0.00750061683 ; //[torr]
		}
	}
}
void CVariable::ChemistryUpdate( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config  )
{

	Cell *Cell_i ;
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		if (plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			TotalNumberDensity[ i ] = 0.0 ;

			for( int k = 0 ; k < config->TotalSpeciesNum ; k++ ){
				*(InputNumberDensity + ( (i)*config->TotalSpeciesNum + k ) ) = U0[ k ][ i ]*Ref_N ;
				//cout<<"InputNumberDensity: "<<*(InputNumberDensity + ( (i)*config->TotalSpeciesNum + k ) )<<endl;
				if ( config->Species[k].Type == BACKGROUND ) TotalNumberDensity[ i ] += U0[ k ][ i ]*Ref_N ;
			}
		}//End plasma cells
	}//End cell loop
	//exit(1);


	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
			for ( int iSpecies = 0 ; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) {

				if( iSpecies == 0 )	{
					/* electron */
					*(InputTemperature + ( (i)*config->TotalSpeciesNum + iSpecies ) ) = T[ iSpecies ][ i ] ;

				}	else {
					/* other species */
					*(InputTemperature + ( (i)*config->TotalSpeciesNum + iSpecies ) ) = T[ iSpecies ][ i ]*11604.52500617 ;
					//cout<<*(InputTemperature + ( (i)*config->TotalSpeciesNum + iSpecies ) )<<endl; exit(1) ;
				}

			}//End iSpecies
		}//End plasma cells
	}//End cell loop



	if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "2Fluid" ) {
		CalMeanEnergy( m ) ;
	} else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "3Fluid" ){

	} else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "Streamer" ){
		//CalMeanEnergy( m ) ;
	} else if( config->PFM_Assumption == "LMEA" ) {
		//cout<<"Here"<<endl;
	    Chemistry.CalTotalPressure( InputNumberDensity, InputTemperature ) ;
	    Chemistry.CalSourceSinkRate_for_Te( InputTemperature ) ;
	    Chemistry.CalSourceSinkRate_for_gas( InputTemperature ) ;

	} else {
		cout<<"Error in variable_structure.cpp, ChemistryUpdate"<<endl;
		cout<<"PFM_Assumption -> "<<config->PFM_Assumption<<endl;
		cout<<"PFM_SubModel   -> "<<config->PFM_SubModel<<endl;
	}


	/*--- Update electron transport coefficients ---*/
	UpdateElectronTransport( m, config ) ;

	if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "Streamer" ){
	//
	} else {
		UpdateIonNeutralTransport( m, config ) ;
	}
	
	/* Update electron transport coefficients */
	if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "2Fluid" ) {

	SourceSink_2Fluid( m, config ) ;

	} else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "3Fluid" ) {

		SourceSink_3Fluid( m, config ) ;

	} else if ( config->PFM_Assumption == "LFA" and config->PFM_SubModel == "Streamer" ){

		SourceSink_Streamer( m, config ) ;
		//SourceSink_PSST_2018( m, config ) ;

	} else {

		Chemistry.SourceSink( InputNumberDensity, InputTemperature ) ;
		ReactionRatePoint	= Chemistry.ptr_source_sink() ;
		EnergySourcePoint	= Chemistry.ptr_energy_loss() ;
	}
}
void CVariable::UltraMPPComputeReducedElectricField()
{
	//Note: ReduceElectricField has unit.

	 
	double E=0.0 ;
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell *cell = plasma.get_cell( i ) ;

		if ( cell_type[ cell->type ] == PLASMA ) {
			Emag[ i ] = sqrt ( Ex[ i ]*Ex[ i ] + Ey[ i ]*Ey[ i ] + Ez[ i ]*Ez[ i ] )*Ref_EField+ZERO ;
			Etd[ i ] = Emag[ i ]/TotalNumberDensity[ i ]/(1.0E-21) ; 
		}
	}
	plasma.syn_parallel_cell_data( VarTag["Etd"] );
	plasma.syn_parallel_cell_data( VarTag["Emag"] );
}
void CVariable::CalMeanEnergy( boost::shared_ptr<CDomain> &m )
{
	double C23 = 2.0/3.0 ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		//U4[ 0 ][ i ] = MeanEnergyTable.GetValueLog( Etd[ i ] )/Ref_EN ;
		U4[ 0 ][ i ] = MeanEnergyTable.GetValue( Emag[ i ] ) ;
		 T[ 0 ][ i ] = U4[ 0 ][ i ]*C23 ;
	}
	U4[ 0 ] = U4[ 0 ] ;
	 T[ 0 ] =  T[ 0 ] ;
}
void CVariable::UpdateSolution( boost::shared_ptr<CDomain> &m )
{

	for( int k = 0 ; k < Chemistry.species_size ; k++ ) {

		for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
			PreU0[ k ][ i ] = U0[ k ][ i ] ;
			PreU1[ k ][ i ] = U1[ k ][ i ] ;
			PreU2[ k ][ i ] = U2[ k ][ i ] ;
			PreT [ k ][ i ] =  T[ k ][ i ] ;
			PreU4[ k ][ i ] = U4[ k ][ i ] ;
		}

		PreU0[ k ] = PreU0[ k ] ;
		PreU1[ k ] = PreU1[ k ] ;
		PreU2[ k ] = PreU2[ k ] ;
		PreT [ k ] = PreT [ k ] ;
		PreU4[ k ] = PreU4[ k ] ;

		if( nDim == 3 ) {
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				PreU3[ k ][ i ] = U3[ k ][ i ] ;
			}
			PreU3[ k ] = PreU3[ k ] ;
		}
	}

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		PreEx[i] = Ex[i] ;
		PreEy[i] = Ey[i] ;
		PreEz[i] = Ez[i] ;
	}
}
void CVariable::UpdateElectronTransport( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	double C1=0.0, C2=0.0, Thermal=0.0, DriftV=0.0 ;
	// double BackGroundNumberDensity = 0.0, Pressure=0.5, BackGroundTemp=300.0 ; 
	// BackGroundNumberDensity = 133.2894 * Pressure / 1.38065E-23 / BackGroundTemp ;
	double collision = 0.0 ;
	double E_Td=0.0, EMag=0.0, P_torr = 0.0, DriftVelocicy=0.0, ReduceMobi=0.0, Etd_Min=0.0 ;
	map< int, CTable>::iterator MobilityIter, DiffusivityIter ;
	Cell *Cell_i ;

	switch ( config->Species[ 0 ].MobilityType ){
		
		case 0://cout<<"Constant mobility." <<endl ;
			
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( cell_type[ Cell_i->type ] == PLASMA ) 
					Mobi[ 0 ][ i ] = config->Species[ 0 ].MobilityValue ;
			}
			break;

		case 1://cout<<"Mobility calculate from collision table." <<endl ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				collision = TotalNumberDensity[i]*CollTable.GetValue( T[ 0 ][ i ] ) ; //unit
				if( cell_type[ Cell_i->type ] == PLASMA ) 
					Mobi[ 0 ][ i ] = (Qe*Ref_Qe)/(Me*Ref_Mass)/collision ;
			}
			break;

		case 2://cout<<"Mobility calculate from electron temperature table." <<endl ;
			MobilityIter = MobilityMap.find( 0 ) ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( cell_type[ Cell_i->type ] == PLASMA ) 
				Mobi[ 0 ][ i ] = MobilityIter->second.GetValue( T[ 0 ][ i ] )/TotalNumberDensity[i] ;
			}
			break;

		case 3://cout<<"Mobility calculate from reduce electric field table." <<endl ;
			MobilityIter = MobilityMap.find( 0 ) ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
				Cell_i  = plasma.get_cell( i ) ;
				if( cell_type[ Cell_i->type ] == PLASMA ) 
					Mobi[ 0 ][ i ] = MobilityIter->second.GetValue( Etd[ i ] )/TotalNumberDensity[i] ;
			}
			break;
		case 4:
			MobilityIter = MobilityMap.find( 0 ) ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( cell_type[ Cell_i->type ] == PLASMA ) {
					Mobi[ 0 ][ i ] = MobilityIter->second.GetValueLog( Etd[ i ] ) ;
				}
			}
			break;

		case 5:// This is for N*De = value ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					Mobi[ 0 ][ i ] = config->Species[ 0 ].MobilityValue/TotalNumberDensity[i] ;
					//cout<<"NTot: "<<TotalNumberDensity[i]<<endl<<endl;
					//cout<<config->Species[ 0 ].MobilityValue<<endl;
				}
			}
			break;
		case 6:
			MobilityIter = MobilityMap.find( 0 ) ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					Mobi[ 0 ][ i ] = MobilityIter->second.GetValueLog( Emag[ i ] ) ;
				}
			}
			break;
		case 7:
			exit(1);
			break;

		case 8://Nishida
			/*
				Fit coefficients y=exp(A+B*ln(x)+C/x+D/x^2+E/x^3)
   				56.84     -0.2917       2.864      -5.498       3.265    
			*/
			
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if (plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

					E_Td = Etd[ i ] ;
					if ( E_Td < 1.0 ) E_Td = 1.0 ;
					else if ( E_Td > 3000.0 ) E_Td = 3000.0 ;
					ReduceMobi = exp( 56.84 + (-0.2917)*log(E_Td) + (2.864)/E_Td + (-5.498)/E_Td/E_Td + (3.265)/E_Td/E_Td/E_Td ) ;
					Mobi[ 0 ][ i ] = ReduceMobi/TotalNumberDensity[i] ;
				}
			}
			break;

		case 9://Ward: pμ = 3x10^5 [cm^2*Torr/V/S]

			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					P_torr = TotalGasPressure[ i ] ;
					Mobi[ 0 ][ i ] = 30.0/P_torr ;
				}
			}
			break;

		case 10://Air
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					Etd_Min = 2.4E+6/TotalNumberDensity[i]/1.E-21 ;
					E_Td = Etd[ i ] ;
					if( E_Td < Etd_Min ) E_Td = Etd_Min ; 
					DriftV = 4.2E+3*pow( E_Td, 0.74 ) ;
					Mobi[ 0 ][ i ] = DriftV/Emag[ i ] ;
				}

			}
			break;

		case 11://Bagheri et. al., “Comparison of six simulation codes for positive streamers in air,” Plasma Sources Science and Technology, vol. 27, no. 9, p. 095002, 2018.
					
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
				Cell_i  = plasma.get_cell( i ) ;
				if( cell_type[ Cell_i->type ] == PLASMA ) {
		 			Mobi[ 0 ][ i ] = 2.3987*pow( Emag[i], -0.26 ) ;
				}
			}
			break;

		case 12://P. Dordizadeh et al. , “Numerical investigation of the formation of Trichel pulses in a needle-plane geometry,” Journal of Physics D: Applied Physics 48(41),415203 (2015).	
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
				Cell_i  = plasma.get_cell( i ) ;
				if( cell_type[ Cell_i->type ] == PLASMA ) {
					E_Td = Etd[ i ] ;
		 			Mobi[ 0 ][ i ] = 5.747e16*pow(E_Td,0.6064)/Emag[ i ] ;
				}
			}
			break;


		default:
				if ( mpi_rank == 0 ){
					cout <<"The Mobility Type of species[0] you choose is illegal. Please check!!!!" << endl ;
					cout <<"0: constant mobility"<<endl;
					cout <<"1: mobility calculate form input \"CollisionFreqFile.inp\""<<endl;
				}
				exit(1);
			    break;
	}
	switch ( config->Species[ 0 ].DiffusivityType ){

		case 0:
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) 
					Diff[ 0 ][ i ] = config->Species[ 0 ].DiffusivityValue ;
			}
			break;
		case 1://Einstein relation
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) Diff[ 0 ][ i ] = T[ 0 ][ i ]*Mobi[ 0 ][ i ] ;
			}
			break;
		case 2://Diffusion calculate from electron temperature table.
				DiffusivityIter = DiffusivityMap.find( 0 ) ;
				for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
					Cell_i  = plasma.get_cell( i ) ;
					if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) Diff[ 0 ][ i ] = DiffusivityIter->second.GetValue( T[ 0 ][ i ] )/TotalNumberDensity[i] ;
					//cout<<"Diff: "<<Diff[ 0 ][ i ]<<endl;
				}
				//exit(1) ;
			break;
		case 3://Diffusion calculate from reduce electric field table.
				DiffusivityIter = DiffusivityMap.find( 0 ) ;
				for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
					Cell_i  = plasma.get_cell( i ) ;
					if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) Diff[ 0 ][ i ] = DiffusivityIter->second.GetValue( Etd[ i ] )/TotalNumberDensity[i] ;
				}
			break;
		case 4:
				DiffusivityIter = DiffusivityMap.find( 0 ) ;
				for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
					Cell_i  = plasma.get_cell( i ) ;
					if(plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" )Diff[ 0 ][ i ] = DiffusivityIter->second.GetValueLog( Etd[ i ] ) ;
				}
			break;
		case 5://this is for N*De = value ;
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					Diff[ 0 ][ i ] = config->Species[ 0 ].DiffusivityValue/TotalNumberDensity[i] ;
				}
			}
			break;
		case 6:
				DiffusivityIter = DiffusivityMap.find( 0 ) ;
				for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
					Cell_i  = plasma.get_cell( i ) ;
					if(plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" )
						Diff[ 0 ][ i ] = DiffusivityIter->second.GetValueLog( Emag[ i ] ) ;
				}
			break;
		case 7:
			exit(1);
			break;

		case 8://Nishida: Einstein relation
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if(plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){
					Diff[ 0 ][ i ] = T[ 0 ][ i ]*Mobi[ 0 ][ i ] ;
					//cout<<Diff[0][i]<<endl;
				}
			}
			break;

		case 9://Ward: Einstein relation
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if(plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){
					Diff[ 0 ][ i ] = T[ 0 ][ i ]*Mobi[ 0 ][ i ] ;
				}
			}
			break;

		case 10:
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				Cell_i  = plasma.get_cell( i ) ;
				if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					Etd_Min = 2.4E+6/TotalNumberDensity[i]/1.E-21 ;
					E_Td = Etd[ i ] ;
					if( E_Td < Etd_Min ) E_Td = Etd_Min ; 
					Diff[ 0 ][ i ] = 9.7E23*pow(E_Td,0.22)/TotalNumberDensity[i] ; 
				}
			}
			break;
		case 11://Bagheri et. al., “Comparison of six simulation codes for positive streamers in air,” Plasma Sources Science and Technology, vol. 27, no. 9, p. 095002, 2018.
      for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
        Cell_i  = plasma.get_cell( i ) ;
        if( cell_type[ Cell_i->type ] == PLASMA ) {
          Diff[ 0 ][ i ] = (4.3628E-3)*pow( Emag[i], 0.22 ) ;
        }
      }
		break;

		case 12://P. Dordizadeh et al. , “Numerical investigation of the formation of Trichel pulses in a needle-plane geometry,” Journal of Physics D: Applied Physics 48(41),415203 (2015).
      for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
        Cell_i  = plasma.get_cell( i ) ;
        if( cell_type[ Cell_i->type ] == PLASMA ) {

        	E_Td = Etd[ i ] ;

	     		if ( E_Td  < 2.0 ) {

						Diff[ 0 ][ i ] = 1.343e6*pow( E_Td*1.0e-21, 0.3441 )*Mobi[ 0 ][ i ] ;

	        }else if ( E_Td >= 2.0 and E_Td < 12.3 ) {

						Diff[ 0 ][ i ] = 1.213e25*pow( E_Td*1.0e-21, 1.2601 )*Mobi[ 0 ][ i ] ;

	        } else {

						Diff[ 0 ][ i ] = 1.519e9*pow( E_Td*1.0e-21, 0.46113 )*Mobi[ 0 ][ i ] ;

	        }
        }
      }
		break;
		default:
			if ( mpi_rank == 0 ){
				cout <<"The Diffusivity Type of species[0] you choose is illegal. Please check!!!!\n" << endl ;
				cout <<"0: constant Diffusivity"<<endl;
				cout <<"1: Diffusivity calculate form input \"CollisionFreqFile.inp\" \n"<<endl;
			}
			exit(1);
			break;
	}
	/*--- Normalize ---*/
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		Cell_i  = plasma.get_cell( i ) ;
		if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
			// cout<<"Mu: "<<Mobi[ 0 ][ i ]<<endl ;
			// cout<<" D: "<<Diff[ 0 ][ i ]<<endl ;
			Diff[ 0 ][ i ] /= Ref_Diff ;
			Mobi[ 0 ][ i ] /= Ref_Mu ;
		}
	}


	Mobi[ 0 ] = Mobi[ 0 ] ;
	Diff[ 0 ] = Diff[ 0 ] ;
	//exit(1) ;
}
void CVariable::UpdateIonNeutralTransport( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	Cell *Cell_i ;
	double EoverP=0.0, E=0.0, P_torr = 0.5, value = 0.0, collision=0.0, Thermal2=0.0, U=0.0, V=0.0, VTOT2=0.0 ;
	double X=0.0, C1=0.0, C2=0.0, Thermal=0.0 ;
	map< int, CTable>::iterator MobilityIter, DiffusivityIter ;

	for ( int iSpecies = 1 ; iSpecies < config->TotalSpeciesNum ; iSpecies ++ ) {

		if ( config->Species[ iSpecies ].Type == ION or config->Species[ iSpecies ].Type == NEUTRAL ){

			switch ( config->Species[ iSpecies ].MobilityType ){
				case 0://Constant Mobility Value
					if( config->Species[ iSpecies ].ConstantMobilityUpdate == false ) {
						for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
							Cell_i  = plasma.get_cell( i ) ;
							if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
								Mobi[ iSpecies ][ i ] = config->Species[ iSpecies ].MobilityValue ;
							}
						}
						//config->Species[ iSpecies ].ConstantMobilityUpdate = true ;
					}
					break;
				case 1:
					cout<<"For ion and neutral the transport coefficients is not support calculate by collision frequency!!!"<<endl;
					exit(1) ;
					break;
				case 2://
					MobilityIter = MobilityMap.find( iSpecies ) ;
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Mobi[ iSpecies ][ i ] = MobilityIter->second.GetValueLog( Emag[ i ] ) ;
						}
					}
					break;
				case 3:
				
					MobilityIter = MobilityMap.find( iSpecies ) ;
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
						 Mobi[ iSpecies ][ i ] = MobilityIter->second.GetValueLog( Etd[ i ] )/TotalNumberDensity[i] ;
						}
					}
					break;

				case 4://From Mobility table
					MobilityIter = MobilityMap.find( iSpecies ) ;
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Mobi[ iSpecies ][ i ] = MobilityIter->second.GetValueLog( Etd[ i ] ) ;
						}
					}
					break;
				case 5:
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Mobi[ iSpecies ][ i ] = config->Species[ iSpecies ].MobilityValue/TotalNumberDensity[i] ;
							//cout<<"NTot: "<<TotalNumberDensity[i]<<endl<<endl;
						}
					}
					break;

				case 6://Calculate mobility from reduce mobility file.
					MobilityIter = MobilityMap.find( iSpecies ) ;
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						value = MobilityIter->second.GetValueLog( Etd[ i ] ) ;
						if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Mobi[ iSpecies ][ i ] = value*(760.0/config->P_back)*(config->T_back/273.16) ;
						}
					}

					break;
				case 7:

					break;
				case 8:
					exit(1);
					break;
				case 9://Ward
					/* 
						Mobility is almost always specified in units of cm2/(V·s). 
						This is different from the SI unit of mobility, m2/(V·s). They are related by 1m2/(V·s) = 10^4cm2/(V·s).
					*/
	          for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
	          	Cell_i  = plasma.get_cell( i ) ;

	            if (plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
	             		P_torr = TotalGasPressure[ i ] ;
		             	E = sqrt ( Ex[i]*Ex[i] + Ey[i]*Ey[i] + Ez[i]*Ez[i] )*0.01 ;
		              EoverP = E/P_torr ; // unit: cm
		          	if ( E/P_torr <=  60.0 ){
									Mobi[ iSpecies ][ i ] = 1.0E3*(1.0-2.22E-3*EoverP )/P_torr/1.0E4  ;
								} else {
									Mobi[ iSpecies ][ i ] = 8.25E3/sqrt(EoverP)*(1.0 - 86.52/pow(EoverP, 1.5) )/P_torr/1.0E4  ;
								}
								if( Mobi[ iSpecies ][ i ] < 0.0 ) cout<<"Mobility Error"<<endl;
	          	}
						}
					break;
				case 10:
				break;

				case 11://JPL paper
					
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
						Cell_i  = plasma.get_cell( i ) ;
						if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							C1 = 6.6E-19*( 0.25*T[0][i] - 0.1 ) ;
							//C2 = 1.0+pow( 0.25*T[0][i], 1.6 ) ;
							Thermal = sqrt( 2.0*Qe*Ref_Qe*T[1][i]/PI/config->Species[ 1 ].Mass_Kg ) ;
							X = sqrt( U1[1][i]*U1[1][i] + U2[1][i]*U2[1][i] )/U0[1][i]/Thermal + ZERO;
							//cout<<X<<endl;
							C1 = exp(-X*X) + ( 2.0*X + 1.0/X )*0.8862269255*erf(X) ;
							//cout<<C1<<endl;
					 		collision = TotalNumberDensity[i]*Thermal*(1.E-18)*C1 ;
				 			Mobi[ iSpecies ][ i ] = (Qe*Ref_Qe)/(Me*Ref_Mass)/collision ;
						}
					}
					//exit(1) ;
				break;
				case 12:
					
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){

						Cell_i  = plasma.get_cell( i ) ;

						if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) 
						{
				 			Mobi[ iSpecies ][ i ] = 5.5E21/TotalNumberDensity[i] ;
						}
					}
					//exit(1) ;
				break;
				default:
						if ( mpi_rank == 0 ){
							cout <<"The Mobility Type of species[0] you choose is illegal. Please check!!!!" << endl ;
							cout <<"0: constant mobility"<<endl;
							cout <<"1: mobility calculate form input \"CollisionFreqFile.inp\""<<endl;
							cout<<"Current input- iSpecies["<<iSpecies<<"]: "<<config->Species[ iSpecies ].MobilityType <<endl;
						}
						exit(1);
					    break;
			}
			switch ( config->Species[ iSpecies ].DiffusivityType ){

				case 0:
					if( config->Species[ iSpecies ].ConstantDiffusivityUpdate == false ){
						for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
							Cell_i  = plasma.get_cell( i ) ;
							if (plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
								Diff[ iSpecies ][ i ] = config->Species[ iSpecies ].DiffusivityValue ;
							}
						}
						//config->Species[ iSpecies ].ConstantDiffusivityUpdate = true ; 
					}
					break;
				case 1:
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Diff[ iSpecies ][ i ] = Mobi[ iSpecies ][ i ]*T[iSpecies][i] ;
						}
					}
					break;
				case 2:
					exit(1);
					break;
				case 3:
					DiffusivityIter = DiffusivityMap.find( iSpecies ) ;
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Diff[ iSpecies ][ i ] = DiffusivityIter->second.GetValueLog( T[ 0 ][ i ] )/TotalNumberDensity[i] ;
						}
					}
					break;
				case 4:
					DiffusivityIter = DiffusivityMap.find( iSpecies ) ;
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Diff[ iSpecies ][ i ] = DiffusivityIter->second.GetValueLog( Etd[ i ] ) ;
						}
					}
					break;
				case 5:
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
							Diff[ iSpecies ][ i ] = config->Species[ iSpecies ].DiffusivityValue/TotalNumberDensity[i] ;
						}
					}
					break;
				case 6:
					exit(1) ;
					break;
				case 7:
					exit(1);
					break;
				case 8:
					exit(1);
					break;
				case 9:
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {//@ 1Torr, 300 k
							P_torr = TotalGasPressure[ i ] ;
							Diff[ iSpecies ][ i ] = (2.0E2/P_torr)*(1.E-4) ;
						}
					}
					break;
				case 10:
					for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
						Cell_i  = plasma.get_cell( i ) ;
						if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {//@ 1Torr, 300 k
							P_torr = TotalGasPressure[ i ] ;
							//cout<<P_torr<<endl;
							Diff[ iSpecies ][ i ] = (82.992/P_torr)*(1.E-4) ;
						}
					}
					//cout << " The mobility_type of e, It does not be prepared.  " << endl;
					break;
				default:
					if ( mpi_rank == 0 ){
						cout <<"The Diffusivity Type of species[0] you choose is illegal. Please check!!!!\n" << endl ;
						cout <<"0: constant mobility"<<endl;
						cout <<"1: mobility calculate form input \"CollisionFreqFile.inp\" \n"<<endl;
					}
					exit(1);
					break;
			}

			/*--- Normalize ---*/
			for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
				if (plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {
					// cout<<"Mu: "<<Mobi[ iSpecies ][ i ]<<endl ;
					// cout<<" D: "<<Diff[ iSpecies ][ i ]<<endl ;
					Diff[ iSpecies ][ i ] /= Ref_Diff ;
					Mobi[ iSpecies ][ i ] /= Ref_Mu ;
				}
			}
			// exit(1);
			Mobi[ iSpecies ] = Mobi[ iSpecies ] ;
			Diff[ iSpecies ] = Diff[ iSpecies ] ;
		}
	}
}
void CVariable::ResetAvgZero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	eAvgEnergyLoss.zero() ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) 
	{
		AvgPotential[ i ] = 0.0 ;
		       AvgEx[ i ] = 0.0 ;
		       AvgEy[ i ] = 0.0 ;
		       AvgEz[ i ] = 0.0 ;
		     AvgE_Td[ i ] = 0.0 ;
		avg_plot_var[i] = 0.0 ;
	}

	if ( plasma.Mesh.ndim == 3 ) {	
		for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
			     AvgEz[ i ] = 0.0 ;
		}
	}

	for ( int iSpecies = 0 ; iSpecies < config->TotalSpeciesNum ; iSpecies++ ){
		AvgT[ iSpecies ].zero() ;
		AvgU0[ iSpecies ].zero() ;
		AvgU1[ iSpecies ].zero() ;
		AvgU2[ iSpecies ].zero() ;
		AvgU3[ iSpecies ].zero() ;
		AvgU4[ iSpecies ].zero() ;
		AvgJouleHeating[ iSpecies ].zero() ;
	}
	//AvgPowerAbs = 0.0 ;
}
void CVariable::AddAverage( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ )  {
		AvgPotential[ i ] += Potential.current[ i ] / config->StepPerCycle ;
		AvgEx[ i ]        +=  Ex[ i ] / config->StepPerCycle ;
		AvgEy[ i ]        +=  Ey[ i ] / config->StepPerCycle ;
		AvgEz[ i ]        +=  Ez[ i ] / config->StepPerCycle ;
		AvgE_Td[ i ]      += Etd[ i ] / config->StepPerCycle ;
		eAvgEnergyLoss[ i ] += eEnergyLoss[ i ]/config->StepPerCycle ;
		avg_plot_var[ i ] += avg_plot_var[ i ]/config->StepPerCycle ;
	}
	for ( int iSpecies = 0 ; iSpecies < config->TotalSpeciesNum ; iSpecies++ ){
		for ( int i = 0 ; i <plasma.Mesh.cell_number  ; i++ ) {
			 AvgT[ iSpecies ][ i ] +=  T[ iSpecies ][ i ]/config->StepPerCycle ;
			AvgU0[ iSpecies ][ i ] += U0[ iSpecies ][ i ]/config->StepPerCycle ;
			AvgU1[ iSpecies ][ i ] += U1[ iSpecies ][ i ]/config->StepPerCycle ;
			AvgU2[ iSpecies ][ i ] += U2[ iSpecies ][ i ]/config->StepPerCycle ;
			AvgU3[ iSpecies ][ i ] += U3[ iSpecies ][ i ]/config->StepPerCycle ;
			AvgU4[ iSpecies ][ i ] += U4[ iSpecies ][ i ]/config->StepPerCycle ;
			AvgJouleHeating[ iSpecies ][ i ] += JouleHeating[ iSpecies ][ i ]/config->StepPerCycle ;
		}
	}
}
void CVariable::ResetAvgZero_PowerAbs( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	AvgPowerAbs = 0.0 ;
}
void CVariable::AddAverage_PowerAbs( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	AvgPowerAbs += PowerAbs/config->StepPerCycle  ;
}
void CVariable::AddAverage_Electrode( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	I_AvgPowerElectrode 	+= I_PowerElectrode_global_sum/config->StepPerCycle ;
	I_AvgGroundElectrode 	+= I_GroundElectrode_global_sum/config->StepPerCycle ;
}
void CVariable::ResetAvgZero_Electrode( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	I_AvgPowerElectrode=0.0 ;
	I_AvgGroundElectrode=0.0 ;
}
void CVariable::CalculateElectrodeCurrent( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	I_PowerElectrode_local 		= 0.0 ;	I_GroundElectrode_local 		= 0.0 ;
	I_PowerElectrode_global_sum = 0.0 ;	I_GroundElectrode_global_sum 	= 0.0 ;

	DispI_PowerElectrode_local	= 0.0 ;	DispI_GroundElectrode_local 	= 0.0 ;
	Disp_PowerElectrode_global_sum	= 0.0 ;	Disp_GroundElectrode_global_sum 	= 0.0 ;

	for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ){

    	CondI_PowerElectrode_local[iSpecies] = 0.0 ;
	    CondI_PowerElectrode_global_sum[iSpecies] = 0.0 ;

    	CondI_GroundElectrode_local[iSpecies] = 0.0 ;	
	    CondI_GroundElectrode_global_sum[iSpecies] = 0.0 ;
	}


	int j=0 ;

	Cell *Cell_i ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		Cell_i->face_number 	 = Cell_i->face_number ;
	
		/*--- Loop over electrode cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" or  plasma.get_cell_typename( Cell_i->data_id ) == "DIELECTRIC" ){

			/*--- Loop over bulk faces ---*/
			for( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->data_id ; 

				if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" ){

					I_PowerElectrode_local  += ( TotalJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[ 0 ]
										 	+ 	 TotalJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[ 1 ] ) ;

					DispI_PowerElectrode_local += ( DispJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 		+   DispJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] );

					for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ){
						CondI_PowerElectrode_local[iSpecies]  += (-1.0) * (CondJD[iSpecies][0][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 			 	 				+  CondJD[iSpecies][1][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;
					}

				} else if ( plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ){

					I_GroundElectrode_local += ( TotalJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 	+  	 TotalJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;

					DispI_GroundElectrode_local += ( DispJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 		+    DispJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;

					for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ){
						CondI_GroundElectrode_local[iSpecies]  += ( CondJD[iSpecies][0][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 			 	 		+   CondJD[iSpecies][1][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;
					}

				}
			}

			/*--- Loop over boundary faces ---*/
			for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {


				if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id) == "POWER" ){

					I_PowerElectrode_local  += ( TotalJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 	             +   TotalJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;

					DispI_PowerElectrode_local += ( DispJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 		+   DispJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;
					
					for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ){
						CondI_PowerElectrode_local[iSpecies]  += ( CondJD[iSpecies][ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 					   +   CondJD[iSpecies][ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;
					}

				}else if( plasma.get_face_typename( Cell_i->face[ k ]->data_id) == "GROUND" ){

					I_GroundElectrode_local += ( TotalJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 		+TotalJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;

					DispI_GroundElectrode_local += ( DispJD[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 		+    DispJD[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;

					for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ){
						CondI_GroundElectrode_local[iSpecies]  += ( CondJD[iSpecies][0][ i ]*m->PFM_CELL[ i ][ k ].Af[0]
										 			 	 			   	+  CondJD[iSpecies][1][ i ]*m->PFM_CELL[ i ][ k ].Af[1] ) ;
					}

				}
			}
		} 
	}//Cell Loop

	MPI_Allreduce(&I_PowerElectrode_local, &I_PowerElectrode_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&I_GroundElectrode_local, &I_GroundElectrode_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(&DispI_PowerElectrode_local, &Disp_PowerElectrode_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&DispI_GroundElectrode_local, &Disp_GroundElectrode_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for ( int iSpecies = 0 ; iSpecies < config->ChargeSpeciesNum ; iSpecies++ ){
		MPI_Allreduce(&CondI_PowerElectrode_local[iSpecies], &CondI_PowerElectrode_global_sum[iSpecies], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&CondI_GroundElectrode_local[iSpecies], &CondI_GroundElectrode_global_sum[iSpecies], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	}
	// cout<<"I_PowerElectrode_global_sum : "<<I_PowerElectrode_global_sum<<endl;
	// cout<<"I_GroundElectrode_global_sum: "<<I_GroundElectrode_global_sum<<endl;
	// cout<<endl;
}
void CVariable::Alpha_Beta( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{

	double Alpha=0.0, T_Ratio=0.0 ;
	double factor = 1.0 ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell *Cell_i  = plasma.get_cell( i ) ;

		if ( cell_type[ Cell_i->type ] != PLASMA ) {

			Beta[ i ] = 0.0 ;		

		} else {

			Alpha 	= 0.0 ;
			T_Ratio = 0.0 ;

			for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

				if ( config->Species[ jSpecies ].Charge < 0 ) {
					Alpha 	+= ( U0[ jSpecies ][ i ]/U0[ 0 ][ i ] ) ;
					T_Ratio  = ( U0[ jSpecies ][ i ]/U0[ 0 ][ i ])*( T[ 0 ][ i ]/T[ jSpecies ][ i ] ) ;
				}//For negative ion

			}//End jSpecies

			Beta[ i ] = factor*(1.0+Alpha)/(1.0+T_Ratio) ;

		}
	}//Loop over all cells

	Beta = Beta ;
}
// void CVariable::Calculate_LSQ_Coeff_Scalar( boost::shared_ptr<CDomain> &m )
// {
// 	// if ( mpi_rank==MASTER_NODE ) cout<<" Staring build least-square gradient coefficients ..."<<endl;

// 	// for ( int k = 0 ; k < 6 ; k++ ) {
// 	// 	LSQ_Cx[ k ].initial( "LSQ_Cx" ) ;			
// 	// 	LSQ_Cy[ k ].initial( "LSQ_Cy" ) ;			
// 	// 	LSQ_Cz[ k ].initial( "LSQ_Cz" ) ;			
// 	// } 

// 	// double dx=0.0, dy=0.0, a11=0.0, a12=0.0, a21=0.0, a22=0.0, det=0.0 ;
// 	// double ia11=0.0, ia12=0.0, ia21=0.0, ia22=0.0 ;
// 	// Cell *Cell_i, *Cell_j ;

// 	// for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

// 	// 	Cell_i  = plasma.get_cell( i ) ;
		

// 	// 	Cell_i->face_number = Cell_i->face_number ;

// 	// 	/*--- Reset Matrix ---*/
// 	// 	a11 = 0.0 ; a12 = 0.0 ;
// 	// 	a21 = 0.0 ; a22 = 0.0 ;

// 	// 	/*--- Loop over neighbor "cells" ---*/
// 	// 	for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

// 	// 		if ( plasma.get_cell_typename( Cell_i->data_id ) != plasma.get_cell_typename( Cell_i->cell[ k ]->data_id ) ) {//For discontinued face

// 	// 			Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
// 	// 			Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

// 	// 			fPf[ 0 ] = DotProduct( Pf, m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[0] ;
// 	// 			fPf[ 1 ] = DotProduct( Pf, m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[1] ;
// 	// 			//dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
// 	// 			//dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;
// 	// 			//
// 	// 			dx = Cell_i->face[ k ]->r[0]  - Cell_i->r[0] ;
// 	// 			dy = Cell_i->face[ k ]->r[1]  - Cell_i->r[1] ;

// 	// 		} else {

// 	// 			dx =Cell_i->cell[ k ]->r[0]  - Cell_i->r[0] ;
// 	// 			dy =Cell_i->cell[ k ]->r[1]  - Cell_i->r[1] ;
// 	// 		}

// 	// 		a11 = a11 + dx*dx ; 
// 	// 		a12 = a12 + dx*dy ;
// 	// 		a21 = a21 + dx*dy ; 
// 	// 		a22 = a22 + dy*dy ;

// 	// 	}

// 	// 	/*--- Loop over domain boundary "faces" ---*/
// 	// 	for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

// 	// 		Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
// 	// 		Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

// 	// 		fPf[ 0 ] = DotProduct( Pf, m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[0] ;
// 	// 		fPf[ 1 ] = DotProduct( Pf, m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[1] ;

// 	// 		//dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
// 	// 		//dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;
// 	// 		//cout<<fPf[ 0 ]<<"\t"<<fPf[ 0 ]<<endl;
// 	// 		dx = Cell_i->face[ k ]->r[0]  - Cell_i->r[0] ;
// 	// 		dy = Cell_i->face[ k ]->r[1]  - Cell_i->r[1] ;

// 	// 		a11 = a11 + dx*dx ; 	
// 	// 		a12 = a12 + dx*dy ;
// 	// 		a21 = a21 + dx*dy ;	
// 	// 		a22 = a22 + dy*dy ;
// 	// 	}

// 	// 	/*--- Cal. LSQ det. ---*/
// 	// 	det = a11*a22 - a12*a21 ;
// 	// 	//if( fabs(det) < 1.E-14 ) cout<<"LSQ det Singular"<<endl;

// 	// 	--- Cal. Inverse Matrix ---
// 	// 	ia11 =  a22/det ;
// 	// 	ia12 = -a21/det ;
// 	// 	ia21 = -a12/det ;
// 	// 	ia22 =  a11/det ;

// 	// 	/*--- Cal. LSQ Coefficient for "cells" ---*/
// 	// 	for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

// 	// 		if ( plasma.get_cell_typename( Cell_i->data_id ) != plasma.get_cell_typename( Cell_i->cell[ k ]->data_id ) ) {

// 	// 			Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
// 	// 			Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
// 	// 			fPf[ 0 ] = DotProduct( Pf,  m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[0] ;
// 	// 			fPf[ 1 ] = DotProduct( Pf,  m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[1] ;
// 	// 			//dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
// 	// 			//dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;
// 	// 			//cout<<fPf[ 0 ]<<"\t"<<fPf[ 0 ]<<endl;
// 	// 			dx =  Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
// 	// 			dy =  Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

// 	// 		} else {

// 	// 			dx = Cell_i->cell[ k ]->r[0]  - Cell_i->r[0] ;
// 	// 			dy = Cell_i->cell[ k ]->r[1]  - Cell_i->r[1] ;
// 	// 			//cout<<"normal dy: "<<dy<<endl;
// 	// 		}
// 	// 		LSQ_Cx[ k ][ i ] = ia11*dx + ia12*dy ;
// 	// 		LSQ_Cy[ k ][ i ] = ia21*dx + ia22*dy ;
// 	// 		//cout<<"LSQ_Cx["<<k<<"]["<<i<<"]: "<<LSQ_Cx[ k ][ i ]<<endl;
// 	// 		//cout<<"LSQ_Cy["<<k<<"]["<<i<<"]: "<<LSQ_Cy[ k ][ i ]<<endl;
// 	// 	}
// 	// 	//if(mpi_rank==0) cout<<"---------------------------------------"<<endl;
// 	// 	/*--- Cal. LSQ Coefficient for domain boundary "faces" ---*/
// 	// 	for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

// 	// 		Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
// 	// 		Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
// 	// 		fPf[ 0 ] = DotProduct( Pf, m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[0] ;
// 	// 		fPf[ 1 ] = DotProduct( Pf, m->PFM_CELL[ i ][ k ].mf )*m->PFM_CELL[ i ][ k ].mf[1] ;
// 	// 		//dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
// 	// 		//dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

// 	// 		dx =  Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
// 	// 		dy =  Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
// 	// 		LSQ_Cx[ k ][ i ] = ia11*dx + ia12*dy ;
// 	// 		LSQ_Cy[ k ][ i ] = ia21*dx + ia22*dy ;
// 	// 		//cout<<"LSQ_Cx["<<k<<"]["<<i<<"]: "<<LSQ_Cx[ k ][ i ]<<endl;
// 	// 		//cout<<"LSQ_Cy["<<k<<"]["<<i<<"]: "<<LSQ_Cy[ k ][ i ]<<endl;
// 	// 	}//End boundaty face
// 	// 	//cout<<endl;
// 	// }//End cell loop
// 	// //MPI_Barrier(MPI_COMM_WORLD); exit(1) ;
// 	// for( int k = 0 ;  k < 6 ; k++ ){
// 	// 	LSQ_Cx[ k ] = LSQ_Cx[ k ] ;
// 	// 	LSQ_Cy[ k ] = LSQ_Cy[ k ] ;
// 	// 	LSQ_Cz[ k ] = LSQ_Cz[ k ] ;
// 	// }
// 	// if ( mpi_rank==MASTER_NODE ) cout<<" Build least-square gradient coefficients finish ..."<<endl;
// 	// //exit(1) ;
// }
void CVariable::SourceSink_PSST_2018( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
//Bagheri et. al., “Comparison of six simulation codes for positive streamers in air,” Plasma Sources Science and Technology, vol. 27, no. 9, p. 095002, 2018.
	double alpha=0.0, alpha_eta=0.0, eta=0.0 ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number  ; i++ ) {

		Cell *Cell_i  = plasma.get_cell( i ) ;

		if ( cell_type[ Cell_i->type ] == PLASMA ) {
		
			alpha = (1.1944E6 + 4.3666E26/pow(Emag[i],3.0) )*exp( (-2.73E7)/Emag[i] ) ;
			eta   = 340.75 ;
			alpha_eta = alpha - eta ;

			*( ReactionRatePoint[ 0 ] + i  ) = (alpha_eta)*Mobi[ 0 ][ i ]*Emag[ i ]*U0[ 0 ][ i ] ;//+ R_stable ;
			*( ReactionRatePoint[ 1 ] + i  ) = (alpha_eta)*Mobi[ 0 ][ i ]*Emag[ i ]*U0[ 0 ][ i ] ;//+ R_stable ;
			LFASourceSink[ 0 ][ i ] 		     = (alpha_eta)*Mobi[ 0 ][ i ]*Emag[ i ]*U0[ 0 ][ i ] ;//+ R_stable ;
			LFASourceSink[ 1 ][ i ] 		     = (alpha_eta)*Mobi[ 0 ][ i ]*Emag[ i ]*U0[ 0 ][ i ] ;//+ R_stable ;
		}
	}
	LFASourceSink[ 0 ] = LFASourceSink[ 0 ] ;
	LFASourceSink[ 1 ] = LFASourceSink[ 1 ] ;
}
void CVariable::SourceSink_2Fluid( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	double alpha=0.0, eta=0.0, Gamma_ep=2.E-13, Gamma_np=2.E-13, Flux2=0.0, E_Td=0.0 ;
	double nu=0.0, nv=0.0 ;
	double U0_e=0.0, U0_pi=0.0, U0_ni=0.0 ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number  ; i++ ) {

		Cell *Cell_i  = plasma.get_cell( i ) ;

		if ( cell_type[ Cell_i->type ] == PLASMA ) {

			E_Td = Etd[ i ] ;

			alpha 	= AlphaTable.GetValue( Emag[ i ] )* TotalNumberDensity[ i ] ;
			cout<<"alpha: "<<alpha<<endl;
			nu = U1[0][ i ] ;
			nv = U2[0][ i ] ;

			U0_e = U0[ 0 ][ i ]*Ref_N ;
			U0_pi= U0[ 1 ][ i ]*Ref_N ;


			Flux2 = sqrt( nu*nu + nv*nv ) ;
			cout<<"Flux: "<<Flux2<<endl;

			*( ReactionRatePoint[ 0 ] + i  ) = 	alpha*Flux2 ;//- Gamma_ep*U0_e*U0_pi ;
			*( ReactionRatePoint[ 1 ] + i  ) =  alpha*Flux2 ;//- Gamma_ep*U0_e*U0_pi ;//- Gamma_np*U0_pi*U0_ni ;

			LFASourceSink[ 0 ][ i ] = (alpha-eta)*Flux2 ;//- Gamma_ep*U0_e*U0_pi ;
			LFASourceSink[ 1 ][ i ] = (alpha    )*Flux2 ;//- Gamma_ep*U0_e*U0_pi ;//- Gamma_np*U0_pi*U0_ni ;

		}
	}

}
void CVariable::SourceSink_3Fluid( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	double alpha=0.0, eta=0.0, Gamma_ep=2.E-13, Gamma_np=2.E-13, Flux2=0.0, E_Td=0.0 ;
	double nu=0.0, nv=0.0 ;
	double U0_e=0.0, U0_pi=0.0, U0_ni=0.0 ;


	for ( int i = 0 ; i < plasma.Mesh.cell_number  ; i++ ) {

		Cell *Cell_i  = plasma.get_cell( i ) ;

		if ( cell_type[ Cell_i->type ] == PLASMA ) {

			E_Td = Etd[ i ] ;
			alpha 	= AlphaTable.GetValue( Emag[ i ] );// * TotalNumberDensity[ i ] ;
			eta 	=   EtaTable.GetValue( Emag[ i ] );// * TotalNumberDensity[ i ] ;

			//nu = config->Species[0].Charge * EField[0][ i ] * Mobi[0][ i ]*U0[0][ i ] ;
			//nv = config->Species[0].Charge * EField[1][ i ] * Mobi[0][ i ]*U0[0][ i ] ;
			nu = U1[0][ i ] ;
			nv = U2[0][ i ] ;

			U0_e = U0[ 0 ][ i ]*Ref_N ;
			U0_pi= U0[ 1 ][ i ]*Ref_N ;
			U0_ni= U0[ 2 ][ i ]*Ref_N ;

			Flux2 = sqrt( nu*nu + nv*nv ) ;

			*( ReactionRatePoint[ 0 ] + i  ) =  (alpha-eta)*Flux2 - Gamma_ep*U0_e*U0_pi ;
			*( ReactionRatePoint[ 1 ] + i  ) =  (alpha    )*Flux2 - Gamma_ep*U0_e*U0_pi - Gamma_np*U0_pi*U0_ni ;
			*( ReactionRatePoint[ 2 ] + i  ) =  (eta      )*Flux2 				        - Gamma_np*U0_pi*U0_ni ;

			LFASourceSink[ 0 ][ i ] = (alpha-eta)*Flux2 - Gamma_ep*U0_e*U0_pi ;
			LFASourceSink[ 1 ][ i ] = (alpha    )*Flux2 - Gamma_ep*U0_e*U0_pi - Gamma_np*U0_pi*U0_ni ;
			LFASourceSink[ 2 ][ i ] = (eta      )*Flux2 					  - Gamma_np*U0_pi*U0_ni ;

		}
	}

}
void CVariable::UltraMPPComputeEnergyLossFromTable( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	double eLoss = 0.0, eTemp = 0.0 ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number  ; i++ ) {

		Cell *Cell_i  = plasma.get_cell( i ) ;

		if ( cell_type[ Cell_i->type ] == PLASMA ) {

			eTemp = T[ 0 ][ i ] ;
			eLoss 	= eEnergyLossTable.GetValueLog( eTemp ) * TotalNumberDensity[ i ]/(Qe*Ref_Qe) ;
			//cout<<eLoss<<endl;
			*( EnergySourcePoint + i  ) =  eLoss ;
		}
	}

}
void CVariable::UltraMPPComputeArgonIonTemperaturePHELPS( int iSpecies )
{

  for ( int i = 0 ; i < plasma.Mesh.cell_number  ; i++ ) {

    Cell *Cell_i  = plasma.get_cell( i ) ;

    if ( cell_type[ Cell_i->type ] == PLASMA ) {

      T[iSpecies][ i ]  = 1.9*pow( Etd[ i ]/1000.0, 1.1) + 0.026 ;

    }//plasma cell

  }//cell loop
 // plasma.syn_parallel_cell_data( T[ iSpecies ].data_id );
	T[iSpecies]=T[iSpecies] ;
}

void CVariable::SourceSink_Streamer( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	double alpha=0.0, eta=0.0 ;
	double U0_e=0.0, E_Td=0.0 ;
	double T0=273.0, T=300.0 ;
	double Etd_Min = 0.0;
	double C0=750.0*T0/T, C1=1.75E+3, C2=1.15E12, C3=-4.0E4,C4=T/T0, C5=(-1.0)*750.0*T0/587.0/T, C6=0.0, psi=0.0 ;
	Cell *Cell_i ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number  ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

			U0_e = U0[0][i] ;
			alpha 	= AlphaTable.GetValue( Emag[ i ]  ) ;
			*( ReactionRatePoint[ 0 ] + i  ) = (alpha)*Mobi[ 0 ][ i ]*Emag[ i ] *U0_e ;//+ R_stable ;
			*( ReactionRatePoint[ 1 ] + i  ) = (alpha)*Mobi[ 0 ][ i ]*Emag[ i ] *U0_e ;//+ R_stable ;
			LFASourceSink[ 0 ][ i ] 		 = (alpha)*Mobi[ 0 ][ i ]*Emag[ i ] *U0_e ;//+ R_stable ;
			LFASourceSink[ 1 ][ i ] 		 = (alpha)*Mobi[ 0 ][ i ]*Emag[ i ] *U0_e ;//+ R_stable ;
		}
	}
	LFASourceSink[ 0 ] = LFASourceSink[ 0 ] ;
	LFASourceSink[ 1 ] = LFASourceSink[ 1 ] ;
}

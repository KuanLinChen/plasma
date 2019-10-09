
#include "solver_fluid_model.hpp"

using namespace std ;
CFluidModel::CFluidModel()
{
}
void CFluidModel::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, int index )
{
	iSpecies = index ;
	SpeciesType = config->Species[ index ].Type ;

	/*--- PETSc Solver ---*/	
		Correction = config->Equation[ SpeciesType ].Correction ;
		if ( mpi_rank == 0 ){
			cout<<"Creat "<<config->Species[iSpecies].Name<<" fluid equation (w/o drift-diffusion), index: "<<index<<", charge: "<<config->Species[ index ].Charge<<", Speciec Type: "<<SpeciesType<<endl ;
			cout<<"Correction: "<<Correction<<endl;
			cout<<"Gamma: "<<config->Species[ iSpecies ].Gamma<<endl;
		} 
		fixTe = false ;
		if ( plasma.Mesh.geometry == "cylindrical_y" ) {
			Omaga = 1.0 ;
		} else {
			Omaga = 0.0 ;
		}

	/*--- Allocate Residue ---*/
		Res = new double*[ 5 ] ;/* N, Ux, Uy, Uz, E */
		for ( int iEqn = 0 ; iEqn < 5 ; iEqn ++ ) {
			Res[iEqn] = new double [ plasma.Mesh.cell_number ] ;
		}

	/*--- Momentum Transfer ---*/
		CollisionFreq = new CScalar [ config->TotalSpeciesNum ] ;
		for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
			CollisionFreq[ jSpecies ].initial( "nu_"+config->Species[ jSpecies ].Name ) ;		
		}

	/*--- Momentum Transfer ---*/
		CollisionIntegral = new double [ plasma.Mesh.cell_number ] ;
		Mx = new double [ plasma.Mesh.cell_number ] ;
		My = new double [ plasma.Mesh.cell_number ] ;
		Mz = new double [ plasma.Mesh.cell_number ] ;
		Pressure.initial( "P_"+config->Species[ iSpecies ].Name ) ;	
		Thermal2.initial( "Thermal2_"+config->Species[ iSpecies ].Name ) ;	

		Grave = false ;
		CrossSection = 50E-20 ;
}
void CFluidModel::Solve_Continuity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{
	if ( config->Equation[ SpeciesType ].Equation == 1 ){
		CalculateTotalEnergy( m, config, variable ) ;
	}

	CalculateTemperature( m, config, variable ) ;

	CalculateThermal2( m, config, variable ) ;

	CalculateIonNeutralCollisionFrequency( m, config, variable ) ;
	//CalculateCollisionIntegral( m, config, variable ) ;
	if ( config->Equation[ SpeciesType ].Equation == 2 ) {
		CalculateKappa( m, config, variable ) ;
	}
	ComputeFlux_HLL( m, config, variable ) ;

	ContinuityIntegral( m, config, variable ) ;
}
void CFluidModel::Solve_Momentum( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{
	MomentumIntegral( m, config, variable ) ;

	if ( config->Equation[ SpeciesType ].Equation == 2 ){
		//Calculate_Gradient_T( m, variable ) ;
		//CalculateIonNeutralCollisionFrequency( m, config, variable ) ;
		//CalculateCollisionIntegral( m, config, variable ) ;
		EnergyDensityIntegral( m, config, variable ) ;
	}
	CalculateCondCurrentDensity( m, config, variable ) ;
}
void CFluidModel::ComputeFlux_HLL( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;
	double vn=0.0, Te=0.0 ;

	double RhoL=0.0, RhoR=0.0  ;
	double	 uL=0.0,   uR=0.0 ; double	 vL=0.0,   vR=0.0 ; double wL=0.0, wR=0.0 ; double  unL=0.0,  unR=0.0 ;
	double   pL=0.0,   pR=0.0 ; double   aL=0.0,   aR=0.0 ; double HL=0.0, HR=0.0 ;
	double  TiL=0.0,  TiR=0.0 ; double  TeL=0.0,  TeR=0.0 ;

	double F0 =0.0, F1 =0.0, F2 =0.0, F3 =0.0, F4 =0.0 ; /* iFlux @ interface */
	double F0L=0.0, F1L=0.0, F2L=0.0, F3L=0.0, F4L=0.0 ; /* iFlux @ left cell */
	double F0R=0.0, F1R=0.0, F2R=0.0, F3R=0.0, F4R=0.0 ; /* iFlux @ Reft cell */
	double  SL=0.0, SR=0.0, SLm=0.0, SRp=0.0 ;
	double LogicalSwitch = 0.0, En = 0.0, GradUvel=0.0, GradU4=0.0 ;
	double IonMass = config->Species[ iSpecies ].Mass_Kg / var->Ref_Mass ;

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {


		Cell_i = plasma.get_cell(i) ;

		for( int k = 0 ; k < 5 ; k++ ) Res[ k ][ i ] = 0.0 ;

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*Cell_i->r[0] ; 

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = plasma.get_cell(j) ;

				if ( plasma.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					/*--- Left state ---*/
					RhoL = var->U0[iSpecies][ i ] ;
					RhoR = var->U0[iSpecies][ j ] ;

					uL 	 = var->U1[iSpecies][ i ]/RhoL ;
					uR 	 = var->U1[iSpecies][ j ]/RhoR ;

					vL 	 = var->U2[iSpecies][ i ]/RhoL ;
					vR 	 = var->U2[iSpecies][ j ]/RhoR ;

					wL 	 = var->U3[iSpecies][ i ]/RhoL ;
					wR 	 = var->U3[iSpecies][ j ]/RhoR ;

					unL  = uL*m->PFM_CELL[ i ][ k ].nf[ 0 ] + vL*m->PFM_CELL[ i ][ k ].nf[ 1 ] + wL*m->PFM_CELL[ i ][ k ].nf[ 2 ] ;
					unR  = uR*m->PFM_CELL[ i ][ k ].nf[ 0 ] + vR*m->PFM_CELL[ i ][ k ].nf[ 1 ] + wR*m->PFM_CELL[ i ][ k ].nf[ 2 ] ;

					pL 	= Pressure[ i ] ;
					pR  = Pressure[ j ] ;

					HL  = var->U4[iSpecies][ i ] +	pL   ;
					HR  = var->U4[iSpecies][ j ] +	pR   ;

					TiL  = var-> T[iSpecies][ i ] ;
					TiR  = var-> T[iSpecies][ j ] ;

					TeL  = var-> T[ 0 ][ i ] ;
					TeR  = var-> T[ 0 ][ j ] ;

					aL  = sqrt( var->Qe*( TiL + var->Beta[ i ]*TeL )/IonMass ) ;
					aR  = sqrt( var->Qe*( TiR + var->Beta[ j ]*TeR )/IonMass ) ;

					/*---  Davis' wave estimate---*/
					SL 	= min( unL - aL, unR - aR ) ;
					SR 	= max( unL + aL, unR + aR ) ;
					SLm = min( SL, 0.0 ) ;
					SRp = max( SR, 0.0 ) ;


					/*---  Left Physical Flux  ---*/
					F0L = RhoL*unL ;
					F1L = RhoL*unL * uL + pL*m->PFM_CELL[ i ][ k ].nf[ 0 ]/IonMass ;
					F2L = RhoL*unL * vL + pL*m->PFM_CELL[ i ][ k ].nf[ 1 ]/IonMass ;
					F3L = RhoL*unL * wL + pL*m->PFM_CELL[ i ][ k ].nf[ 2 ]/IonMass ;
					F4L =   HL*unL ;


					/*---  Right Physical Flux  ---*/
					F0R = RhoR*unR ;
					F1R = RhoR*unR * uR + pR*m->PFM_CELL[ i ][ k ].nf[ 0 ]/IonMass ;
					F2R = RhoR*unR * vR + pR*m->PFM_CELL[ i ][ k ].nf[ 1 ]/IonMass ;
					F3R = RhoR*unR * wR + pR*m->PFM_CELL[ i ][ k ].nf[ 2 ]/IonMass ;
					F4R =   HR*unR ;

					/*---  Numerical Flux, Note: outward flux  ---*/
				    F0 = ( SRp*F0L - SLm*F0R + (SLm*SRp)*(var->U0[iSpecies][ j ]-var->U0[iSpecies][ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F1 = ( SRp*F1L - SLm*F1R + (SLm*SRp)*(var->U1[iSpecies][ j ]-var->U1[iSpecies][ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F2 = ( SRp*F2L - SLm*F2R + (SLm*SRp)*(var->U2[iSpecies][ j ]-var->U2[iSpecies][ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F3 = ( SRp*F3L - SLm*F3R + (SLm*SRp)*(var->U3[iSpecies][ j ]-var->U3[iSpecies][ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F4 = ( SRp*F4L - SLm*F4R + (SLm*SRp)*(var->U4[iSpecies][ j ]-var->U4[iSpecies][ i ]) )/(SRp-SLm);//+1.e-19) ;

					GradUvel = fabs(uL - uR )/m->PFM_CELL[ i ][ k ].dDist ;
					GradU4   = fabs(HL - HR )/m->PFM_CELL[ i ][ k ].dDist ;

				    /*---  Adding to residue  ---*/
				    Res[ 0 ][ i ] += LogicalSwitch*F0*m->PFM_CELL[ i ][ k ].dArea ;
				    Res[ 1 ][ i ] += LogicalSwitch*F1*m->PFM_CELL[ i ][ k ].dArea ;
				    Res[ 2 ][ i ] += LogicalSwitch*F2*m->PFM_CELL[ i ][ k ].dArea ;
				    Res[ 3 ][ i ] += LogicalSwitch*F3*m->PFM_CELL[ i ][ k ].dArea ;
				    Res[ 4 ][ i ] += LogicalSwitch*F4*m->PFM_CELL[ i ][ k ].dArea ;

	 			} else {/*--- For discontuity face ---*/

	 				switch ( config->Species[ iSpecies ].Type ){

						case 0:/*--- Electron ---*/

						break;

						case 1:/*--- Ion ---*/

							/*--- Left state ---*/
							RhoL = var->U0[iSpecies][ i ] ;

							uL 	 = var->U1[iSpecies][ i ]/RhoL ;
							vL 	 = var->U2[iSpecies][ i ]/RhoL ;
							wL 	 = var->U3[iSpecies][ i ]/RhoL ;	
							unL  = uL*m->PFM_CELL[ i ][ k ].nf[ 0 ] + vL*m->PFM_CELL[ i ][ k ].nf[ 1 ] + wL*m->PFM_CELL[ i ][ k ].nf[ 2 ] ;

							//pL 	 = (config->Species[ iSpecies ].Gamma-1.0)*( var->U4[iSpecies][ i ]-0.5*config->Species[ iSpecies ].Mass_Kg*RhoL*(uL*uL+vL*vL+wL*wL) ) ;
							pL 	= Pressure[ i ] ;
							HL 	= var->U4[iSpecies][ i ] ;//+	pL ;

							/*---  Physical Flux  ---*/
							unL = max ( 0.0, unL ) ;
							F0 = RhoL*unL ;
							F1 = RhoL*unL * uL ;//+ pL*m->PFM_CELL[ i ][ k ].nf[ 0 ]/config->Species[ iSpecies ].Mass_Kg ;
							F2 = RhoL*unL * vL ;//+ pL*m->PFM_CELL[ i ][ k ].nf[ 1 ]/config->Species[ iSpecies ].Mass_Kg ;
							F3 = RhoL*unL * wL ;//+ pL*m->PFM_CELL[ i ][ k ].nf[ 2 ]/config->Species[ iSpecies ].Mass_Kg ;
							F4 =   HL*unL ;

							/*---  Adding to residue  ---*/
				    		Res[ 0 ][ i ] += LogicalSwitch*F0*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 1 ][ i ] += LogicalSwitch*F1*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 2 ][ i ] += LogicalSwitch*F2*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 3 ][ i ] += LogicalSwitch*F3*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 4 ][ i ] += LogicalSwitch*F4*m->PFM_CELL[ i ][ k ].dArea ;


						break;

							
						case 2:	/*--- Neutral, Diffusion flux ---*/

						break;

						default:
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id) == "NEUMANN" ){
	 				//do nothing?????????????????????????????????????????
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						case 0:/*--- Electron ---*/

						break;

						case 1:/*--- Ion ---*/

							/*--- Left state ---*/
							RhoL = var->U0[iSpecies][ i ] ;
							uL 	 = var->U1[iSpecies][ i ]/RhoL ;
							vL 	 = var->U2[iSpecies][ i ]/RhoL ;
							wL 	 = var->U3[iSpecies][ i ]/RhoL ;

							unL  = uL*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
								 + vL*m->PFM_CELL[ i ][ k ].nf[ 1 ] 
								 + wL*m->PFM_CELL[ i ][ k ].nf[ 2 ] ;
							//unL = unL + GradUvel*m->PFM_CELL[ i ][ k ].dPPf ;

							//pL 	 = (config->Species[ iSpecies ].Gamma-1.0)*( var->U4[iSpecies][ i ]-0.5*config->Species[ iSpecies ].Mass_Kg*RhoL*(uL*uL+vL*vL+wL*wL) ) ;
							pL = Pressure[ i ] ;
							HL = var->U4[iSpecies][ i ] ;//+	pL ;
							//HL = HL + GradU4*m->PFM_CELL[ i ][ k ].dPPf ;

							/*---  Physical Flux  ---*/
							unL = max ( 0.0, unL ) ;
							F0 = RhoL*unL ;
							F1 = RhoL*unL * uL ;//+ pL*m->PFM_CELL[ i ][ k ].nf[ 0 ]/config->Species[ iSpecies ].Mass_Kg ;
							F2 = RhoL*unL * vL ;//+ pL*m->PFM_CELL[ i ][ k ].nf[ 1 ]/config->Species[ iSpecies ].Mass_Kg ;
							F3 = RhoL*unL * wL ;//+ pL*m->PFM_CELL[ i ][ k ].nf[ 2 ]/config->Species[ iSpecies ].Mass_Kg ;
							F4 =   HL*unL ;

							/*---  Adding to residue  ---*/
				    		Res[ 0 ][ i ] += LogicalSwitch*F0*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 1 ][ i ] += LogicalSwitch*F1*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 2 ][ i ] += LogicalSwitch*F2*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 3 ][ i ] += LogicalSwitch*F3*m->PFM_CELL[ i ][ k ].dArea ;
				    		Res[ 4 ][ i ] += LogicalSwitch*F4*m->PFM_CELL[ i ][ k ].dArea ;

						break;

						/*--- Neutral, Diffusion flux ---*/	
						case 2:	
						break;

						default:
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}
	 		var->Momentum_Term[0][ i ] = Res[ 1 ][ i ] ;

	 		var->Energy_Term[ 0 ][ i ] = Res[ 4 ][ i ] ;
	 		var->Energy_Term[ 4 ][ i ] = Pressure[ i ] ;
	 		var->Energy_Term[ 5 ][ i ] = var->U1[iSpecies][ i ]/var->U0[iSpecies][ i ] ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		Res[ 0 ][ i ] = 0.0 ;
			Res[ 1 ][ i ] = 0.0 ;
			Res[ 2 ][ i ] = 0.0 ;
			Res[ 3 ][ i ] = 0.0 ;
			Res[ 4 ][ i ] = 0.0 ;
	 	}
	 	//cout<<"R2: "<<Res[ 2 ][ i ]<<endl;
	}//Cell Loop
	//exit(1) ;
}
void CFluidModel::ContinuityIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{

	double LogicalSwitch = 0.0, SourceSink=0.0, Source=0.0 ;

	Cell *Cell_i ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*Cell_i->r[0] ; 

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

	 		/*--- Source/Sink term ---*/
	 		SourceSink = (double)*( var->ReactionRatePoint[iSpecies] + i  )/var->Ref_SS ;
	 		Source = (-1.0)*SourceSink*Cell_i->volume ;
	 		Res[ 0 ][ i ] += LogicalSwitch*Source ;
	 		Res[ 0 ][ i ]  = Res[ 0 ][ i ]/Cell_i->volume/LogicalSwitch ;

			/*--- Integrate ---*/
			var->U0[ iSpecies ][ i ]  = var->PreU0[iSpecies][ i ] - var->Dt * Res[ 0 ][ i ] ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		var->U0[iSpecies][ i ] = 0.0 ;
	 	}
	}//Cell Loop
	var->U0[ iSpecies ] = var->U0[ iSpecies ] ;
}
void CFluidModel::MomentumIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{

	double U=0.0, V=0.0, W=0.0, VTOT2=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0 ;
	double IonMass = config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass ;
	double LogicalSwitch=0.0 ;
	double PlasmaParmeter = var->Ref_N*var->Ref_L*var->Ref_L*var->Ref_L ;

	Cell *Cell_i ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*Cell_i->r[0] ; 

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			U0 = var->PreU0[iSpecies][ i ] ;
			U1 = var->U1[iSpecies][ i ] ;
			U2 = var->U2[iSpecies][ i ] ;
			U3 = var->U3[iSpecies][ i ] ;
			U4 = var->U4[iSpecies][ i ] ;
			//P = (config->Species[ iSpecies ].Gamma-1.0)*( U4 - 0.5*IonMass*( U1*U1+U2*U2+U3*U3 )/U0 ) ;

			/*--- Force term { (1/m)qnE } ---*/
        	Res[ 1 ][ i ] +=  LogicalSwitch * (-config->Species[ iSpecies ].Charge*var->Qe*U0*var->EField[ 0 ][ i ]/IonMass) * Cell_i->volume ;
        	Res[ 2 ][ i ] +=  LogicalSwitch * (-config->Species[ iSpecies ].Charge*var->Qe*U0*var->EField[ 1 ][ i ]/IonMass) * Cell_i->volume ;
        	Res[ 3 ][ i ] +=  LogicalSwitch * (-config->Species[ iSpecies ].Charge*var->Qe*U0*var->EField[ 2 ][ i ]/IonMass) * Cell_i->volume ;
	 		var->Momentum_Term[1][ i ] = LogicalSwitch * (-config->Species[ iSpecies ].Charge*var->Qe*U0*var->EField[ 0 ][ i ]/IonMass) * Cell_i->volume ;
			/*--- Collision term ---*/  
			U = U1/U0 ;
			V = U2/U0 ;
			W = U3/U0 ;
			VTOT2 = U*U + V*V + W*W ;

            Res[ 1 ][ i ]   += LogicalSwitch * (-Mx[ i ]*PlasmaParmeter/IonMass) * Cell_i->volume ;
			Res[ 2 ][ i ]   += LogicalSwitch * (-My[ i ]*PlasmaParmeter/IonMass) * Cell_i->volume ;
			Res[ 3 ][ i ]   += LogicalSwitch * (-Mz[ i ]*PlasmaParmeter/IonMass) * Cell_i->volume ;
			var->Momentum_Term[ 2 ][ i ] = LogicalSwitch * (-Mx[ i ]/IonMass) * Cell_i->volume ;
			/*--- Y-Axis-symmetric ---*/
	        Res[ 1 ][ i ] += (-1.0) * Omaga * (Pressure[ i ]/IonMass) * Cell_i->volume ;

	 		/*--- Integrate ---*/
	 		Res[ 1 ][ i ]  = Res[ 1 ][ i ]/Cell_i->volume/LogicalSwitch ;
	 		Res[ 2 ][ i ]  = Res[ 2 ][ i ]/Cell_i->volume/LogicalSwitch ;
	 		Res[ 3 ][ i ]  = Res[ 3 ][ i ]/Cell_i->volume/LogicalSwitch ;
			var->U1[ iSpecies ][ i ]  = var->PreU1[iSpecies][ i ] - var->Dt * Res[ 1 ][ i ] ;
			var->U2[ iSpecies ][ i ]  = var->PreU2[iSpecies][ i ] - var->Dt * Res[ 2 ][ i ] ;
			var->U3[ iSpecies ][ i ]  = var->PreU3[iSpecies][ i ] - var->Dt * Res[ 3 ][ i ] ;


	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		var->U1[iSpecies][ i ] = 0.0 ;
	 		var->U2[iSpecies][ i ] = 0.0 ;
	 		var->U3[iSpecies][ i ] = 0.0 ;
	 	}
	}//Cell Loop

	var->U1[ iSpecies ] = var->U1[ iSpecies ] ;
	var->U2[ iSpecies ] = var->U2[ iSpecies ] ;
	var->U3[ iSpecies ] = var->U3[ iSpecies ] ;
}
void CFluidModel::EnergyDensityIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{

	double U=0.0, V=0.0, W=0.0, VTOT2=0.0, JHeating=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0, F4=0.0, CapitalE=0.0 ;
	double RhoR=0.0, RhoL=0.0, uL=0.0, uR=0.0, vL=0.0, vR=0.0, wL=0.0, wR=0.0, unL=0.0, unR=0.0, pL=0.0, pR=0.0, aL=0.0, aR=0.0, TiL=0.0, TiR=0.0, TeL=0.0, TeR=0.0, eL=0.0, eR=0.0 ;
	double qz=0.0, qx=0.0, qy=0.0 ;

	double LogicalSwitch=0.0, Ad_dPN=0.0 ;

	double CrossSection = 50E-20 ; //Grave
	Cell *Cell_i, *Cell_j ;
	int j=0 ;
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*Cell_i->r[0] ; 

		Cell_i = plasma.get_cell(i) ;

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Central state ---*/
			U0 = var->U0[iSpecies][ i ] ;
			U1 = var->U1[iSpecies][ i ] ;
			U2 = var->U2[iSpecies][ i ] ;
			U3 = var->U3[iSpecies][ i ] ;
			U4 = var->U4[iSpecies][ i ] ;
			U = U1/U0 ;
			V = U2/U0 ;
			W = U3/U0 ;
			//cout<<"U3: "<<U3<<endl;
			VTOT2 = U*U + V*V + W*W ;
			/*--- Divergence U ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = plasma.get_cell(j) ;

				if ( plasma.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					qx = 0.5 * ( var->Kappa[ i ]*var->Qe*var->GradT[ iSpecies ][ 0 ][ i ]
							 	+var->Kappa[ j ]*var->Qe*var->GradT[ iSpecies ][ 0 ][ j ] ) ;

					qy = 0.5 * ( var->Kappa[ i ]*var->Qe*var->GradT[ iSpecies ][ 1 ][ i ]
							 	+var->Kappa[ j ]*var->Qe*var->GradT[ iSpecies ][ 1 ][ j ] ) ;
					qz = 0.0 ;

					//F4 = qx*m->PFM_CELL[ i ][ k ].nf[ 0 ] + qy*m->PFM_CELL[ i ][ k ].nf[ 1 ] + qz*m->PFM_CELL[ i ][ k ].nf[ 2 ] ;
					Ad_dPN = m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
					F4 = -Ad_dPN*var->T[ iSpecies ][ i ]*var->Qe + Ad_dPN*var->T[ iSpecies ][ j ]*var->Qe ;
					//qn = 0.0 ;
					Res[ 4 ][ i ] += LogicalSwitch*F4*m->PFM_CELL[ i ][ k ].dArea ;
					var->Energy_Term[ 1 ][ i ] =  LogicalSwitch*F4*m->PFM_CELL[ i ][ k ].dArea ;

	 			}
	 		}

			/*--- Joule Heating ---*/  
			//UNSTEADY + CONVECTION - Joule Heating - collision = 0 
			JHeating = var->Qe*config->Species[iSpecies].Charge*( var->EField[ 0 ][ i ]*U1 + var->EField[ 1 ][ i ]*U2 ) ;
			Res[ 4 ][ i ] +=  LogicalSwitch * (-JHeating) * Cell_i->volume ;
			var->Energy_Term[ 2 ][ i ] =  LogicalSwitch*(-JHeating)*Cell_i->volume ;


			//Kunshner's
			Res[ 4 ][ i ] += LogicalSwitch*CollisionIntegral[ i ]*Cell_i->volume ;
			var->Energy_Term[ 3 ][ i ] =  LogicalSwitch*CollisionIntegral[ i ]*Cell_i->volume ;


	 		/*--- Integrate ---*/
	 		Res[ 4 ][ i ]  = Res[ 4 ][ i ]/Cell_i->volume/LogicalSwitch ;
			var->U4[ iSpecies ][ i ]  = var->PreU4[iSpecies][ i ] - var->Dt * Res[ 4 ][ i ] ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		var->U4[ iSpecies ][ i ] = 0.0 ;
	 	}
	}//Cell Loop
	//exit(1) ;
	var->U4[ iSpecies ] = var->U4[ iSpecies ] ;
}
void CFluidModel::CalculateTemperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0 ;
	double IonMass = config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass ;

	Cell *Cell_i ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;
	
		if ( plasma.get_cell_typename( Cell_i->data_id ) != "PLASMA" ){

			var->T[iSpecies][ i ] = 0.0 ;
			Pressure[ i ] = 0.0 ;

		} else {

			U0 = var->U0[iSpecies][ i ] ;
			U1 = var->U1[iSpecies][ i ] ;
			U2 = var->U2[iSpecies][ i ] ;
			U3 = var->U3[iSpecies][ i ] ;
			U4 = var->U4[iSpecies][ i ] ;

			Pressure[ i ] = (config->Species[ iSpecies ].Gamma-1.0)*( U4 - 0.5*IonMass*( U1*U1+U2*U2+U3*U3 )/U0 ) ;
			var->T[iSpecies][ i ] =  Pressure[ i ]/U0/var->Qe ;
			if ( var->T[iSpecies][ i ] < 0.0  ) var->T[iSpecies][ i ] = 0.0256 ;
		}
	}//Loop over all cells
	Pressure = Pressure ;
	var->T[iSpecies] = var->T[iSpecies] ;
}
void CFluidModel::CalculateTotalEnergy( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	double P=0.0, U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0 ;
	double IonMass = config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass ;

	Cell *Cell_i ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		
		Cell_i  = plasma.get_cell( i ) ;

		if ( plasma.get_cell_typename( Cell_i->data_id ) != "PLASMA" ){

			var->T[iSpecies][ i ] = 0.0 ;

		} else {

			U0 = var->U0[iSpecies][ i ] ;
			U1 = var->U1[iSpecies][ i ] ;
			U2 = var->U2[iSpecies][ i ] ;
			U3 = var->U3[iSpecies][ i ] ;
            P = U0*var->Qe*config->Species[iSpecies].InitialTemperature ;
			var->U4[iSpecies][ i ] = P/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*IonMass*( U1*U1+U2*U2+U3*U3)/U0 ;
		}
	}//Loop over all cells
	var->U4[iSpecies] = var->U4[iSpecies] ;
	//exit(1);
}
void CFluidModel::Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	int j=0, NeighborCellIndex=0 ;
	double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0 ;
	Cell *Cell_i, *Cell_j ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		Gx = 0.0 ; Gy = 0.0 ;

		if ( plasma.get_cell_typename( Cell_i->data_id ) != "PLASMA" ){

			var->GradU0[ iSpecies ][ 0 ][ i ] = 0.0 ;
			var->GradU0[ iSpecies ][ 1 ][ i ] = 0.0 ;	

		} else {

			/*--- Loop over neighbor "faces" ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = plasma.get_cell(j) ;
				
				if ( plasma.get_cell_typename( Cell_j->data_id ) != "PLASMA" ) {//For discontinued face, apply neumann

					GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
					GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;

					BC_Value = var->U0[iSpecies][ i ] + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

					dVar = BC_Value - var->U0[iSpecies][ i ] ;

				}else{

					dVar = var->U0[iSpecies][ j ] - var->U0[iSpecies][ i ] ;

				}

				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	     		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

			}//End iCell

			/*--- Loop over boundary faces ---*/
			for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
				GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;

				BC_Value = var->U0[iSpecies][ i ] + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

				dVar = BC_Value - var->U0[ iSpecies ][ i ] ;

				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	      		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

			}
			var->GradU0[ iSpecies ][ 0 ][ i ] = Gx ;
			var->GradU0[ iSpecies ][ 1 ][ i ] = Gy ;
		}
	}//Loop over all cells

	/*--- Update ghost cells ---*/
	var->GradU0[ iSpecies ][ 0 ] = var->GradU0[ iSpecies ][ 0 ] ;
	var->GradU0[ iSpecies ][ 1 ] = var->GradU0[ iSpecies ][ 1 ] ;

}
void CFluidModel::Calculate_Gradient_T( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	int j=0, NeighborCellIndex=0 ;
	double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0 ;

	Cell *Cell_i, *Cell_j ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		Gx = 0.0 ; Gy = 0.0 ;


		if ( plasma.get_cell_typename( Cell_i->data_id ) != "PLASMA" ){

			var->GradT[ iSpecies ][ 0 ][ i ] = 0.0 ;
			var->GradT[ iSpecies ][ 1 ][ i ] = 0.0 ;	

		} else {

			/*--- Loop over neighbor "faces" ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = plasma.get_cell(j) ;
				
				if ( plasma.get_cell_typename( Cell_j->data_id ) != "PLASMA" ) {//For discontinued face, apply neumann

					GVarP[ 0 ] = var->GradT[ iSpecies ][ 0 ][ i ] ;
					GVarP[ 1 ] = var->GradT[ iSpecies ][ 1 ][ i ] ;

					BC_Value = var->T[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;
					dVar = BC_Value - var->T[iSpecies][ i ] ;

				}else{

					dVar = var->T[iSpecies][ j ] - var->T[iSpecies][ i ] ;

				}

				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	     		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

			}//End iCell

			/*--- Loop over boundary faces ---*/
			for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				GVarP[ 0 ] = var->GradT[ iSpecies ][ 0 ][ i ] ;
				GVarP[ 1 ] = var->GradT[ iSpecies ][ 1 ][ i ] ;

				BC_Value = var->T[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

				dVar = BC_Value - var->T[ iSpecies ][ i ] ;

				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	      		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

			}
			var->GradT[ iSpecies ][ 0 ][ i ] = Gx ;
			var->GradT[ iSpecies ][ 1 ][ i ] = Gy ;
		}
	}//Loop over all cells

	/*--- Update ghost cells ---*/
	var->GradT[ iSpecies ][ 0 ] = var->GradT[ iSpecies ][ 0 ] ;
	var->GradT[ iSpecies ][ 1 ] = var->GradT[ iSpecies ][ 1 ] ;

}
void CFluidModel::Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		var->GradU0[ iSpecies ][ 0 ][ i ] = 0.0 ;
		var->GradU0[ iSpecies ][ 1 ][ i ] = 0.0 ;
	}//Loop over all cells

	/*--- Update ghost cells ---*/
	var->GradU0[ iSpecies ][ 0 ] = var->GradU0[ iSpecies ][ 0 ] ;
	var->GradU0[ iSpecies ][ 1 ] = var->GradU0[ iSpecies ][ 1 ] ;
}
void CFluidModel::CalculateSurfaceCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = plasma.get_cell(j) ;

				if ( plasma.get_cell_typename( Cell_j->data_id ) == "DIELECTRIC" ) {

					m->PFM_CELL[ i ][ k ].SurfaceCharge += var->Dt*var->Qe*config->Species[iSpecies].Charge
					*fabs( var->U1[ iSpecies ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
						+  var->U2[ iSpecies ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

				}//discontiuity face

	 		}//End bulk 

	 	/*--- Loop over DIELECTRIC cells ---*/
	 	} else if( plasma.get_cell_typename( Cell_i->data_id ) == "DIELECTRIC" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = plasma.get_cell(j) ;

				if ( plasma.get_cell_typename( Cell_j->data_id ) == "PLASMA" ) {

					m->PFM_CELL[ i ][ k ].SurfaceCharge += var->Dt*var->Qe*config->Species[iSpecies].Charge
					*fabs( var->U1[ iSpecies ][ j ]*m->PFM_CELL[ i ][ k ].nf[ 0 ]*(-1.0) 
						+  var->U2[ iSpecies ][ j ]*m->PFM_CELL[ i ][ k ].nf[ 1 ]*(-1.0) ) ;

				}//discontiuity face

	 		}//End bulk 
	 	}
	}//Cell Loop
}
void CFluidModel::CalculateCondCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	Cell *Cell_i ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			var->CondJD[iSpecies][0][i] =  config->Species[ iSpecies ].Charge*var->Qe*var->U1[ iSpecies ][ i ] ;
			var->CondJD[iSpecies][1][i] =  config->Species[ iSpecies ].Charge*var->Qe*var->U2[ iSpecies ][ i ] ;
			var->CondJD[iSpecies][2][i] =  config->Species[ iSpecies ].Charge*var->Qe*var->U3[ iSpecies ][ i ] ;

			var->TotalJD[0][i] += var->CondJD[iSpecies][0][ i ] ;
			var->TotalJD[1][i] += var->CondJD[iSpecies][1][ i ] ;
			var->TotalJD[2][i] += var->CondJD[iSpecies][2][ i ] ;
					
		 	/*--- Loop over SOLID cells ---*/
	 	} else {

	 		var->CondJD[iSpecies][0][i] = 0.0 ;
	 		var->CondJD[iSpecies][1][i] = 0.0 ;
	 		var->CondJD[iSpecies][2][i] = 0.0 ;

	 	}
	}//Cell Loop
	/*--- Update ghost cells ---*/
	//may not need to update ghost cells.
}
void CFluidModel::CalculateThermal2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	double IonMass = config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass ;

	Cell *Cell_i ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){
        	Thermal2[ i ] = 8.0*var->Qe*var->T[iSpecies][ i ]/var->PI/IonMass ;
	 	 
		} else{
			Thermal2[ i ] = 0.0 ;
		}
	}//Cell Loop
		Thermal2 = Thermal2 ;
}
void CFluidModel::CalculateIonNeutralCollisionFrequency( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{

	double U=0.0, V=0.0, W=0.0, VTOT2=0.0, ReducedMass=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0, A=0.0, B=0.0, C=0.0, SourceSink=0.0 ;

	Cell *Cell_i ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Mx[ i ] = 0.0 ;
		My[ i ] = 0.0 ;
		Mz[ i ] = 0.0 ;
		CollisionIntegral[ i ] = 0.0 ;

		Cell_i  = plasma.get_cell( i ) ;


		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			U0 = var->U0[iSpecies][ i ] ;
			U1 = var->U1[iSpecies][ i ] ;
			U2 = var->U2[iSpecies][ i ] ;
			U3 = var->U3[iSpecies][ i ] ;
			U4 = var->U4[iSpecies][ i ] ;

			U = U1/U0 ;
			V = U2/U0 ;
			W = U3/U0 ;
			VTOT2 = U*U + V*V + W*W ;

	        for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

	            if ( 
	            	//config->Species[ jSpecies ].Type == NEUTRAL or 
	            	config->Species[ jSpecies ].Type == BACKGROUND ){

					if ( Grave ){

						CollisionFreq[ jSpecies ][ i ] = var->U0[jSpecies][ i ]*(CrossSection)*sqrt( VTOT2 ) ;
						Mx[ i ]   += -0.5*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)*var->PI*U1*CollisionFreq[ jSpecies ][ i ] ;
						My[ i ]   += -0.5*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)*var->PI*U2*CollisionFreq[ jSpecies ][ i ] ;
						Mz[ i ]   += -0.5*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)*var->PI*U3*CollisionFreq[ jSpecies ][ i ] ;

					}else{

						ReducedMass = ( (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) * (config->Species[ jSpecies ].Mass_Kg/var->Ref_Mass) )
									/ ( (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) + (config->Species[ jSpecies ].Mass_Kg/var->Ref_Mass) ) ;

						CollisionFreq[ jSpecies ][ i ] = var->U0[jSpecies][ i ]*CrossSection*sqrt( 1.777777778*Thermal2[ i ] + VTOT2 ) ;
						Mx[ i ]   += (-1.0)*ReducedMass*U1*CollisionFreq[ jSpecies ][ i ] ;
						My[ i ]   += (-1.0)*ReducedMass*U2*CollisionFreq[ jSpecies ][ i ] ;
						Mz[ i ]   += (-1.0)*ReducedMass*U3*CollisionFreq[ jSpecies ][ i ] ;
					}	
					//cout<<"iSpecies"
					/*--- Collision Integral ---*/
					A = 2.0*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)/(config->Species[ jSpecies ].Mass_Kg/var->Ref_Mass)*CollisionFreq[ jSpecies ][ i ] ;
					B = 0.5*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)*U0*VTOT2 ;
					B = 0.0 ;
					C = 1.5*U0*var->Qe*(var->T[iSpecies][i] - var->T[jSpecies][i]) ;
					CollisionIntegral[ i ] += A*(B+C) ;
					//CollisionIntegral[ i ] += A*U4 ;
					/*
					ReducedMass = 2.0*( config->Species[ iSpecies ].Mass_Kg * config->Species[ jSpecies ].Mass_Kg )
									/ ( config->Species[ iSpecies ].Mass_Kg + config->Species[ jSpecies ].Mass_Kg )
									/ ( config->Species[ iSpecies ].Mass_Kg + config->Species[ jSpecies ].Mass_Kg ) ;
					A = U0 * var->U0[jSpecies][ i ] * CrossSection * ReducedMass ;
					B = 0.5*config->Species[ iSpecies ].Mass_Kg*pow(1.351283845*Thermal2[i]+VTOT2, 1.5) ;
					C = sqrt(Thermal2[i]+VTOT2)*(3.0/2.0)*var->Qe*(var->T[jSpecies][i] - var->T[iSpecies][i]) ;
					CollisionIntegral[ i ] += A*(B-C) ;
					*/
	            } 
	        }//End j-Species

	 	} 
	}//Cell Loop
	//exit(1) ;
	for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
		CollisionFreq[ jSpecies ] = CollisionFreq[ jSpecies ] ;
	}
	//exit(1) ;
}
void CFluidModel::CalculateCollisionIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{

	double U=0.0, V=0.0, W=0.0, VTOT2=0.0, ReducedMass=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0, A=0.0, B=0.0, C=0.0, SourceSink=0.0 ;

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		CollisionIntegral[ i ] = 0.0 ;

		Cell_i  = plasma.get_cell( i ) ;

		/*--- Loop over PLASMA cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			U0 = var->U0[iSpecies][ i ] ;
			U1 = var->U1[iSpecies][ i ] ;
			U2 = var->U2[iSpecies][ i ] ;
			U3 = var->U3[iSpecies][ i ] ;
			U4 = var->U4[iSpecies][ i ] ;

			U = U1/U0 ;
			V = U2/U0 ;
			W = U3/U0 ;
			VTOT2 = U*U + V*V + W*W ;

	        for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

	            if ( 
	            	//config->Species[ jSpecies ].Type == NEUTRAL or 
	            	config->Species[ jSpecies ].Type == BACKGROUND ){

					if ( Grave ){

						CollisionFreq[ jSpecies ][ i ] = var->U0[jSpecies][ i ]*(CrossSection)*sqrt( VTOT2 ) ;


					}else{

						ReducedMass = ( (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) * (config->Species[ jSpecies ].Mass_Kg/var->Ref_Mass) )
									/ ( (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) + (config->Species[ jSpecies ].Mass_Kg/var->Ref_Mass) ) ;

						CollisionFreq[ jSpecies ][ i ] = var->U0[jSpecies][ i ]*CrossSection*sqrt( 1.777777778*Thermal2[ i ] + VTOT2 ) ;

					}	
					//cout<<"iSpecies"
					/*--- Collision Integral ---*/
					A = 2.0*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)/(config->Species[ jSpecies ].Mass_Kg/var->Ref_Mass)*CollisionFreq[ jSpecies ][ i ] ;
					B = 0.5*(config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass)*U0*VTOT2 ;
					B = 0.0 ;
					C = 1.5*U0*var->Qe*(var->T[iSpecies][i] - var->T[jSpecies][i]) ;
					CollisionIntegral[ i ] += A*(B+C) ;
					//CollisionIntegral[ i ] += A*U4 ;
					/*
					ReducedMass = 2.0*( config->Species[ iSpecies ].Mass_Kg * config->Species[ jSpecies ].Mass_Kg )
									/ ( config->Species[ iSpecies ].Mass_Kg + config->Species[ jSpecies ].Mass_Kg )
									/ ( config->Species[ iSpecies ].Mass_Kg + config->Species[ jSpecies ].Mass_Kg ) ;
					A = U0 * var->U0[jSpecies][ i ] * CrossSection * ReducedMass ;
					B = 0.5*config->Species[ iSpecies ].Mass_Kg*pow(1.351283845*Thermal2[i]+VTOT2, 1.5) ;
					C = sqrt(Thermal2[i]+VTOT2)*(3.0/2.0)*var->Qe*(var->T[jSpecies][i] - var->T[iSpecies][i]) ;
					CollisionIntegral[ i ] += A*(B-C) ;
					*/
	            } 
	        }//End j-Species

	 	} 
	}//Cell Loop


}
void CFluidModel::CalculateKappa( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{

	double TotalCollisionFreq=0.0 ;
	double IonMass = config->Species[iSpecies].Mass_Kg/var->Ref_Mass ;
	Cell *Cell_i ;
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		/*--- Loop over PLASMA cells ---*/
		Cell_i  = plasma.get_cell( i ) ;

		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

			TotalCollisionFreq = 0.0 ;
			
	        for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

	            if ( config->Species[ jSpecies ].Type == NEUTRAL ){
	            	TotalCollisionFreq += CollisionFreq[ jSpecies ][ i ] ;
	            }else if( config->Species[ jSpecies ].Type == BACKGROUND ) {
	            	TotalCollisionFreq += CollisionFreq[ jSpecies ][ i ] ;
	            }

	        }//End Species
	        var->Kappa[ i ] = -2.5*Pressure[ i ]/IonMass/TotalCollisionFreq ;
	        //cout<<var->Kappa[ i ]<<endl;
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
			var->Kappa[ i ] = 0.0 ;
	 	}
	 	//var->DUDT[ 0 ][ i ] = TotalCollisionFreq ;
	}//Cell Loop
	var->Kappa = var->Kappa ;
}


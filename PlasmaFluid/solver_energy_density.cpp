#include "solver_energy_density.hpp"

using namespace std ;
CEnergyDensity::CEnergyDensity()
{
}
void CEnergyDensity::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, int index )
{
	iSpecies = index ;
	C53 = 5.0/3.0 ;
	C43 = 4.0/3.0 ;
	if ( mpi_id == 0 ){
		cout<<"Creat "<<config->Species[iSpecies].Name<<" continuity solver, index: "<<index<<", charge: "<<config->Species[ index ].Charge<<endl;
	} 
	
	/*--- PETSc Solver ---*/
		s.set_local_unknown_number ( m->local_cell_number ) ;
		d_nnz =	new int  [ m->local_cell_number ] ;
		o_nnz =	new int  [ m->local_cell_number ] ;
		for ( int i = 0 ; i < m->local_cell_number ; i++ ) {
		 	d_nnz[ i ] = 5 ;
		 	o_nnz[ i ] = 5 ;
		}
		s.set_nz( 5, d_nnz, o_nnz ) ;
		s.init() ;
		//s.reuse_preconditioner( true ) ;
		NormalizeCoeff = new double [ m->local_cell_number ] ;
	/*--- PETSc Solver ---*/	
	Correction 		 = config->Equation[ config->Species[ index ].Type ].Correction ;
	WallType = config->Equation[ config->Species[ index ].Type ].WallBoundaryType ;
	
	fixTe = false ;
	//TG = true ;
	TG = false ;
	Insert = false ;
	eLOSS = 1 ;//1: from chemistry module, 2: from table.
	Reflec = config->ReflectionCoefficient ;
	its=0 ;
}
void CEnergyDensity::Solver( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{
	Zero_Gradient( m, variable ) ;

	/*--- 1st solve the equation w/o cross diffusion term ---*/
//cout<<"CEnergyDensity: A"<<endl;
		switch( WallType ){
			case 0: //default
				Bulid_A_B_1st_default( m, config, variable ) ;
				break ;

			case 1://neumann
				exit(1) ; //not support yet
				break;

			case 2://zero number density
				Bulid_A_B_1st_zero( m, config, variable ) ;
				break;
			default:
				break;
		}

		s.solve() ; 
		its = s.its ;

		for ( int i = 0 ; i < m->local_cell_number ; i++ ) {
			variable->U4[iSpecies][i] = s.x[i] ;
		}
		variable->U4[iSpecies] = variable->U4[iSpecies] ;
//cout<<"CEnergyDensity: B"<<endl;

	/*--- 2nd solve the equation w/ cross diffusion term ---*/
//cout<<"CEnergyDensity: C"<<endl;
		for ( int k = 0 ; k < Correction ; k++ ) {

//cout<<"CEnergyDensity: C1"<<endl;
			Calculate_Gradient( m, variable ) ;
//cout<<"CEnergyDensity: C2"<<endl;
			Bulid_A_B_2nd( m, config, variable ) ;
//cout<<"CEnergyDensity: C3"<<endl;
		  	s.solve(); 
		  	its += s.its ;
//cout<<"CEnergyDensity: C4"<<endl;
			for ( int i = 0 ; i < m->local_cell_number ; i++ ) {
				variable->U4[iSpecies][i] = s.x[i] ;
			}
//cout<<"CEnergyDensity: C5"<<endl;
			variable->U4[iSpecies] = variable->U4[iSpecies] ;
		}
//cout<<"CEnergyDensity: D"<<endl;
		//Calculate_Gradient( m, variable ) ;
//cout<<"CEnergyDensity: E"<<endl;
	/*--- calculate electron temperature ---*/
		CalculateTemperature( m, variable ) ;
//cout<<"CEnergyDensity: F"<<endl;
}
void CEnergyDensity::Bulid_A_B_1st_default( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0 ;
	double Diff = 0.0, Mobi=0.0, SourceSink=0.0, JdotE=0.0,TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	for( int i = 0 ; i < nCell ; i++ ) {

		iFace 	 = m->cell[ i ].face_number ;
		iCell 	 = m->cell[ i ].cell_number ;

		/*--- Reset  ---*/
		row 	 = m->cell[ i ].id ;
		col[ 0 ] = row ;
		ncol 	 = 1 ;

		Source 	 = 0.0 ;

		for( int k = 0 ; k < 5 ; k++ ) C[ k ] = 0.0 ;
		/*--- Loop over PLASMA cells ---*/
		if ( m->cell[ i ].type == PLASMA ){

			/*--- Time ---*/
			C[ 0 ] = m->cell[ i ].volume ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = m->Cell[ i ][ k ].NeighborCellId ;

				if ( m->cell[ j ].type == PLASMA ){

					/*--- S-G ---*/
					dL = m->Cell[ i ][ k ].dNPf / m->Cell[ i ][ k ].dDist ;
					dR = m->Cell[ i ][ k ].dPPf / m->Cell[ i ][ k ].dDist ;
					/*--- S-G ---*/
					// U = 0.5*( config->Species[ iSpecies ].Charge * var->EField[ 0 ][ i ] * var->Mobi[iSpecies][ i ] 
					// 		+ config->Species[ iSpecies ].Charge * var->EField[ 0 ][ j ] * var->Mobi[iSpecies][ j ] ) ;

					// V = 0.5*( config->Species[ iSpecies ].Charge  * var->EField[ 1 ][ i ] * var->Mobi[iSpecies][ i ] 
					// 		+ config->Species[ iSpecies ].Charge  * var->EField[ 1 ][ j ] * var->Mobi[iSpecies][ j ] ) ;
					U = dL*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ j ] * var->Mobi[iSpecies][ j ]  ;

					vn = U*m->Cell[ i ][ k ].nf[ 0 ] + V*m->Cell[ i ][ k ].nf[ 1 ] ;

					//Diff = 0.5*( var->Diff[iSpecies][ i ] + var->Diff[iSpecies][ j ] ) ; 
					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ; 

					Pe = vn*m->Cell[ i ][ k ].dDist/Diff ;

					col[ ncol ] = m->Cell[ i ][ k ].NeighborGlobalCellId ;

					if ( Pe < -ZERO  ){

						C[ 0 ] += C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;
						C[ncol] = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;

					} else if ( Pe > ZERO  ){

						C[ 0 ] += C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;
						C[ncol] = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;

					} else {

						Diff = - ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;  
						C[ 0 ] 	  += -C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist*var->Dt ;
						C[ncol]    =  C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist*var->Dt ;

					}

					// if ( fabs(Pe) < ZERO )
					// {
					// 	//Diff = -0.5*( var->Diff[iSpecies][ i ] + var->Diff[iSpecies][ j ] );
					// 	Diff = -( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );
					// 	C[ 0 ] 	  += -C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist*var->Dt ;
					// 	C[ncol]    =  C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist*var->Dt ;

					// } else {

					// 	f1 = Pe /( exp( Pe )-1.0 ) ;
					// 	f2 = f1 + Pe ;
					// 	C[ 0 ] 	  += -C53*Diff/m->Cell[ i ][ k ].dDist*(-f2)*m->Cell[ i ][ k ].dArea*var->Dt ;
					// 	C[ncol]    = -C53*Diff/m->Cell[ i ][ k ].dDist*( f1)*m->Cell[ i ][ k ].dArea*var->Dt ;
					// }
					// if ( TG ){
					// 	Mobi = 0.5*( var->Mobi[iSpecies][ i ] + var->Mobi[iSpecies][ j ] ) ;  
					// 	TempGradient = var->Qe*( var->T[ 0 ][ j ] - var->T[ 0 ][ i ] ) /m->Cell[ i ][ k ].dDist ;
     				//  C[ 0  ] += C53*0.5*( -Mobi*TempGradient/var->Qe )*m->Cell[ i ][ k ].dArea*var->Dt ;
     				//  C[ncol] += C53*0.5*( -Mobi*TempGradient/var->Qe )*m->Cell[ i ][ k ].dArea*var->Dt ;
					// }

		 			ncol ++ ;

	 			} else {//Discontinue face

	 				switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 0 ][ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 1 ][ i ] ;
	 						vn = C43*max( 0.0, U*m->Cell[ i ][ k ].nf[ 0 ]+V*m->Cell[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	 						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->EField[ 0 ][ i ]*0.5*m->Cell[ i ][ k ].dDist)/var->T[0][i]) ;
	 						
	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						Te_sec = config->SecondaryElectronEmissionEnergy ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 0 ][ i ]*m->Cell[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 1 ][ i ]*m->Cell[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 						}
	 						if ( m->cell[ j ].type == EMITTER ){
	 							SecondaryElectronEmission += var->EmitterFlux( i ) ;
	 						} 
	 						Source += 2.0*Te_sec*SecondaryElectronEmission*m->Cell[ i ][ k ].dArea*var->Dt ;

	 						C[ 0 ] +=  vn*m->Cell[ i ][ k ].dArea*var->Dt ;
	 						
						break;

						/*--- Ion, drift flux ---*/	
						case 1:
							cout<<"Not Support Ion energy yet"<<endl;exit(1) ;
						break;

						/*--- Neutral, Diffusion flux ---*/	
						case 2:	
							cout<<"Not Support Neutral energy yet"<<endl;exit(1) ;
						break;

						default:
							if( mpi_id == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}//For discontuity face
	 		}//End bulk face

			/*--- Loop over boundary faces ---*/
	 		for( int k = iCell ; k < iFace ; k++ ) {

	 			if( m->cell[ i ].face[ k ]->type == NEUMANN ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 0 ][ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 1 ][ i ] ;
	 						vn = C43*max( 0.0, U*m->Cell[ i ][ k ].nf[ 0 ]+V*m->Cell[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	 						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->EField[ 0 ][ i ]*0.5*m->Cell[ i ][ k ].dDist)/var->T[0][i]) ;
	 						
	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						Te_sec = config->SecondaryElectronEmissionEnergy ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 0 ][ i ]*m->Cell[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 1 ][ i ]*m->Cell[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 						}
	 						if ( m->cell[ i ].face[ k ]->type == EMITTER ) {
	 							SecondaryElectronEmission += var->EmitterFlux( i ) ;	
	 						}
	 					
	 						Source += 2.0*Te_sec*SecondaryElectronEmission*m->Cell[ i ][ k ].dArea*var->Dt ;
	 						C[ 0 ] +=  vn*m->Cell[ i ][ k ].dArea*var->Dt ;
	 						
						break;

						/*--- Ion, drift flux ---*/	
						case 1:
							cout<<"Not Support Ion energy yet"<<endl;exit(1) ;
						break;

						/*--- Neutral, Diffusion flux ---*/	
						case 2:	
							cout<<"Not Support Neutral energy yet"<<endl;exit(1) ;
						break;

						default:
							if( mpi_id == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		Source += var->PreU4[iSpecies][ i ] * m->cell[ i ].volume ;

			/*--- Joule heating, Note: since d_Te is in eV, there is no need to multiply the elementary charge ---*/
			JdotE = config->Species[iSpecies].Charge*( var->EField[ 0 ][ i ]*var->U1[ iSpecies ][ i ] 
	 												 + var->EField[ 1 ][ i ]*var->U2[ iSpecies ][ i ] ) ;

	 		/*--- energy loss term ---*/
	 		// switch ( eLOSS ) {
	 		// 	case 1:
	 				SourceSink = *(var->EnergySourcePoint + (i) ) ;
				// 	break;
				// case 2:
	 		// 		SourceSink = var->eEnergyLossTable.GetValueLog( var->T[iSpecies][ i ] ) * var->TotalNumberDensity[ i ]*var->U0[ iSpecies ][ i ] ;
				// 	break;
				// default:
	 		// 		SourceSink = *(var->EnergySourcePoint + (i) ) ;
		  //   		break;
	 		// }
	 		//cout<<"eLoss: "<<SourceSink<<endl;
	 		//exit(1) ;
	 		//Source += (-SourceSink)*m->cell[ i ].volume*var->Dt ;
	 		Source += (JdotE-SourceSink/var->Ref_ES)*m->cell[ i ].volume*var->Dt ;

	 		var->eEnergyLoss[ i ] = SourceSink ;
	 		var->JouleHeating[iSpecies][i] = JdotE*var->Qe ;
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		C[ 0 ] = 1.0 ;
	 		var->JouleHeating[iSpecies][i] = 0.0 ;
	 		Source = 0.0 ;
	 	}
	 // 	NormalizeCoeff[ i ] = C[ 0 ] ;
	 // 	for( int k = 0 ; k < 5 ; k++ ) C[ k ] /= NormalizeCoeff[i] ;
		// Source /= NormalizeCoeff[i] ;
	
	 	s.push_matrix( row, ncol, col, C ) ;
	 	s.push_source( row, Source) ;
	}//Cell Loop
}
void CEnergyDensity::Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0 ;
	double Diff = 0.0, Mobi=0.0, SourceSink=0.0, JdotE=0.0,TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	for( int i = 0 ; i < nCell ; i++ ) {

		iFace 	 = m->cell[ i ].face_number ;
		iCell 	 = m->cell[ i ].cell_number ;

		/*--- Reset  ---*/
		row 	 = m->cell[ i ].id ;
		col[ 0 ] = row ;
		ncol 	 = 1 ;

		Source 	 = 0.0 ;

		for( int k = 0 ; k < 5 ; k++ ) C[ k ] = 0.0 ;
		/*--- Loop over PLASMA cells ---*/
		if ( m->cell[ i ].type == PLASMA ){

			/*--- Time ---*/
			C[ 0 ] = m->cell[ i ].volume/var->Dt ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = m->Cell[ i ][ k ].NeighborCellId ;

				if ( m->cell[ j ].type == PLASMA ){

					/*--- S-G ---*/
					dL = m->Cell[ i ][ k ].dNPf / m->Cell[ i ][ k ].dDist ;
					dR = m->Cell[ i ][ k ].dPPf / m->Cell[ i ][ k ].dDist ;

					/*--- S-G ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ j ] * var->Mobi[iSpecies][ j ]  ;

					vn = U*m->Cell[ i ][ k ].nf[ 0 ] + V*m->Cell[ i ][ k ].nf[ 1 ] ;

					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ; 

					Pe = vn*m->Cell[ i ][ k ].dDist/Diff ;

					col[ ncol ] = m->Cell[ i ][ k ].NeighborGlobalCellId ;

					if ( Pe < -ZERO  ){

						C[ 0 ] += C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->Cell[ i ][ k ].dArea ;
						C[ncol] = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->Cell[ i ][ k ].dArea ;

					} else if ( Pe > ZERO  ){

						C[ 0 ] += C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->Cell[ i ][ k ].dArea ;
						C[ncol] = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->Cell[ i ][ k ].dArea ;

					} else {

						Diff = -(1.0)*( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;  
						C[ 0 ] 	  += -C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist ;
						C[ncol]    =  C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist ;

					}

		 			ncol ++ ;

	 			} else {//Discontinue face

					Diff 	= -var->Diff[iSpecies][ i ] ;  
					C[ 0 ] += -C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist ;

	 			}//For discontuity face
	 		}//End bulk face

			/*--- Loop over boundary faces ---*/
	 		for( int k = iCell ; k < iFace ; k++ ) {

	 			if( m->cell[ i ].face[ k ]->type == NEUMANN ){
	 				//do nothing
	 			} else {

					Diff 	= -var->Diff[iSpecies][ i ] ;  
					C[ 0 ] += -C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist ;
	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		Source += var->PreU4[iSpecies][ i ] * m->cell[ i ].volume/var->Dt ;

			/*--- Joule heating, Note: since d_Te is in eV, there is no need to multiply the elementary charge ---*/
			JdotE = config->Species[iSpecies].Charge*( var->EField[ 0 ][ i ]*var->U1[ iSpecies ][ i ] 
	 												 + var->EField[ 1 ][ i ]*var->U2[ iSpecies ][ i ] ) ;

	 		/*--- energy loss term ---*/
	 		// switch ( eLOSS ) {
	 		// 	case 1:
	 				SourceSink = *(var->EnergySourcePoint + (i) ) ;
				// 	break;
				// case 2:
	 		// 		SourceSink = var->eEnergyLossTable.GetValueLog( var->T[iSpecies][ i ] ) * var->TotalNumberDensity[ i ]*var->U0[ iSpecies ][ i ] ;
				// 	break;
				// default:
	 		// 		SourceSink = *(var->EnergySourcePoint + (i) ) ;
		  //   		break;
	 		// }
	 		//cout<<"eLoss: "<<SourceSink<<endl;
	 		//exit(1) ;
	 		//Source += (-SourceSink)*m->cell[ i ].volume*var->Dt ;
	 		Source += (JdotE-SourceSink/var->Ref_ES)*m->cell[ i ].volume ;
	 		var->eEnergyLoss[ i ] = SourceSink ;
	 		//exit(1) ;
	 		var->JouleHeating[iSpecies][i] = JdotE*var->Qe ;
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		C[ 0 ] = 1.0 ;
	 		var->JouleHeating[iSpecies][i] = 0.0 ;
	 		Source = 0.0 ;
	 	}
	
	 	s.push_matrix( row, ncol, col, C ) ;
	 	s.push_source( row, Source) ;
	}//Cell Loop
}
void CEnergyDensity::Bulid_A_B_2nd( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int iFace=0, iCell=0, j=0 ;
	// int nCell = m->local_cell_number ;
	// double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0 ;
	// double Diff = 0.0, SourceSink=0.0, JouleHeating=0.0 ;
	// double P=0.0, N=0.0 ;

	// for( int i = 0 ; i < nCell ; i++ ) {

	// 	iFace 	 = m->cell[ i ].face_number ;
	// 	iCell 	 = m->cell[ i ].cell_number ;

	// 	/*--- Reset  ---*/
	// 	Source 	 = 0.0 ;
	// 	row 	 = m->cell[ i ].id ;
	// 	/*--- Loop over PLASMA cells ---*/
	// 	if ( m->cell[ i ].type == PLASMA ){

	// 		/*--- Loop over bulk faces ---*/
	// 		for ( int k = 0 ; k < iCell ; k++ ){

	// 			j = m->Cell[ i ][ k ].NeighborCellId ;

	// 			if ( m->cell[ j ].type == PLASMA ){

	// 				/*--- S-G ---*/
	// 				U = 0.5*( config->Species[ iSpecies ].Charge * var->EField[ 0 ][ i ] * var->Mobi[iSpecies][ i ] 
	// 						+ config->Species[ iSpecies ].Charge * var->EField[ 0 ][ j ] * var->Mobi[iSpecies][ j ] ) ;

	// 				V = 0.5*( config->Species[ iSpecies ].Charge  * var->EField[ 1 ][ i ] * var->Mobi[iSpecies][ i ] 
	// 						+ config->Species[ iSpecies ].Charge  * var->EField[ 1 ][ j ] * var->Mobi[iSpecies][ j ] ) ;

	// 				vn = U*m->Cell[ i ][ k ].nf[ 0 ] + V*m->Cell[ i ][ k ].nf[ 1 ] ;

	// 				Diff = 0.5*( var->Diff[iSpecies][ i ] + var->Diff[iSpecies][ j ] ) ; 

	// 				Pe = vn*m->Cell[ i ][ k ].dDist/Diff ;

	// 				if ( Pe < -ZERO  ){

	// 					P = C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;
	// 					N = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;

	// 				} else if ( Pe > ZERO  ){

	// 					P = C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;
	// 					N = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->Cell[ i ][ k ].dArea*var->Dt ;

	// 				} else {

	// 					Diff = - 0.5*( var->Diff[iSpecies][ i ] + var->Diff[iSpecies][ j ] ) ;  
	// 					P = -C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist*var->Dt ;
	// 					N =  C53*Diff*m->Cell[ i ][ k ].dArea/m->Cell[ i ][ k ].dDist*var->Dt ;
	// 				}

	// 				GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
	// 				GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;

	// 				GVarN[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ j ] ;
	// 				GVarN[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ j ] ;

	// 				Source += (-1.0)*( DotProduct( GVarN, m->Cell[ i ][ k ].NNP)*N 
	// 								 + DotProduct( GVarP, m->Cell[ i ][ k ].PPP)*P ) ;


	//  			} else {/*--- For discontuity face ---*/

	//  				switch ( config->Species[ iSpecies ].Type ){

	// 					/*--- Electron, thermal flux ---*/
	// 					case 0:
	//  						/*--- Drift term ---*/
	// 						U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 0 ][ i ] ;
	//  						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 1 ][ i ] ;
	//  						vn = C43*max( 0.0, U*m->Cell[ i ][ k ].nf[ 0 ]+V*m->Cell[ i ][ k ].nf[ 1 ] ) ;

	//  						/*--- Thermal flux term ---*/
	//  						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	//  						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) ) ;
	 						
	//  						/*--- Secondary electron emission ---*/
	//  						SecondaryElectronEmission = 0.0 ;
	//  						Te_sec = 0.5 ;
	//  						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	//  							if (config->Species[ jSpecies ].Type == ION ){
	//  								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 0 ][ i ]*m->Cell[ i ][ k ].nf[ 0 ] 
	//  												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 1 ][ i ]*m->Cell[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	//  								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	//  							}
	//  						}
	//  						Source += 2.0*Te_sec*SecondaryElectronEmission*m->Cell[ i ][ k ].dArea*var->Dt ;

	//  						P =  vn*m->Cell[ i ][ k ].dArea*var->Dt ;

	//  						GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
	// 						GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;

	// 						Source += (-1.0)*( 0.0 + DotProduct( GVarP, m->Cell[ i ][ k ].PPP)*P ) ;
	 						
	// 					break;

	// 					/*--- Ion, drift flux ---*/	
	// 					case 1:
	// 						cout<<"Not Support Ion energy yet"<<endl;exit(1) ;
	// 					break;

	// 					/*--- Neutral, Diffusion flux ---*/	
	// 					case 2:	
	// 						cout<<"Not Support Neutral energy yet"<<endl;exit(1) ;
	// 					break;

	// 					default:
	// 						if( mpi_id == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
	// 						exit(1);
	// 		    		break;
	// 				}//End switch
	//  			}//For discontuity face
	//  		}//End bulk face

	// 		/*--- Loop over boundary faces ---*/
	//  		for( int k = iCell ; k < iFace ; k++ ) {

	//  			if( m->cell[ i ].face[ k ]->type == NEUMANN ){
	//  				//do nothing
	//  			}else{

	// 				switch ( config->Species[ iSpecies ].Type ){

	// 					/*--- Electron, thermal flux ---*/
	// 					case 0:
	//  						/*--- Drift term ---*/
	// 						U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 0 ][ i ] ;
	//  						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 1 ][ i ] ;
	//  						vn = C43*max( 0.0, U*m->Cell[ i ][ k ].nf[ 0 ]+V*m->Cell[ i ][ k ].nf[ 1 ] ) ;

	//  						/*--- Thermal flux term ---*/
	//  						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	//  						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) ) ;
	 						
	//  						/*--- Secondary electron emission ---*/
	//  						SecondaryElectronEmission = 0.0 ;
	//  						Te_sec = 0.5 ;
	//  						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	//  							if (config->Species[ jSpecies ].Type == ION ){
	//  								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 0 ][ i ]*m->Cell[ i ][ k ].nf[ 0 ] 
	//  												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 1 ][ i ]*m->Cell[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	//  								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	//  							}
	//  						}
	//  						Source += 2.0*Te_sec*SecondaryElectronEmission*m->Cell[ i ][ k ].dArea*var->Dt ;

	//  						P =  vn*m->Cell[ i ][ k ].dArea*var->Dt ;

	//  						GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
	// 						GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;

	// 						Source += (-1.0)*( 0.0 + DotProduct( GVarP, m->Cell[ i ][ k ].PPP)*P ) ;
	 						
	// 					break;

	// 					/*--- Ion, drift flux ---*/	
	// 					case 1:
	// 						cout<<"Not Support Ion energy yet"<<endl;exit(1) ;
	// 					break;

	// 					/*--- Neutral, Diffusion flux ---*/	
	// 					case 2:	
	// 						cout<<"Not Support Neutral energy yet"<<endl;exit(1) ;
	// 					break;

	// 					default:
	// 						if( mpi_id == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
	// 						exit(1);
	// 		    		break;
	// 				}//End switch
	//  			}
	//  		}

	//  		/*--- Previous solution ---*/
	//  		Source += var->PreU4[iSpecies][ i ] * m->cell[ i ].volume ;

	// 		/*--- Joule heating, Note: since d_Te is in eV, there is no need to multiply the elementary charge ---*/
	// 		JouleHeating = config->Species[iSpecies].Charge*( var->EField[ 0 ][ i ]*var->U1[ iSpecies ][ i ] 
	//  														+ var->EField[ 1 ][ i ]*var->U2[ iSpecies ][ i ] ) ;

	//  		/*--- energy loss term ---*/
	//  		SourceSink = *(var->EnergySourcePoint + (i) ) ;

	//  		//Source += (-SourceSink)*m->cell[ i ].volume*var->Dt ;
	//  		Source += (JouleHeating-SourceSink)*m->cell[ i ].volume*var->Dt ;

	//  	/*--- Loop over SOLID cells ---*/
	//  	} else {
	//  		Source = 0.0 ;
	//  	}
	//  	// Source /= NormalizeCoeff[ i ] ;
	//  	s.push_source( row, Source) ;
	// }//Cell Loop
}
void CEnergyDensity::Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	int nCell = m->local_cell_number ;
	for ( int i = 0 ; i < nCell ; i++ ) {
		var->GradU4[ iSpecies ][ 0 ][ i ] = 0.0 ;
		var->GradU4[ iSpecies ][ 1 ][ i ] = 0.0 ;
	}//Loop over all cells
	var->GradU4[ iSpecies ][ 0 ] = var->GradU4[ iSpecies ][ 0 ] ;
	var->GradU4[ iSpecies ][ 1 ] = var->GradU4[ iSpecies ][ 1 ] ;
}
void CEnergyDensity::Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	int nCell = m->local_cell_number ;
	int iFace=0, iCell=0, j=0, NeighborCellIndex=0 ;
	double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0 ;

	for ( int i = 0 ; i < nCell ; i++ ) {

		Gx = 0.0 ; Gy = 0.0 ;
		iCell = m->cell[ i ].cell_number ;  
		iFace = m->cell[ i ].face_number ; 
		
		if ( m->cell[ i ].type != PLASMA ){

			var->GradU4[ iSpecies ][ 0 ][ i ] = 0.0 ;
			var->GradU4[ iSpecies ][ 1 ][ i ] = 0.0 ;

		} else {

			/*--- Loop over neighbor "faces" ---*/
			for ( int k = 0 ; k < iCell ; k++ ) {

				j = m->Cell[ i ][ k ].NeighborCellId ;

				if ( m->cell[ j ].type != PLASMA ) {//For discontinued face, apply neumann

					GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
					GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;
					BC_Value = var->U4[iSpecies][ i ] + DotProduct( GVarP, m->Cell[ i ][ k ].PPP)  ;
					dVar = BC_Value - var->U4[iSpecies][ i ] ;

				}else{
					dVar = var->U4[iSpecies][ j ] - var->U4[iSpecies][ i ] ;
				}
				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	     		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

			}

			/*--- Loop over boundary faces ---*/
			for ( int k = iCell ; k < iFace ; k++ ) {

				// if 		( m->cell[ i ].face[ k ]->type ==   POWER ) BC_Value = Power_voltage ;
				// else if ( m->cell[ i ].face[ k ]->type ==  GROUND )	BC_Value = 0.0 ;
				// else if ( m->cell[ i ].face[ k ]->type == NEUMANN )	{
				GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
				GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;
				BC_Value = var->U4[iSpecies][ i ] + DotProduct( GVarP, m->Cell[ i ][ k ].PPP)  ;
				//} else cout<<"error"<<endl ;
				
				dVar = BC_Value - var->U4[iSpecies][ i ] ;

				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	      		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;


			}
			var->GradU4[ iSpecies ][ 0 ][ i ] = Gx ;
			var->GradU4[ iSpecies ][ 1 ][ i ] = Gy ;
		}
	}//Loop over all cells
	var->GradU4[ iSpecies ][0] = var->GradU4[ iSpecies ][0] ;
	var->GradU4[ iSpecies ][1] = var->GradU4[ iSpecies ][1] ;
}
void CEnergyDensity::CalculateTemperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	int nCell = m->local_cell_number  ;
	double val=0.0, C23=2.0/3.0 ;

	for ( int i = 0 ; i < nCell ; i++ ) {
		
		if ( m->cell[ i ].type != PLASMA ){

			var->T[iSpecies][ i ] = 0.0 ;

		} else {

			val = var->U4[ iSpecies ][ i ]/var->U0[ iSpecies ][ i ]*C23 ;
			if ( val < 0.05 ){

				val = 0.05 ;
				//cout<<"zero"<<endl;
//			}else if( isnan(val) ){
//				val = 1.e-5 ;
			
			} else if( val > 50.0 ){
				val = 50.0 ;
			}
			//if(var->U0[ iSpecies ][ i ] == 0.0 ) {cout<<"N: Zero"<<endl; exit(1) ;}
			var->T[iSpecies][ i ] = val ;
			// cout<<"i: "<<i<<"\t"<<setprecision(15)<<"U4: "<<var->U4[ iSpecies ][ i ]<<endl;
			// cout<<"i: "<<i<<"\t"<<setprecision(15)<<"Te: "<<var->T[iSpecies][ i ]<<endl;
			// cout<<endl;
		}
	}//Loop over all cells
	
	var->T[iSpecies] = var->T[iSpecies] ;
}
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
	if ( mpi_rank == 0 ){
		cout<<"Creat "<<config->Species[iSpecies].Name<<" continuity solver, index: "<<index<<", charge: "<<config->Species[ index ].Charge<<endl;
	} 
	
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

	iMatrix = 1 + config->TotalSpeciesNum  + 1 ;
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

		plasma.get_solution( iMatrix, variable->U4[iSpecies].data ) ;
		variable->U4[iSpecies] = variable->U4[iSpecies] ;

	/*--- calculate electron temperature ---*/
		CalculateTemperature( m, variable ) ;
}
void CEnergyDensity::Bulid_A_B_1st_default( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0 ;
	double Diff = 0.0, Mobi=0.0, SourceSink=0.0, JdotE=0.0,TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	Cell *Cell_i, *Cell_j ;
	plasma.before_matrix_construction( iMatrix ) ;
	plasma.before_source_term_construction( iMatrix ) ;

	for( int i = 0 ; i < nCell ; i++ ) {

		Cell_i = plasma.get_cell( i ) ;
		iFace  = Cell_i->face_number ;
		iCell  = Cell_i->cell_number ;

		/*--- Loop over PLASMA cells ---*/
		if ( Cell_i->type == PLASMA ){

			/*--- Unsteady term ---*/
			plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, Cell_i->volume ) ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = Cell_i->cell[ k ]->local_id ;
				Cell_j = plasma.get_cell( j ) ;


				if ( Cell_j->type == PLASMA ){

					/*--- S-G ---*/
					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;
					/*--- S-G ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ j ] * var->Mobi[iSpecies][ j ]  ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					
					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ; 

					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					
					if ( Pe < -ZERO  ){

						//C[ 0 ] += C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						//C[ncol] = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_j->id, C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;

					} else if ( Pe > ZERO  ){

						//C[ 0 ] += C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						//C[ncol] = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_j->id, C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;

					} else {

						Diff = - ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;  
						//C[ 0 ] 	  += -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ;
						//C[ncol]    =  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ) ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_j->id,  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ) ;

					}

	 			} else {//Discontinue face

	 				switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 0 ][ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 1 ][ i ] ;
	 						vn = C43*max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ]+V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	 						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->EField[ 0 ][ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;
	 						
	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						Te_sec = config->SecondaryElectronEmissionEnergy ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 						}

	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
	 						//Source += 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						plasma.add_entry_in_source_term(iMatrix, i, 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
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
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}//For discontuity face
	 		}//End bulk face

			/*--- Loop over boundary faces ---*/
	 		for( int k = iCell ; k < iFace ; k++ ) {

	 			if( Cell_i->face[ k ]->type == NEUMANN ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 0 ][ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->EField[ 1 ][ i ] ;
	 						vn = C43*max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ]+V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	 						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->EField[ 0 ][ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;
	 						
	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						Te_sec = config->SecondaryElectronEmissionEnergy ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 0 ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->EField[ 1 ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 						}

	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
	 						//Source += 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						plasma.add_entry_in_source_term(iMatrix, i, 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
	 						
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
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		//Source += var->PreU4[iSpecies][ i ] * Cell_i->volume ;
	 		plasma.add_entry_in_source_term(iMatrix, i,  var->PreU4[iSpecies][ i ] * Cell_i->volume ) ;

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
	 		//Source += (-SourceSink)*Cell_i->volume*var->Dt ;
	 		//Source += (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume*var->Dt ;
	 		plasma.add_entry_in_source_term(iMatrix, i, (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume*var->Dt ) ;

	 		var->eEnergyLoss[ i ] = SourceSink ;
	 		var->JouleHeating[iSpecies][i] = JdotE*var->Qe ;
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		//C[ 0 ] = 1.0 ;
	 		plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, 1.0 ) ;
	 		var->JouleHeating[iSpecies][i] = 0.0 ;
	 		//Source = 0.0 ;
	 	}
	}//Cell Loop
	plasma.finish_matrix_construction( iMatrix ) ;
	plasma.finish_source_term_construction( iMatrix ) ;
		
	MPI_Barrier(MPI_COMM_WORLD) ;
}
void CEnergyDensity::Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0 ;
	double Diff = 0.0, Mobi=0.0, SourceSink=0.0, JdotE=0.0,TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	Cell *Cell_i, *Cell_j ;
	plasma.before_matrix_construction( iMatrix ) ;
	plasma.before_source_term_construction( iMatrix ) ;

	for( int i = 0 ; i < nCell ; i++ ) {

		Cell_i = plasma.get_cell(i);

		iFace 	 = Cell_i->face_number ;
		iCell 	 = Cell_i->cell_number ;

		/*--- Loop over PLASMA cells ---*/
		if ( Cell_i->type == PLASMA ){

			/*--- Unsteady term ---*/
			//C[ 0 ] = Cell_i->volume/var->Dt ;
			plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, Cell_i->volume/var->Dt ) ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = Cell_i->cell[k]->local_id ;

				Cell_j = plasma.get_cell(j) ;

				if ( Cell_j->type == PLASMA ){

					/*--- S-G ---*/
					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					/*--- S-G ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 0 ][ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->EField[ 1 ][ j ] * var->Mobi[iSpecies][ j ]  ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ; 

					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					
					if ( Pe < -ZERO  ){

						//C[ 0 ] += C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						//C[ncol] = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_j->id, C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

					} else if ( Pe > ZERO  ){

						//C[ 0 ] += C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						//C[ncol] = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_j->id, C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

					} else {

						Diff = -(1.0)*( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;  
						//C[ 0 ] 	  += -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
						//C[ncol]    =  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
						plasma.add_entry_in_matrix( iMatrix, i,  Cell_j->id,  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
					}

	 			} else {//Discontinue face

					Diff 	= -var->Diff[iSpecies][ i ] ;  
					//C[ 0 ] += -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
					plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
	 			}//For discontuity face
	 		}//End bulk face

			/*--- Loop over boundary faces ---*/
	 		for( int k = iCell ; k < iFace ; k++ ) {

	 			if( Cell_i->face[ k ]->type == NEUMANN ){
	 				//do nothing
	 			} else {

					Diff 	= -var->Diff[iSpecies][ i ] ;  
					//C[ 0 ] += -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
					plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id,-C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		//Source += var->PreU4[iSpecies][ i ] * Cell_i->volume/var->Dt ;
	 		plasma.add_entry_in_source_term( iSpecies+1, i, var->PreU4[iSpecies][ i ] * Cell_i->volume/var->Dt ) ;
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
	 		//Source += (-SourceSink)*Cell_i->volume*var->Dt ;
	 		//Source += (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume ;
	 		plasma.add_entry_in_source_term( iSpecies+1, i, (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume ) ;
	 		var->eEnergyLoss[ i ] = SourceSink ;
	 		var->JouleHeating[iSpecies][i] = JdotE*var->Qe ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		plasma.add_entry_in_matrix( iMatrix, i,  Cell_i->id, 1.0 ) ;
	 		var->JouleHeating[iSpecies][i] = 0.0 ;
	 	}
	}//Cell Loop
	plasma.finish_matrix_construction( iMatrix ) ;
	plasma.finish_source_term_construction( iMatrix ) ;
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
	Cell *Cell_i, *Cell_j ;

	for ( int i = 0 ; i < nCell ; i++ ) {

		Cell_i = plasma.get_cell( i ) ;
		iFace 	 = Cell_i->face_number ;
		iCell 	 = Cell_i->cell_number ;
		
		Gx = 0.0 ; Gy = 0.0 ;

		if ( Cell_i->type != PLASMA ){

			var->GradU4[ iSpecies ][ 0 ][ i ] = 0.0 ;
			var->GradU4[ iSpecies ][ 1 ][ i ] = 0.0 ;

		} else {

			/*--- Loop over neighbor "faces" ---*/
			for ( int k = 0 ; k < iCell ; k++ ) {

				j = Cell_i->cell[k]->local_id ;

				if ( Cell_j->type != PLASMA ) {//For discontinued face, apply neumann

					GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
					GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;
					BC_Value = var->U4[iSpecies][ i ] + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)  ;
					dVar = BC_Value - var->U4[iSpecies][ i ] ;

				}else{
					dVar = var->U4[iSpecies][ j ] - var->U4[iSpecies][ i ] ;
				}
				Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	     	Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

			}

			/*--- Loop over boundary faces ---*/
			for ( int k = iCell ; k < iFace ; k++ ) {

				// if 		( Cell_i->face[ k ]->type ==   POWER ) BC_Value = Power_voltage ;
				// else if ( Cell_i->face[ k ]->type ==  GROUND )	BC_Value = 0.0 ;
				// else if ( Cell_i->face[ k ]->type == NEUMANN )	{
				GVarP[ 0 ] = var->GradU4[ iSpecies ][ 0 ][ i ] ;
				GVarP[ 1 ] = var->GradU4[ iSpecies ][ 1 ][ i ] ;
				BC_Value = var->U4[iSpecies][ i ] + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)  ;
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

	Cell *Cell_i ;

	for ( int i = 0 ; i < nCell ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		if ( Cell_i->type != PLASMA ){

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
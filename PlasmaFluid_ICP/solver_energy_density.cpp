#include "solver_energy_density.hpp"
#define Debug_EE_Bulid_A_B_1st_zero 0
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
		cout<<"Creat "<<config->Species[iSpecies].Name<<" energy density solver, index: "<<index<<", charge: "<<config->Species[ index ].Charge<<endl;
	} 
	
	Correction = config->Equation[ config->Species[ index ].Type ].Correction ;
	WallType   = config->Equation[ config->Species[ index ].Type ].WallBoundaryType ;
	
	fixTe = false ;
	//TG = true ;
	TG = false ;
	Insert = false ;
	eLOSS = 1 ;//1: from chemistry module, 2: from table.
	Reflec = config->ReflectionCoefficient ;
	its=0 ;

	energy_density.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;
	energy_density.load_mesh( &plasma) ;
}
void CEnergyDensity::Solver( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{
	//Zero_Gradient( m, variable ) ;

	switch( WallType ){

		case 0: //default
			Bulid_A_B_1st_default( m, config, variable ) ;
			break ;

		case 1://neumann
			cout<<"CEnergyDensity-Neumann, not support yet."<<endl;
			exit(1) ;
			break;

		case 2://zero number density
			Bulid_A_B_1st_zero( m, config, variable ) ;
			break;
		default:
			break;

	}//Wall boundary type.

	energy_density.get_solution( variable->U4[iSpecies].data ) ;
	variable->U4[iSpecies] = variable->U4[iSpecies] ;
	its = energy_density.get_iteration_number() ;
	/*--- calculate electron temperature ---*/
		CalculateTemperature( m, variable ) ;
		if( iSpecies == ELECTRON ) CalculatePowerAbs( m, variable  );
}
void CEnergyDensity::Bulid_A_B_1st_default( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;

	double vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0 ;
	double Diff = 0.0, Mobi=0.0, SourceSink=0.0, JdotE=0.0,TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	Cell *Cell_i, *Cell_j ;
	energy_density.before_matrix_construction() ;
	energy_density.before_source_term_construction() ;

	for( int i = 0 ; i < energy_density.Mesh.cell_number ; i++ ) {

		Cell_i = energy_density.get_cell( i ) ;
		iFace  = Cell_i->face_number ;
		iCell  = Cell_i->cell_number ;

		/*--- Loop over PLASMA cells ---*/
		if (  energy_density.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Unsteady term ---*/
			energy_density.add_entry_in_matrix( i,  Cell_i->id, Cell_i->volume ) ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = Cell_i->cell[ k ]->local_id ;
				Cell_j = energy_density.get_cell( j ) ;


				if ( energy_density.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					/*--- S-G ---*/
					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;
					/*--- S-G ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ]  ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					
					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ; 

					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					
					if ( Pe < -ZERO  ){

						//C[ 0 ] += C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						//C[ncol] = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						energy_density.add_entry_in_matrix( i,  Cell_i->id, C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
						energy_density.add_entry_in_matrix( i,  Cell_j->id, C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;

					} else if ( Pe > ZERO  ){

						//C[ 0 ] += C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						//C[ncol] = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
						energy_density.add_entry_in_matrix( i,  Cell_i->id, C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
						energy_density.add_entry_in_matrix( i,  Cell_j->id, C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;

					} else {

						Diff = - ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;  
						//C[ 0 ] 	  += -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ;
						//C[ncol]    =  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ;
						energy_density.add_entry_in_matrix( i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ) ;
						energy_density.add_entry_in_matrix( i,  Cell_j->id,  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist*var->Dt ) ;

					}

	 			} else {//Discontinue face

	 				switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->Ex[ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->Ey[ i ] ;
	 						vn = C43*max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ]+V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	 						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->Ex[ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;
	 						
	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						Te_sec = config->SecondaryElectronEmissionEnergy ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 						}

	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						energy_density.add_entry_in_matrix( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
	 						//Source += 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						energy_density.add_entry_in_source_term( i, 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
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

	 			if( energy_density.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->Ex[ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ] * var->Ey[ i ] ;
	 						vn = C43*max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ]+V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ; if( fixTe ) Te = 0.5 ;
	 						vn += C43*0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->Ex[ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;
	 						
	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						Te_sec = config->SecondaryElectronEmissionEnergy ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){
	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ; ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 						}

	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						energy_density.add_entry_in_matrix( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
	 						//Source += 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ;
	 						energy_density.add_entry_in_source_term( i, 2.0*Te_sec*SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea*var->Dt ) ;
	 						
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
	 		energy_density.add_entry_in_source_term( i,  var->PreU4[iSpecies][ i ] * Cell_i->volume ) ;

			/*--- Joule heating, Note: since d_Te is in eV, there is no need to multiply the elementary charge ---*/
			JdotE = config->Species[iSpecies].Charge*( var->Ex[ i ]*var->U1[ iSpecies ][ i ] 
	 												 										 + var->Ey[ i ]*var->U2[ iSpecies ][ i ] ) ;

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
	 		energy_density.add_entry_in_source_term( i, (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume*var->Dt ) ;

	 		var->eEnergyLoss[ i ] = SourceSink ;
	 		var->JouleHeating[iSpecies][i] = JdotE*var->Qe ;
			
			//For ICP power absorption
			energy_density.add_entry_in_source_term( i, var->Power_Absorption_plasma[ i ]*Cell_i->volume ) ;
			
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		//C[ 0 ] = 1.0 ;
	 		energy_density.add_entry_in_matrix( i,  Cell_i->id, 1.0 ) ;
	 		var->JouleHeating[iSpecies][i] = 0.0 ;
	 		//Source = 0.0 ;
	 	}
	}//Cell Loop
	energy_density.finish_matrix_construction() ;
	energy_density.finish_source_term_construction() ;

	MPI_Barrier(MPI_COMM_WORLD) ;
}
void CEnergyDensity::Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, iCell=0, j=0 ;

	double vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0, Te_sec=0.0, Source=0.0 ;
	double Diff = 0.0, Mobi=0.0, SourceSink=0.0, JdotE=0.0,TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	Cell *Cell_i, *Cell_j ;
	energy_density.before_matrix_construction() ;
	energy_density.before_source_term_construction() ;

	for( int i = 0 ; i < energy_density.Mesh.cell_number ; i++ ) {

		Cell_i = energy_density.get_cell( i ) ;
		iFace  = Cell_i->face_number ;
		iCell  = Cell_i->cell_number ;

		/*--- Loop over PLASMA cells ---*/
		if (  energy_density.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Unsteady term ---*/
			energy_density.add_entry_in_matrix( i,  Cell_i->id, Cell_i->volume/var->Dt ) ;

			#if Debug_EE_Bulid_A_B_1st_zero
				C[ 0 ] = Cell_i->volume/var->Dt ;
			#endif

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = Cell_i->cell[ k ]->local_id ;
				Cell_j = energy_density.get_cell( j ) ;


				if ( energy_density.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					/*--- S-G ---*/
					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;
					/*--- S-G ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ]  ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					
					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ; 

					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					
					if ( Pe < -ZERO  ){

						energy_density.add_entry_in_matrix( i,  Cell_i->id, C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						energy_density.add_entry_in_matrix( i,  Cell_j->id, C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

						#if Debug_EE_Bulid_A_B_1st_zero
							C[ 0 ] += C53*vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
							C[ncol] = C53*vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						#endif

					} else if ( Pe > ZERO  ){

						energy_density.add_entry_in_matrix( i,  Cell_i->id, C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						energy_density.add_entry_in_matrix( i,  Cell_j->id, C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						#if Debug_EE_Bulid_A_B_1st_zero
							C[ 0 ] += C53*vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
							C[ncol] = C53*vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						#endif
					} else {

						Diff = - ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;  
						energy_density.add_entry_in_matrix( i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
						energy_density.add_entry_in_matrix( i,  Cell_j->id,  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
						#if Debug_EE_Bulid_A_B_1st_zero						
							C[ 0 ] += -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
							C[ncol] =  C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
						#endif
					}
					ncol ++ ;

	 			} else {//Discontinue face

	 				switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:

							Diff = -var->Diff[iSpecies][ i ] ;
							energy_density.add_entry_in_matrix( i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
							#if Debug_EE_Bulid_A_B_1st_zero
							C[0]+= -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
							#endif
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

	 			if( energy_density.get_face_typename( Cell_i->face[ k ]->data_id )  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						/*--- Electron, thermal flux ---*/
						case 0:
							Diff = -var->Diff[iSpecies][ i ] ;
							energy_density.add_entry_in_matrix( i,  Cell_i->id, -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
							#if Debug_EE_Bulid_A_B_1st_zero
								C[0]+= -C53*Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
							#endif
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
	 		energy_density.add_entry_in_source_term( i,  var->PreU4[iSpecies][ i ] * Cell_i->volume/var->Dt ) ;
	 		#if Debug_EE_Bulid_A_B_1st_zero
				Source +=  var->PreU4[iSpecies][ i ] * Cell_i->volume/var->Dt ;
	 		#endif
			/*--- Joule heating, Note: since d_Te is in eV, there is no need to multiply the elementary charge ---*/
			JdotE = config->Species[iSpecies].Charge*( 
														 var->Ex[ i ]*var->U1[ iSpecies ][ i ] //);
	 												 + var->Ey[ i ]*var->U2[ iSpecies ][ i ] ) ;

	 		/*--- energy loss term ---*/
			SourceSink = *(var->EnergySourcePoint + (i) ) ;

	 		energy_density.add_entry_in_source_term( i, (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume ) ;

	 		#if Debug_EE_Bulid_A_B_1st_zero
	 			Source +=  (JdotE-SourceSink/var->Ref_ES)*Cell_i->volume ;
	 		#endif


	 		var->eEnergyLoss[ i ] = SourceSink ;
	 		var->JouleHeating[iSpecies][i] = JdotE*var->Qe ;
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		//C[ 0 ] = 1.0 ;
	 		energy_density.add_entry_in_matrix( i,  Cell_i->id, 1.0 ) ;
	 		#if Debug_EE_Bulid_A_B_1st_zero
	 		C[0] = 1.0 ;
	 		#endif
	 		var->JouleHeating[iSpecies][i] = 0.0 ;
	 		//Source = 0.0 ;
	 	}
	 	#if Debug_EE_Bulid_A_B_1st_zero
			cout<<"i: "<<i<<endl;
			for ( int icol = 0 ; icol < ncol ; icol++ ){
				cout<<scientific<<setprecision(15)<<"C["<<icol<<"]: "<<C[icol]<<endl;
			}
			cout<<scientific<<setprecision(15)<<Source<<endl;
			cout<<"U4: "<<var->PreU4[iSpecies][ i ]<<endl;
			cout<<"U1: "<<var->PreU0[iSpecies][ i ]<<endl;
			cout<<"Te: "<<var->PreU4[iSpecies][ i ]/var->PreU0[iSpecies][ i ]*2.0/3.0<<endl;
			cout<<"T2: "<<var->PreT[iSpecies][ i ]<<endl;
			cout<<endl;
		#endif
	}//Cell Loop
	energy_density.finish_matrix_construction() ;
	energy_density.finish_source_term_construction() ;
	MPI_Barrier(MPI_COMM_WORLD) ;
}
void CEnergyDensity::CalculateTemperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{

	double val=0.0, C23=2.0/3.0 ;

	Cell *Cell_i ;

	for ( int i = 0 ; i < energy_density.Mesh.cell_number ; i++ ) {

		Cell_i = energy_density.get_cell(i) ;

		if (  energy_density.get_cell_typename( Cell_i->data_id ) != "PLASMA" ){

			var->T[iSpecies][ i ] = 0.0 ;

		} else {

			val = var->U4[ iSpecies ][ i ]/var->U0[ iSpecies ][ i ]*C23 ;
			if ( val < 0.05 ){

				val = 0.05 ;

			} else if( val > 50.0 ) {

				val = 50.0 ;
				
			}
			//if(var->U0[ iSpecies ][ i ] == 0.0 ) {cout<<"N: Zero"<<endl; exit(1) ;}
			var->T[iSpecies][ i ] = val ;
			//cout<<"i: "<<i<<", Te_diff: "<<var->T[iSpecies][ i ]-val<<endl;
			// cout<<"i: "<<i<<"\t"<<setprecision(15)<<"U4: "<<var->U4[ iSpecies ][ i ]<<endl;
			// cout<<"i: "<<i<<"\t"<<setprecision(15)<<"Te: "<<var->T[iSpecies][ i ]<<endl;
			// cout<<endl;
		}
	}//Loop over all cells
	
	var->T[iSpecies] = var->T[iSpecies] ;
}
void CEnergyDensity::CalculatePowerAbs( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{

	double val=0.0, C23=2.0/3.0 ;

	Cell *Cell_i ;
	double power_local=0.0 ;
	var->PowerAbs = 0.0 ;

	for ( int i = 0 ; i < energy_density.Mesh.cell_number ; i++ ) {
		Cell_i = energy_density.get_cell(i) ;
		if (  energy_density.get_cell_typename( Cell_i->data_id ) != "PLASMA" ){
			//no power
		} else {
			power_local += var->JouleHeating[iSpecies][i]*Cell_i->volume ;
		}
	}//Loop over all cells
	var->PowerAbs = energy_density.parallel_sum( &power_local ) ;
	MPI_Barrier(MPI_COMM_WORLD) ;
}

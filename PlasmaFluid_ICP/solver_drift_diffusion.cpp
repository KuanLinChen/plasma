
#include "solver_drift_diffusion.hpp"
#define Debug_Bulid_A_B_1st_zero 0
using namespace std ;
CDriftDiffusion::CDriftDiffusion()
{
}
void CDriftDiffusion::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, int index )
{
	iSpecies = index ;
	SpeciesType = config->Species[ index ].Type ;

	Correction 		 = config->Equation[ SpeciesType ].Correction ;
	WallType 		 = config->Equation[ SpeciesType ].WallBoundaryType ;

	if ( mpi_rank == 0 ){
		cout<<"Creat "<<config->Species[iSpecies].Name<<" continuity solver, index: "<<index<<", charge: "<<config->Species[ index ].Charge<<", change sign: "<<config->Species[ index ].sgn<<", Speciec Type: "<<SpeciesType<<endl;
		cout<<"Correction: "<<Correction<<endl;
	} 
	fixTe = false ;
	TG = false ;
	Insert = false ;
	Reflec = config->ReflectionCoefficient ;
	its=0 ;
	drift_diffusion.initial( gargc2, gargv2, &mpi_rank, &mpi_size ) ;
	drift_diffusion.set_linear_solver_library("PETSC");
	drift_diffusion.load_mesh( &plasma) ;
}
void CDriftDiffusion::Solve_Diffusion( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{

}
void CDriftDiffusion::Solve( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{
	if( config->Species[iSpecies].Activate ==0) return ; 

	/*--- 1st solve the equation w/o cross diffusion term ---*/
		switch( WallType ){
			case 0: //default
				Bulid_A_B_1st_default( m, config, variable ) ;
				break ;

			case 1://neumann
				Bulid_A_B_1st_neumann( m, config, variable ) ;
				break;

			case 2://zero number density
				Bulid_A_B_1st_zero( m, config, variable ) ;
				break;

			case 3://Zero-Dimension
				Bulid_A_B_1st_0D( m, config, variable ) ;
			break;

			default:
				cout<<"wall type"<<WallType<<endl;exit(1);
				break;
		}
		drift_diffusion.get_solution( variable->U0[iSpecies].data ) ;
		its = drift_diffusion.get_iteration_number() ;
		variable->U0[iSpecies] = variable->U0[iSpecies] ;

	/*--- calculate drift-diffustion flux ---*/
		switch( WallType ){
			case 0: //default
				CalculateAvgDDFlux_default( m, config, variable ) ;
				break ;

			case 1://neumann
				CalculateAvgDDFlux_neumann( m, config, variable ) ;
				break;

			case 2://zero number density
				CalculateAvgDDFlux_zero( m, config, variable ) ;
				break;
			default:
				break;
		}
		if ( config->Species[iSpecies].Charge != 0.0 ){
			CalculateCondCurrentDensity ( m, config, variable ) ;
			CalculateSurfaceCharge 		( m, config, variable ) ;
		} 
}
void CDriftDiffusion::Bulid_A_B_1st_default( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0 ;
	double Diff=0.0, Mobi=0.0, SourceSink=0.0, TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;

	Cell *Cell_i, *Cell_j ;
	drift_diffusion.before_matrix_construction() ;
	drift_diffusion.before_source_term_construction() ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {


		Cell_i = drift_diffusion.get_cell(i) ;

		Source 	 = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( cell_type[ Cell_i->type ] == PLASMA ){

			/*--- Unsteady term ---*/
			drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, Cell_i->volume/var->Dt ) ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ; 
				Cell_j = drift_diffusion.get_cell(j) ;

				if ( cell_type[ Cell_j->type ] == PLASMA ){

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ]  ;


					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					Diff = ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );

					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;


					if ( Pe < -ZERO ) {

						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id, vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

					} else if ( Pe > ZERO ) {

						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id, vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

					} else {

						Diff = (-1.0)*( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );
						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id,  Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;

					}

	 			} else {/*--- For discontuity face ---*/

	 				switch ( config->Species[ iSpecies ].Type ){

						case 0:/*--- Electron ---*/

	 						/*--- Drift term ---*/
							U  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] );

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ;	if ( fixTe ) Te = 0.5 ;
	 						vn += 0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->Ex[ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;

	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 							
	 						}
	 						//Source += SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
							drift_diffusion.add_entry_in_matrix     ( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;
							drift_diffusion.add_entry_in_source_term( i, SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea ) ;
						break;

						case 1:/*--- Ion ---*/

							/*--- Drift term ---*/
							U = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

							/*--- Thermal flux term ---*/
	 						//vn += 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;
	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
	 						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;
						break;

						/*--- Neutral, Diffusion flux ---*/	
						case 2:	
							Diff = -var->Diff[iSpecies][ i ] ;
	  					//C[ 0 ] += -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dPPf ;//*var->Dt ;
	  					drift_diffusion.add_entry_in_matrix     ( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dPPf ) ;

	  					/*--- Thermal flux term ---*/
	 						vn = 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;
	 						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;

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

	 			if( face_type[ Cell_i->face[ k ]->type ]  == NEUMANN ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						case 0:/*--- Electron ---*/

	 						/*--- Drift term ---*/
							U  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] );

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ;	if ( fixTe ) Te = 0.5 ;
	 						vn += 0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->Ex[ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;

	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 							
	 						} 
	 						
	 						//Source += SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
							drift_diffusion.add_entry_in_matrix     ( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;
							drift_diffusion.add_entry_in_source_term( i, SecondaryElectronEmission*m->PFM_CELL[ i ][ k ].dArea ) ;
						break;

						case 1:/*--- Ion ---*/

							/*--- Drift term ---*/
							U = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

							/*--- Thermal flux term ---*/
	 						//vn += 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;
	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
	 						drift_diffusion.add_entry_in_matrix     ( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;
						break;

						/*--- Neutral, Diffusion flux ---*/	
						case 2:	

							Diff = -var->Diff[iSpecies][ i ] ;
	  					//C[ 0 ] += -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dPPf ;//*var->Dt ;
	  					drift_diffusion.add_entry_in_matrix     ( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dPPf  ) ;
	  						/*--- Thermal flux term ---*/
	 						vn = 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;
	 						//C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea ;//*var->Dt ;
	 						drift_diffusion.add_entry_in_matrix     ( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;
						break;

						default:
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		//Source += (var->PreU0[iSpecies][ i ])*Cell_i->volume/var->Dt ;
	 		drift_diffusion.add_entry_in_source_term( i, (var->PreU0[iSpecies][ i ])*Cell_i->volume/var->Dt ) ;

	 		/*--- Source/Sink term ---*/
	 		if ( config->PFM_Assumption == "LFA" ) {
	 			SourceSink = var->LFASourceSink[ iSpecies ][ i ]/var->Ref_SS ;
	 		} else {
	 			SourceSink = (double)*( var->ReactionRatePoint[iSpecies] + i  )/var->Ref_SS ;
	 		}
	 		var->ProductionRate[iSpecies][ i ] = SourceSink ;

	 		drift_diffusion.add_entry_in_source_term( i, SourceSink*Cell_i->volume ) ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {

	 		drift_diffusion.add_entry_in_matrix( i, Cell_i->id, 1.0 ) ;
	 		var->ProductionRate[iSpecies][ i ] = 0.0 ;
	 		
	 	}//End plasma Cell.
	}//Cell Loop
	
	drift_diffusion.finish_matrix_construction() ;
	drift_diffusion.finish_source_term_construction() ;
		
	MPI_Barrier(MPI_COMM_WORLD) ;

	
}
void CDriftDiffusion::Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0 ;
	double Diff=0.0, Mobi=0.0, SourceSink=0.0, TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;
	Cell *Cell_i, *Cell_j ;


	drift_diffusion.before_matrix_construction() ;
	drift_diffusion.before_source_term_construction() ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell(i);

		#if Debug_Bulid_A_B_1st_zero
			ncol = 1 ;
			Source = 0.0 ;
			for( int k = 0 ; k < 5 ; k++ ) C[ k ] = 0.0 ;	
		#endif

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Unsteady term ---*/
			drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, Cell_i->volume/var->Dt ) ;

			#if Debug_Bulid_A_B_1st_zero
				C[ 0 ] = Cell_i->volume/var->Dt ;
			#endif

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ;

				Cell_j = drift_diffusion.get_cell(j) ;

				if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ]  ;


					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					Diff = ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );
					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;


					if ( Pe < -ZERO ) {

						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id, vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

						#if Debug_Bulid_A_B_1st_zero
							C[ 0 ] += vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
							C[ncol] = vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						#endif

					} else if ( Pe > ZERO ) {

						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id, vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

						#if Debug_Bulid_A_B_1st_zero
							C[ 0 ] += vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
							C[ncol] = vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						#endif

					} else {

						Diff = (-1.0)*( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );
						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id,  Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;	
						#if Debug_Bulid_A_B_1st_zero
							C[ 0 ] += -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
							C[ncol] =  Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
						#endif

					}
					ncol ++ ;

	 			} else {/*--- For discontuity face ---*/

					Diff 	= -var->Diff[iSpecies][ i ] ;
	  			drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
	  			#if Debug_Bulid_A_B_1st_zero
	  				C[ 0 ] +=  vn*m->PFM_CELL[ i ][ k ].dArea ;
	  			#endif
	 			}
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					Diff = -var->Diff[iSpecies][ i ] ;
					drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
					#if Debug_Bulid_A_B_1st_zero
						C[ 0 ] += -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ;
					#endif
	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		drift_diffusion.add_entry_in_source_term( i, var->PreU0[iSpecies][ i ]*Cell_i->volume/var->Dt ) ;
	 		// Source += var->PreU0[iSpecies][ i ]*Cell_i->volume/var->Dt ;

	 		/*--- Source/Sink term ---*/
	 		if ( config->PFM_Assumption == "LFA" ) {

	 			SourceSink = var->LFASourceSink[ iSpecies ][ i ]/var->Ref_SS ;

	 		} else {

	 			SourceSink = (double)*( var->ReactionRatePoint[iSpecies] + i  )/var->Ref_SS ;

	 		}
	 		var->ProductionRate[iSpecies][ i ] = SourceSink ;
	 		drift_diffusion.add_entry_in_source_term( i, SourceSink*Cell_i->volume ) ;
	 		#if Debug_Bulid_A_B_1st_zero
	 			Source += SourceSink*Cell_i->volume ;
	 		#endif

	 	/*--- Loop over SOLID cells ---*/
	 	} else {

	 		drift_diffusion.add_entry_in_matrix( i, Cell_i->id, 1.0 ) ;
	 		#if Debug_Bulid_A_B_1st_zero
	 			C[0] = 1.0 ;
	 		#endif
	 		var->ProductionRate[iSpecies][ i ] = 0.0 ;	
	 	}

	 	#if Debug_Bulid_A_B_1st_zero
			cout<<"i: "<<i<<endl;
			for ( int icol = 0 ; icol < ncol ; icol++ ){
				cout<<scientific<<"C["<<icol<<"]: "<<C[icol]<<endl;
			}
			cout<<Source<<endl;
			cout<<endl;
	 	#endif

	}//Cell Loop
	drift_diffusion.finish_matrix_construction() ;
	drift_diffusion.finish_source_term_construction() ;

	#if Debug_Bulid_A_B_1st_zero
		exit(1) ;
	#endif
}
void CDriftDiffusion::Bulid_A_B_1st_neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double Source=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0, IonFlux=0.0 ;
	double Diff=0.0, Mobi=0.0, SourceSink=0.0, TempGradient=0.0, f1=0.0, f2=0.0, dL=0.0, dR=0.0 ;
	Cell *Cell_i, *Cell_j ;

	drift_diffusion.before_matrix_construction() ;
	drift_diffusion.before_source_term_construction() ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell(i) ;

		Source 	 = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Unsteady term ---*/
			drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, Cell_i->volume/var->Dt ) ;

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ;

				Cell_j = drift_diffusion.get_cell(j) ;

				if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ]  ;


					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					Diff = ( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );

					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;


					if ( Pe < -ZERO ) {

						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id, vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

					} else if ( Pe > ZERO ) {

						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id, vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ) ;

					} else {

						Diff = (-1.0)*( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] );
						drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, -Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;
						drift_diffusion.add_entry_in_matrix( i,  Cell_j->id,  Diff*m->PFM_CELL[ i ][ k ].dArea/m->PFM_CELL[ i ][ k ].dDist ) ;	

					}

	 			} else {/*--- For discontuity face ---*/

 					/*--- Drift term ---*/
					U  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
 					V  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
 					vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] );
 					drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea ) ;

	 			}
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ) {

	 				//do nothing

	 			} else {

					/*--- Drift term ---*/
					U  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
					V  = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
					vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] );
					drift_diffusion.add_entry_in_matrix( i,  Cell_i->id, vn*m->PFM_CELL[ i ][ k ].dArea  ) ;

	 			}
	 		}

	 		/*--- Previous solution ---*/
	 		//Source += (var->PreU0[iSpecies][ i ])*Cell_i->volume/var->Dt ;
	 		drift_diffusion.add_entry_in_source_term( i, (var->PreU0[iSpecies][ i ])*Cell_i->volume/var->Dt ) ;

	 		/*--- Source/Sink term ---*/
	 		if ( config->PFM_Assumption == "LFA" ) {
	 			SourceSink = var->LFASourceSink[ iSpecies ][ i ]/var->Ref_SS ;
	 		} else {
	 			SourceSink = (double)*( var->ReactionRatePoint[iSpecies] + i  )/var->Ref_SS ;
	 		}

	 		var->ProductionRate[iSpecies][ i ] = SourceSink ;
	 		drift_diffusion.add_entry_in_source_term( i, SourceSink*Cell_i->volume ) ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {

	 		drift_diffusion.add_entry_in_matrix( i, Cell_i->id,  1.0 ) ;
	 		var->ProductionRate[iSpecies][ i ] = 0.0 ;
	 		
	 	}
	}//Cell Loop
	drift_diffusion.finish_matrix_construction() ;
	drift_diffusion.finish_source_term_construction() ;
}
void CDriftDiffusion::Bulid_A_B_1st_0D( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double Source=0.0 ;
	double SourceSink=0.0 ;
	Cell *Cell_i, *Cell_j ;

	drift_diffusion.before_matrix_construction() ;
	drift_diffusion.before_source_term_construction() ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		Source 	 = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Unsteady term ---*/
			//C[ 0 ] = Cell_i->volume/var->Dt ;
			drift_diffusion.add_entry_in_matrix( i, Cell_i->id, Cell_i->volume/var->Dt ) ;

	 		/*--- Previous solution ---*/
	 		//Source += (var->PreU0[iSpecies][ i ])*Cell_i->volume/var->Dt ;
	 		drift_diffusion.add_entry_in_source_term( i, (var->PreU0[iSpecies][ i ])*Cell_i->volume/var->Dt ) ;

	 		/*--- Source/Sink term ---*/
	 		if ( config->PFM_Assumption == "LFA" ) {
	 			SourceSink = var->LFASourceSink[ iSpecies ][ i ]/var->Ref_SS ;
	 		} else {
	 			SourceSink = (double)*( var->ReactionRatePoint[iSpecies] + i  )/var->Ref_SS ;
	 		}
	 		var->ProductionRate[iSpecies][ i ] = SourceSink ;

	 		//Source += SourceSink*Cell_i->volume ;//*var->Dt ;
	 		drift_diffusion.add_entry_in_source_term( i, SourceSink*Cell_i->volume ) ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {

	 		drift_diffusion.add_entry_in_matrix( i, Cell_i->id, Cell_i->volume/var->Dt ) ;
	 		var->ProductionRate[iSpecies][ i ] = 0.0 ;
	 		
	 	}
	}//Cell Loop

	drift_diffusion.finish_matrix_construction() ;
	drift_diffusion.finish_source_term_construction() ;
}
void CDriftDiffusion::CalculateAvgDDFlux_default( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double xFlux=0.0, yFlux=0.0, faceFlux=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0 ;
	double P=0.0, N=0.0, PV=0.0, NV=0.0, IonFlux=0.0 ;
	double E = 0.0, Diff=0.0, Mobi=0.0, TempGradient=0.0, dL=0.0, dR=0.0, f1=0.0, f2=0.0 ;

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell( i ) ;

		xFlux = 0.0 ;
		yFlux = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				j = Cell_i->cell[k]->local_id ;
				Cell_j = drift_diffusion.get_cell( j ) ;

				if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					/*--- S-G Scheme ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ] ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ;
					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					/*--- S-G ---*/
					if ( Pe < -ZERO ) {

						P = vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						N = vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;

					} else if ( Pe > ZERO ) {

						P = vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						N = vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;

					} else {

						Diff = -( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;
						P = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;
						N =  Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;

					}
					GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
					GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;

					GVarN[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ j ] ;
					GVarN[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ j ] ;

					PV = var->U0[iSpecies][ i ] ;
					NV = var->U0[iSpecies][ j ] ;

					faceFlux = ( P*PV + N*NV )*m->PFM_CELL[ i ][ k ].dPPf ;

					xFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
					yFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

	 			} else if( drift_diffusion.get_cell_typename( Cell_j->data_id ) != "PLASMA" ) {//Discontinue face

	 				if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ) {
	 					//do nothing
	 				} else {

						switch ( config->Species[ iSpecies ].Type ) {

							case 0:/*--- Electron ---*/

		 						/*--- Drift term ---*/
								U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
		 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
		 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

		 						/*--- Thermal flux term ---*/
		 						Te = var->T[ 0 ][ i ] ;
		 						if ( fixTe ) Te = 0.5 ;
		 						vn += 0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->Ex[ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;

		 						/*--- Secondary electron emission ---*/
		 						SecondaryElectronEmission = 0.0 ;
		 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

		 							if (config->Species[ jSpecies ].Type == ION ){
		 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
		 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ;
		 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
		 							}
		 							
		 						}//if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) != DIELECTRIC ) SecondaryElectronEmission = 0.0 ;


								GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
								GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;

								PV = var->U0[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;


		 						faceFlux = ( vn*PV - SecondaryElectronEmission)*m->PFM_CELL[ i ][ k ].dArea*m->PFM_CELL[ i ][ k ].dPPf ;

		 						xFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
								yFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

							break;

							case 1:/*--- Ion ---*/

								/*--- Drift term ---*/
								U = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
		 						V = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
		 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

								/*--- Thermal flux term ---*/
		 						//vn += 0.25*sqrt( 8.0*var->Qe*var->T[ iSpecies ][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;

		 						GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
								GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;
								PV = var->U0[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

		 						faceFlux=  ( vn*PV )*m->PFM_CELL[ i ][ k ].dArea*m->PFM_CELL[ i ][ k ].dPPf ;

		 						xFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
								yFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

							break;
							
							case 2:/*--- Neutral, Diffusion flux ---*/	

								/*--- Diffusion flux ---*/
								Diff = -var->Diff[iSpecies][ i ] ; 

		  						P = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;

		  					GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
								GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;
								PV = var->U0[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

		  					faceFlux= ( P*PV )*m->PFM_CELL[ i ][ k ].dPPf ;

		 						xFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
								yFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

							break;

							default:
								if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
								exit(1);
				    		break;
						}//End switch
					}
	 			}//For discontuity face
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						case 0:/*--- Electron ---*/

	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ;
	 						if ( fixTe ) Te = 0.5 ;
	 						vn += 0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);//*exp(-fabs(var->Ex[ i ]*0.5*m->PFM_CELL[ i ][ k ].dDist)/var->T[0][i]) ;

	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 							
	 						}//if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) != DIELECTRIC ) SecondaryElectronEmission = 0.0 ;

							GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
							GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;

							PV = var->U0[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;
							//NV = var->U0[iSpecies][ j ] + DotProduct( GVarN, m->PFM_CELL[ i ][ k ].NNP ) ;

	 						faceFlux = ( vn*PV - SecondaryElectronEmission)*m->PFM_CELL[ i ][ k ].dArea*m->PFM_CELL[ i ][ k ].dPPf ;

	 						xFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
							yFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
							// cout<<"i: "<<i<<endl;
							// cout<<"xBc flux: "<<faceFlux <<endl ;//* m->PFM_CELL[ i ][ k ].nf[ 0 ]<<endl ;
							// cout<<"yBc flux: "<<faceFlux <<endl ;//* m->PFM_CELL[ i ][ k ].nf[ 1 ]<<endl<<endl ;

						break;

						case 1:/*--- Ion ---*/

							/*--- Drift term ---*/
							U = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

							/*--- Thermal flux term ---*/
	 						//vn += 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;

	 						GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
							GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;
							PV = var->U0[iSpecies][ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

	 						faceFlux=  ( vn*PV )*m->PFM_CELL[ i ][ k ].dArea*m->PFM_CELL[ i ][ k ].dPPf ;

	 						xFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
							yFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

						break;

						
						case 2:/*--- Neutral ---*/	

							/*--- Diffusion flux ---*/	
							Diff = -var->Diff[iSpecies][ i ] ; 

	  						P = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;

	  						GVarP[ 0 ] = var->GradU0[ iSpecies ][ 0 ][ i ] ;
							GVarP[ 1 ] = var->GradU0[ iSpecies ][ 1 ][ i ] ;
							PV = var->U0[iSpecies][ i ] + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP ) ;

	  						faceFlux= ( P*PV )*m->PFM_CELL[ i ][ k ].dPPf ;

	 						xFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
							yFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

						break;

						default:
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}


	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		xFlux = 0.0 ;
	 		yFlux = 0.0 ;
	 	}
	 	var->U1[ iSpecies ][ i ] = xFlux/Cell_i->volume ;
		var->U2[ iSpecies ][ i ] = yFlux/Cell_i->volume ;

	}//Cell Loop

	/*--- Update ghost cells ---*/
	var->U1[ iSpecies ] = var->U1[ iSpecies ] ;
	var->U2[ iSpecies ] = var->U2[ iSpecies ] ;
}
void CDriftDiffusion::CalculateAvgDDFlux_neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double xFlux=0.0, yFlux=0.0, faceFlux=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0 ;
	double P=0.0, N=0.0, PV=0.0, NV=0.0, IonFlux=0.0 ;
	double E = 0.0, Diff=0.0, Mobi=0.0, TempGradient=0.0, dL=0.0, dR=0.0, f1=0.0, f2=0.0 ;

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell( i ) ;

		xFlux = 0.0 ;
		yFlux = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j      = Cell_i->cell[ k ]->local_id ;
				Cell_j = drift_diffusion.get_cell( j ) ;

				if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					/*--- S-G Scheme ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ] ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ;
					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					/*--- S-G ---*/
					if ( Pe < -ZERO ) {

						P = vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						N = vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;

					} else if ( Pe > ZERO ) {

						P = vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						N = vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;

					} else {

						Diff = -( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;
						P = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;
						N =  Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;

					}

					PV = var->U0[iSpecies][ i ] ;
					NV = var->U0[iSpecies][ j ] ;

					faceFlux = ( P*PV + N*NV )*m->PFM_CELL[ i ][ k ].dPPf ;

					xFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
					yFlux += faceFlux*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

	 			} else if( drift_diffusion.get_cell_typename( Cell_j->data_id ) != "PLASMA" ) {//Discontinue face

	 				if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ) {
	 					//do nothing
	 				} else {

 						/*--- Drift term ---*/
						U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;

 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;
						PV = var->U0[iSpecies][ i ] ;

 						faceFlux = ( vn*PV )*m->PFM_CELL[ i ][ k ].dArea*m->PFM_CELL[ i ][ k ].dPPf ;

 						xFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
						yFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					}
	 			}//For discontuity face
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					/*--- Drift term ---*/
					U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
					V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
					vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

					PV = var->U0[iSpecies][ i ] ;

					faceFlux = ( vn*PV )*m->PFM_CELL[ i ][ k ].dArea*m->PFM_CELL[ i ][ k ].dPPf ;

					xFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
					yFlux += faceFlux* m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

	 			}
	 		}


	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		xFlux = 0.0 ;
	 		yFlux = 0.0 ;
	 	}
	 	var->U1[ iSpecies ][ i ] = xFlux/Cell_i->volume ;
		var->U2[ iSpecies ][ i ] = yFlux/Cell_i->volume ;

	}//Cell Loop

	/*--- Update ghost cells ---*/
	var->U1[ iSpecies ] = var->U1[ iSpecies ] ;
	var->U2[ iSpecies ] = var->U2[ iSpecies ] ;
}
void CDriftDiffusion::CalculateAvgDDFlux_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double xFlux=0.0, yFlux=0.0, faceFlux=0.0, vn=0.0, U=0.0, V=0.0, Pe=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0 ;
	double P=0.0, N=0.0, PV=0.0, NV=0.0, IonFlux=0.0 ;
	double E = 0.0, Diff=0.0, Mobi=0.0, TempGradient=0.0, dL=0.0, dR=0.0, f1=0.0, f2=0.0 ;

	Cell *Cell_i, *Cell_j ;
	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {
		// cout<<"i: "<<i<<endl;
		Cell_i = drift_diffusion.get_cell(i) ;

		xFlux = 0.0 ;
		yFlux = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ;
				Cell_j = drift_diffusion.get_cell(j) ;

				if (  drift_diffusion.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					/*--- S-G Scheme ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ] ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					Diff = dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ;
					Pe = vn*m->PFM_CELL[ i ][ k ].dDist/Diff ;

					/*--- S-G ---*/
					if ( Pe < -ZERO ) {

						P = vn*(     - 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						N = vn*( 1.0 + 1.0/( exp(-Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;

					} else if ( Pe > ZERO ) {

						P = vn*( 1.0 + 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;
						N = vn*(     - 1.0/( exp( Pe)-1.0) )*m->PFM_CELL[ i ][ k ].dArea ;

					} else {

						Diff = -( dL*var->Diff[iSpecies][ i ] + dR*var->Diff[iSpecies][ j ] ) ;
						P = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;
						N =  Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;

					}
					//if( m->PFM_CELL[ i ][ k ].nf[ 1 ] != 0.0 and iSpecies == 0 ) cout<<Pe<<endl;
					PV = var->U0[iSpecies][ i ] ;
					NV = var->U0[iSpecies][ j ] ;

					faceFlux = ( P*PV + N*NV )*m->PFM_CELL[ i ][ k ].dPPf ;

					xFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
					yFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					// cout<<"faceFlux: "<<faceFlux<<endl;
					// cout<<"dPPF:     "<<m->PFM_CELL[ i ][ k ].dPPf<<endl;
					// cout<<"nf[0]: "<<m->PFM_CELL[ i ][ k ].nf[ 0 ]<<endl ;
					// cout<<"nf[1]: "<<m->PFM_CELL[ i ][ k ].nf[ 1 ]<<endl ;
					// cout<<"x0["<<k<<"]: "<< xFlux <<endl;
	 			// 	cout<<"y0["<<k<<"]: "<< yFlux <<endl;
					//if( iSpecies == 0 ) cout<<yFlux/Cell_i->volume<<endl;

	 			} else if( drift_diffusion.get_cell_typename( Cell_j->data_id ) != "PLASMA" ) {//Discontinue face

					/*--- Diffusion flux ---*/
					Diff = -var->Diff[iSpecies][ i ] ; 
					P    = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;
					PV   = var->U0[iSpecies][ i ] ;

					faceFlux= ( P*PV - 0.0 )*m->PFM_CELL[ i ][ k ].dPPf ;
						
					xFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
					yFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					//cout<<"x1["<<k<<"]: "<< xFlux <<endl;
	 				//cout<<"y1["<<k<<"]: "<< yFlux <<endl;

	 			}//For discontuity face
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					/*--- Diffusion flux ---*/
					Diff = -var->Diff[iSpecies][ i ] ; 
					P = -Diff/m->PFM_CELL[ i ][ k ].dDist*m->PFM_CELL[ i ][ k ].dArea ;
					PV = var->U0[iSpecies][ i ] ;
						
					faceFlux= ( P*PV - 0.0 )*m->PFM_CELL[ i ][ k ].dPPf ;
						
					xFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 0 ] ;
					yFlux += faceFlux * m->PFM_CELL[ i ][ k ].nf[ 1 ] ;
					//cout<<"x2["<<k<<"]: "<< xFlux <<endl;
	 				//cout<<"y2["<<k<<"]: "<< yFlux <<endl;
					//if( iSpecies == 0 ) cout<<yFlux/Cell_i->volume<<endl;

	 			}
	 		}


	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		xFlux = 0.0 ;
	 		yFlux = 0.0 ;
	 	}
	 	// cout<<"xFlux: "<<xFlux<<endl;
	 	// cout<<"yFlux: "<<yFlux<<endl<<endl;
	 	//if( iSpecies == 0 ) cout<<yFlux/Cell_i->volume<<endl;
	 	var->U1[ iSpecies ][ i ] = xFlux/Cell_i->volume ;
		var->U2[ iSpecies ][ i ] = yFlux/Cell_i->volume ;

	}//Cell Loop
	// exit(1) ;

	/*--- Update ghost cells ---*/
	var->U1[ iSpecies ] = var->U1[ iSpecies ] ;
	var->U2[ iSpecies ] = var->U2[ iSpecies ] ;
}
void CDriftDiffusion::CalculateDDConvection( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;
	double faceFlux=0.0, vn=0.0, U=0.0, V=0.0, ThermalVel=0.0, Te=0.0, SecondaryElectronEmission=0.0 ;
	double P=0.0, N=0.0, PV=0.0, NV=0.0, IonFlux=0.0 ;
	double E = 0.0, Diff=0.0, Mobi=0.0, dL=0.0, dR=0.0, f1=0.0, f2=0.0 ;

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell(i) ;

		var->DD_Convection[iSpecies][ i ] = 0.0 ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ;

				Cell_j = drift_diffusion.get_cell(j) ;

				if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) == "PLASMA" ){

					dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
					dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

					/*--- Upwind ---*/
					U = dL*config->Species[ iSpecies ].Charge * var->Ex[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ex[ j ] * var->Mobi[iSpecies][ j ] ;

					V = dL*config->Species[ iSpecies ].Charge * var->Ey[ i ] * var->Mobi[iSpecies][ i ] 
					  + dR*config->Species[ iSpecies ].Charge * var->Ey[ j ] * var->Mobi[iSpecies][ j ] ;

					vn = U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ;

					
					if ( vn > 0.0 ) {

						P = vn ;
						N = 0.0;

					} else if ( vn < 0.0 ) {

						P = 0.0;
						N = vn ;

					}
					PV = var->U0[iSpecies][ i ] ;
					NV = var->U0[iSpecies][ j ] ;

					faceFlux = ( P*PV + N*NV )*m->PFM_CELL[ i ][ k ].dArea ;

					var->DD_Convection[iSpecies][ i ] += faceFlux ;

	 			} else if( drift_diffusion.get_cell_typename( Cell_j->data_id ) != "PLASMA" ) {//Discontinue face

	 				if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ) {
	 					//do nothing
	 				} else {

						switch ( config->Species[ iSpecies ].Type ) {

							case 0:/*--- Electron ---*/

		 						/*--- Drift term ---*/
								U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
		 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
		 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

		 						/*--- Thermal flux term ---*/
		 						Te = var->T[ 0 ][ i ] ;
		 						if ( fixTe ) Te = 0.5 ;
		 						vn += 0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);

		 						/*--- Secondary electron emission ---*/
		 						SecondaryElectronEmission = 0.0 ;
		 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

		 							if (config->Species[ jSpecies ].Type == ION ){
		 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
		 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ;
		 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
		 							}
		 							
		 						} //if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) != DIELECTRIC ) SecondaryElectronEmission = 0.0 ;

								PV = var->U0[iSpecies][ i ] ;

		 						faceFlux = ( vn*PV - SecondaryElectronEmission)*m->PFM_CELL[ i ][ k ].dArea ;

								var->DD_Convection[iSpecies][ i ] += faceFlux ;

							break;

							case 1:/*--- Ion ---*/

								/*--- Drift term ---*/
								U = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
		 						V = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
		 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

								/*--- Thermal flux term ---*/
								//vn += 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;
								PV = var->U0[iSpecies][ i ] ;

		 						faceFlux=  ( vn*PV )*m->PFM_CELL[ i ][ k ].dArea ;

								var->DD_Convection[iSpecies][ i ] += faceFlux ;

							break;
							
							case 2:/*--- Neutral, Diffusion flux ---*/	

							break;

							default:
								if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
								exit(1);
				    		break;
						}//End switch
					}
	 			}//For discontuity face
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	 			if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){
	 				//do nothing
	 			}else{

					switch ( config->Species[ iSpecies ].Type ){

						case 0:/*--- Electron ---*/

	 						/*--- Drift term ---*/
							U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

	 						/*--- Thermal flux term ---*/
	 						Te = var->T[ 0 ][ i ] ;
	 						if ( fixTe ) Te = 0.5 ;
	 						vn += 0.25*sqrt( 8.0*var->Qe*Te / var->PI / (config->Species[ 0 ].Mass_Kg/var->Ref_Mass) )*(1.0-Reflec);

	 						/*--- Secondary electron emission ---*/
	 						SecondaryElectronEmission = 0.0 ;
	 						for ( int jSpecies = 1 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ){

	 							if (config->Species[ jSpecies ].Type == ION ){
	 								IonFlux = max( 0.0, config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ex[ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
	 												  + config->Species[jSpecies].Charge * var->Mobi[jSpecies][ i ]* var->Ey[ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] )*var->U0[jSpecies][ i ] ;
	 								SecondaryElectronEmission += config->SecondaryElectronEmissionCoeff*IonFlux ;
	 							}
	 							
	 						} //if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) != DIELECTRIC ) SecondaryElectronEmission = 0.0 ;


							PV = var->U0[iSpecies][ i ] ;

	 						faceFlux = ( vn*PV - SecondaryElectronEmission)*m->PFM_CELL[ i ][ k ].dArea ;

							var->DD_Convection[iSpecies][ i ] += faceFlux ;

						break;

						case 1:/*--- Ion ---*/

							/*--- Drift term ---*/
							U = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 						V = config->Species[ iSpecies ].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;
	 						vn = max( 0.0, U*m->PFM_CELL[ i ][ k ].nf[ 0 ] + V*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

							/*--- Thermal flux term ---*/
	 						//vn += 0.25*sqrt( 8.0*var->Qe*var->T[iSpecies][ i ] / var->PI / (config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass) ) ;
							PV = var->U0[iSpecies][ i ] ;

	 						faceFlux=  ( vn*PV )*m->PFM_CELL[ i ][ k ].dArea ;

							var->DD_Convection[iSpecies][ i ] += faceFlux ;

						break;

						
						case 2:/*--- Neutral ---*/	

						break;

						default:
							if( mpi_rank == 0 ) cout << "Continuity boundary condition error, Pls contact K.-L. Chen " << endl;
							exit(1);
			    		break;
					}//End switch
	 			}
	 		}


	 	/*--- Loop over SOLID cells ---*/
	 	} else {
			var->DD_Convection[iSpecies][ i ] = 0.0 ;
	 	}
	}//Cell Loop
}
void CDriftDiffusion::CalculateSurfaceCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int j=0 ;

	Cell *Cell_j, *Cell_i ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell( i ) ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = Cell_i->cell[k]->local_id ;

				Cell_j = drift_diffusion.get_cell(j) ; 

				if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) == "DIELECTRIC" ) {

					m->PFM_CELL[ i ][ k ].SurfaceCharge += var->Dt*var->Qe*config->Species[iSpecies].Charge
					*fabs( var->U1[ iSpecies ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 0 ] 
						+  var->U2[ iSpecies ][ i ]*m->PFM_CELL[ i ][ k ].nf[ 1 ] ) ;

				}//discontiuity face

	 		}//End bulk 

	 	/*--- Loop over DIELECTRIC cells ---*/
	 	} else if( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "DIELECTRIC" ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ){

				j = m->PFM_CELL[ i ][ k ].NeighborCellId ;

				if ( drift_diffusion.get_cell_typename( Cell_j->data_id ) == "PLASMA" ) {

					m->PFM_CELL[ i ][ k ].SurfaceCharge += var->Dt*var->Qe*config->Species[iSpecies].Charge
					*fabs( var->U1[ iSpecies ][ j ]*m->PFM_CELL[ i ][ k ].nf[ 0 ]*(-1.0) 
						+  var->U2[ iSpecies ][ j ]*m->PFM_CELL[ i ][ k ].nf[ 1 ]*(-1.0) ) ;

				}//discontiuity face

	 		}//End bulk 
	 	}
	}//Cell Loop
}

void CDriftDiffusion::CalculateCondCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{

	Cell *Cell_i, *Cell_j ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell( i ) ;

		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){
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
void CDriftDiffusion::Semi_Empirical_Temperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	double mass = 0.0, DriftVelocity = 0.0 ;

	double BackGroundMass = config->BGMass_Kg/var->Ref_Mass  ;
	double IonMass = config->Species[ iSpecies ].Mass_Kg/var->Ref_Mass ;
	double BackGroundTemp = config->BGTemperature_K ;
	double VTOT=0.0, U=0.0, V=0.0 ;
	Cell *Cell_i, Cell_j ;

	for( int i = 0 ; i < drift_diffusion.Mesh.cell_number ; i++ ) {

		Cell_i = drift_diffusion.get_cell(i) ;

		mass =  (BackGroundMass+IonMass)/(5.0*IonMass + 3.0*BackGroundMass )*BackGroundMass ;
		
		/*--- Loop over PLASMA cells ---*/
		if ( drift_diffusion.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

			U = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ex[ i ] ;
	 		V = config->Species[iSpecies].Charge * var->Mobi[iSpecies][ i ]* var->Ey[ i ] ;

			VTOT = ( U*U + V*V ) ;

	 		var->T[iSpecies][ i ] = BackGroundTemp + (IonMass+BackGroundMass)/( 5.0*IonMass + 3.0*BackGroundMass )*BackGroundMass*VTOT/var->Kb ;
	 		//var->T[iSpecies][ i ] /= var->Kb ;
	 		var->T[iSpecies][ i ] /= var->K2eV ;

	  	}
	}//Cell Loop
	MPI_Barrier(MPI_COMM_WORLD) ;
}
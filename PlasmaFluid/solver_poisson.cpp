
#include "solver_poisson.hpp"

using namespace std ;
CPoisson::CPoisson()
{
}
void CPoisson::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	Correction = config->Equation[ POISSON ].Correction ;
	if ( mpi_rank == 0 ){
		cout<<"Creat POISSON"<<endl;
		cout<<"Correction: "<<Correction<<endl;
	}
	its=0 ;
}
void CPoisson::Solve( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  )
{
	int nCell = m->local_cell_number ;
	Zero_Gradient( m, variable ) ;
	UpdateElectricalMap( config, variable ) ;

	/*--- Calculate the net charge density for poisson's source term ---*/
	if ( config->Equation[ POISSON ].Equation < 3 ) {
		Calculate_NetCharge( m, config, variable ) ;
	}


	/*--- Calculate the net charge density for poisson's source term ---*/
		switch ( config->Equation[ POISSON ].Equation ) {

			/*--- Orignal eqn. ---*/
			case 0:
				CalculatePermitt( m, config, variable ) ;
				break;
			/*--- Semi-Implicit accroding K. M. Lin et al. CPC 183 (2012) 1225–1236. ---*/			
			case 1:
				CalculateEffectivePermitt( m, config, variable ) ;
				break;
			/*--- Semi-Implicit poisson w/ ion prediction.  ---*/			
			case 2:
				CalculateEffectivePermittEleOnly( m, config, variable ) ;
				break;
			/*--- Semi-Implicit with diffusion term  ---*/			
			case 3:
				break;
			/*--- Error print out and exit progeam.  ---*/
			default:
				if( mpi_rank == 0 ) cout << "Poisson eqn. type error, Pls contact K.-L. Chen " << endl;
				exit(1);
			break;
		}//End switch

	/*--- 1st solve the equation w/o cross diffusion term ---*/
		if ( Correction == -1 ){
//cout<<"phi->A"<<endl;
			Bulid_A_B_0th( m, config, variable ) ;
//cout<<"phi->B"<<endl;
		} else {
			//cout<<"1st"<<endl;
			//Bulid_A_B_1st( m, config, variable ) ;
			Bulid_A_B_Orthogonal2( m, config, variable ) ;
		}
		//s.solve() ; 
		plasma.get_solution( 0, variable->Phi.data ) ;
		variable->Phi = variable->Phi ;


	/*--- 2nd solve the equation w/ cross diffusion term ---*/
		// for ( int k = 0 ; k < Correction ; k++ ) {

		// 	/*---Note: Calculate twice gradient for nuemann boundary condition. ---*/
		// 	Calculate_Gradient( m, config, variable ) ;
		// 	Calculate_Gradient( m, config, variable ) ;

		// 	variable->EField[0] = variable->EField[0] ; 
		// 	variable->EField[1] = variable->EField[1] ;
		//   	Bulid_A_B_2nd( m, config, variable ) ;
		//   	//s.solve(); 

		//   	its += s.its ;
		//   	for ( int i = 0 ; i < m->local_cell_number ; i++ ) {
		//   		variable->Phi[i] = s.x[i] ;
		//   		//cout<<variable->Phi[i]<<endl;
		//   	}
		//   	variable->Phi = variable->Phi ;
		// }
	//s.reuse_preconditioner( false ) ;

	/*--- Minus ---*/
		if ( Correction < 0 ){
			Calculate_Gradient_Neumann2( m, config, variable ) ;
			//Calculate_Gradient( m, config, variable ) ;
		} else {
			//Calculate_Gradient_GG( m, config, variable ) ;
			Calculate_Gradient_Neumann2( m, config, variable ) ;
		}
//cout<<"phi->G"<<endl;
		//CalculateReconstructionGradient( m, config, variable ) ;
		for ( int i = 0 ; i < nCell ; i++ ){
			variable->EField[ 0 ][ i ] = (-1.0)*variable->EField[ 0 ][ i ] ; 
			variable->EField[ 1 ][ i ] = (-1.0)*variable->EField[ 1 ][ i ] ; 
		} 
		variable->EField[0] = variable->EField[0] ; 
		variable->EField[1] = variable->EField[1] ;

		variable->CalReducedElectricField( m ) ;
	/*--- Displacement current density ---*/
		CalculateDispCurrentDensity( m, config, variable ) ;

}
void CPoisson::Bulid_A_B_Orthogonal( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int iFace=0, jFace=0, iCell=0, jCell=0, j=0, P=0, N=0, Pj=0, Nj=0, kk=0, col=0, jj ;
	// int nCell = m->local_cell_number ;
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0, ShapeFunction=0.0, CP=0.0, CN=0.0, BC_Value=0.0 ;

	// Cell *Cell_i, *Cell_j ;
	// //Face *Face_i ;

	// // 0 stand for poisson equation.
	// plasma.before_matrix_construction( 0 ) ;
	// plasma.before_source_term_construction( 0 ) ;

	// for( int i = 0 ; i < nCell ; i++ ) {

	// 	Cell_i  = plasma.get_cell( i ) ;

	// 	iFace 	 = Cell_i->face_number ;
	// 	iCell 	 = Cell_i->cell_number ;

	// 	/*--- Loop over electrode cells ---*/
	// 	if ( Cell_i->type == POWER or Cell_i->type == GROUND ) {

	// 		plasma.add_entry_in_matrix     ( 0, i, Cell_i->id, 1.0 ) ;
	// 		plasma.add_entry_in_source_term( 0, i, SineVoltage( Cell_i->type, config, var ) ) ;

	// 	} else {

	// 		for ( int k = 0 ; k < iCell ; k++ ) {

	// 			/*--- Orthogonal term ---*/
	// 			j = Cell_i->cell[ k ]->local_id ;

	// 			Cell_j = plasma.get_cell( j ) ;

	// 			jCell = Cell_j->cell_number ;
	// 			jFace = Cell_j->face_number ;

	// 			HarmonicMean = var->Eps[ i ]*var->Eps[ j ]
	// 						/( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
	// 						+  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;

	// 			Ad_dPN = HarmonicMean * ( m->PFM_CELL[ i ][ k ].dArea ) ;

	// 			plasma.add_entry_in_matrix( 0, i, Cell_i->id, -Ad_dPN ) ;
	// 			plasma.add_entry_in_matrix( 0, i, Cell_j->id,  Ad_dPN ) ;

	// 			//plasma.add_entry_in_source_term( i, SineVoltage( Cell_i->type, config, var ) ) ;
	// 			//Source += (-1.0)*(m->PFM_CELL[ i ][ k ].SurfaceCharge)*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf 
	// 			// 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) * m->PFM_CELL[ i ][ k ].dArea ;	
	// 		}//Loop over neighbor cells

	// 	/*--------------------------------------------------------------*/
	// 		for( int k = iCell ; k < iFace ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( Cell_i->face[ k ]->type == NEUMANN ){

	// 			} else if ( Cell_i->face[ k ]->type == GROUND or Cell_i->face[ k ]->type == POWER ){

	// 				electrode_voltage = SineVoltage( Cell_i->face[ k ]->type, config, var ) ;
	// 				//s.Add_Entries( row, row, -Ad_dPN ) ;
	// 				plasma.add_entry_in_source_term( 0, i, -Ad_dPN ) ;
	// 				//Source += -(electrode_voltage)*Ad_dPN ;
	// 				plasma.add_entry_in_source_term( i, -(electrode_voltage)*Ad_dPN ) ;

	// 			}else if ( Cell_i->face[ k ]->type == DIELECTRIC ) {

	// 				cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
	// 				exit(1) ;

	// 			}else{
	// 				cout<< Cell_i->face[ k ]->Typename<<endl;
	// 				cout<< Cell_i->face[ k ]->type<<endl;
	// 				cout<<"error-\" solver_poisson.cpp-3\""<<endl;
	// 			}
	// 		}//Loop over boundary face cells
	// 	/*--------------------------------------------------------------*/
	// 		//cout<<"Ad_dPN: "<<Ad_dPN<<endl;
	// 	}
	// 	//Source += -var->NetQ[ i ]*(m->cell[ i ].volume) ;
	// 	plasma.add_entry_in_source_term( 0, i, -var->NetQ[ i ]*Cell_i->volume ) ;
	// 	//s.push_source ( row, Source) ;
	// }//Cell Loop

	// plasma.finish_matrix_construction( 0 ) ;
	// plasma.finish_source_term_construction( 0 ) ;

}
void CPoisson::Bulid_A_B_Orthogonal2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int iFace=0, jFace=0, iCell=0, jCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double Ad_dPN=0.0, HarmonicMean=0.0, electrode_voltage=0.0, ShapeFunction=0.0 ;

	Cell *Cell_i, *Cell_j ;

	// 0 stand for poisson equation.
	plasma.before_matrix_construction( 0 ) ;
	plasma.before_source_term_construction( 0 ) ;

	for( int i = 0 ; i < nCell ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		iFace 	 = Cell_i->face_number ;
		iCell 	 = Cell_i->cell_number ;

		/*--- Loop over electrode cells ---*/
		if ( Cell_i->type == POWER or Cell_i->type == GROUND ) {

			plasma.add_entry_in_matrix     ( 0, i, Cell_i->id, 1.0 ) ;
			plasma.add_entry_in_source_term( 0, i, SineVoltage( Cell_i->type, config, var ) ) ;

		} else {

			for ( int k = 0 ; k < iCell ; k++ ) {
				/*--- Orthogonal term ---*/
				j = Cell_i->cell[ k ]->local_id ;

				Cell_j = plasma.get_cell( j ) ;

				jCell = Cell_j->cell_number ;
				jFace = Cell_j->face_number ;


				HarmonicMean = var->Eps[ i ]*var->Eps[ j ]
							/( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
							+  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;

				Ad_dPN = HarmonicMean * ( m->PFM_CELL[ i ][ k ].dArea )/var->Eps[ i ] ;

				plasma.add_entry_in_matrix( 0, i, Cell_i->id, -Ad_dPN ) ;
				plasma.add_entry_in_matrix( 0, i, Cell_j->id,  Ad_dPN ) ;

			}//Loop over neighbor cells

		/*--------------------------------------------------------------*/
			for( int k = iCell ; k < iFace ; k++ ) {

				Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea)/var->Eps[ i ] ;

				if ( Cell_i->face[ k ]->type == NEUMANN ){

				} else if ( Cell_i->face[ k ]->type == GROUND or Cell_i->face[ k ]->type == POWER ){

					electrode_voltage = SineVoltage( Cell_i->face[ k ]->type, config, var ) ;

					plasma.add_entry_in_matrix( 0, i,  Cell_i->id, -Ad_dPN ) ;

					plasma.add_entry_in_source_term( 0, i, -(electrode_voltage)*Ad_dPN ) ;

				}else if ( Cell_i->face[ k ]->type == DIELECTRIC ) {

					cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
					exit(1) ;

				}else{
					cout<< Cell_i->face[ k ]->Typename<<endl;
					cout<< Cell_i->face[ k ]->type<<endl;
					cout<<"error-\" solver_poisson.cpp-3\""<<endl;
				}

			}//Loop over boundary face cells
		/*--------------------------------------------------------------*/
		}
		plasma.add_entry_in_source_term( 0, i, -var->NetQ[ i ]*Cell_i->volume/var->Eps[ i ] ) ;
	}//Cell Loop
	plasma.finish_matrix_construction( 0 ) ;
	plasma.finish_source_term_construction( 0 ) ;
}
void CPoisson::Bulid_A_B_0th( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int iFace=0, jFace=0, iCell=0, jCell=0, j=0, P=0, N=0, Pj=0, Nj=0, kk=0, col=0, jj ;
	// int nCell = m->local_cell_number ;//+ m->ghost_cell_number_level_1 ; //Include overloap cell
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0, ShapeFunction=0.0, CP=0.0, CN=0.0, BC_Value=0.0 ;

	// s.ZeroEntriesMat() ;
	// Cell *Cell_i, *Cell_j, *Cell_jj, *Cell_ii ;
	// for( int i = 0 ; i < nCell ; i++ ) {

	// 	Cell_i = plasma.get_cell(i) ;

	// 	iFace 	 = Cell_i->face_number ;
	// 	iCell 	 = Cell_i->cell_number ;

	// 	/*--- Reset  ---*/
	// 	row 	 = m->cell[ i ].id ;//Global id.
	// 	P = row ;
	// 	Source 	 = 0.0 ;

	// 	/*--- Loop over electrode cells ---*/
	// 	if ( Cell_i->type == POWER or Cell_i->type == GROUND ){

	// 		s.Add_Entries( row, P, 1.0 ) ;
	// 		Source += (SineVoltage( Cell_i->type, config, var )) ;

	// 	} else {

	// 		for ( int k = 0 ; k < iCell ; k++ ) {
	// 			/*--- Orthogonal term ---*/
	// 			j 	  = m->PFM_CELL[ i ][ k ].NeighborCellId ;

	// 			jCell = Cell_j->cell_number ;
	// 			jFace = Cell_j->face_number ;

	// 			N = m->PFM_CELL[ i ][ k ].NeighborGlobalCellId ;

	// 			HarmonicMean = var->Eps[ i ]*var->Eps[ j ]
	// 						/( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
	// 						+  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;

	// 			Ad_dPN = HarmonicMean * ( m->PFM_CELL[ i ][ k ].dArea ) ;

	// 			s.Add_Entries( row, P, -Ad_dPN ) ; 
	// 			s.Add_Entries( row, N,  Ad_dPN ) ; 

	// 			Source += (-1.0)*(m->PFM_CELL[ i ][ k ].SurfaceCharge)*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf 
	// 			 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) * m->PFM_CELL[ i ][ k ].dArea ;				 		

	// 			/*--- Non-Orthogonal term 
	// 				For ∇P gradient term @ P point, dVar = V_j - V_i		
	// 			*/
	// 			for ( int kk = 0 ; kk < iCell ; kk++ ) {

	// 				col = m->PFM_CELL[ i ][ kk ].NeighborGlobalCellId ;
	// 				jj 	= m->PFM_CELL[ i ][ kk ].NeighborCellId ;
	// 				//V_i term - P point
	// 				ShapeFunction =  var->LSQ_Cx[ kk ][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 0 ] 
	// 				 				+var->LSQ_Cy[ kk ][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 1 ] ; //(∇ dot P'P)

	// 				s.Add_Entries( row, P, -Ad_dPN*ShapeFunction*(-1.0) ) ;

	// 				//V_j term - N point
	// 				if ( Cell_i->type == m->cell[ jj ].type ){

	// 					s.Add_Entries( row, col, -Ad_dPN*ShapeFunction*(1.0) ) ;

	// 				} else {

	// 					if( m->cell[ jj ].type == POWER or m->cell[ jj ].type == GROUND ){
	// 						//calculate gradient using value in the surface.
	// 						electrode_voltage = SineVoltage( m->cell[ jj ].type, config, var ) ;
	// 						Source += (-1.0)*-Ad_dPN*ShapeFunction*electrode_voltage ;

	// 					} else {
	// 						//plasma-dielectric interface.
	// 						BC_Value = m->PFM_CELL[ i ][ kk ].dNPf*m->PFM_CELL[ i ][ kk ].dPPf * m->PFM_CELL[ i ][ kk ].SurfaceCharge
	// 							  / ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ kk ].dNPf ) ;
	// 						Source += (-1.0)*-Ad_dPN*ShapeFunction*BC_Value ;
	// 						/*
	// 							BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
	// 							/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );
	// 						*/
	// 						CN = ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf ) / ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ kk ].dNPf );
	// 						CP = ( var->Eps[ i  ]*m->PFM_CELL[ i ][ kk ].dNPf ) / ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ kk ].dNPf );
	// 						s.Add_Entries( row, col, (-Ad_dPN)*ShapeFunction*CN ) ;
	// 						s.Add_Entries( row,   P, (-Ad_dPN)*ShapeFunction*CP ) ;

	// 					}
	// 				}//End P point, interior
	// 			}//End ∇P.

	// 			for ( int kk = iFace ; kk < iCell ; kk++ ) {

	// 				//col = m->PFM_CELL[ i ][ kk ].NeighborGlobalCellId ;
	// 				//jj 	= m->PFM_CELL[ i ][ kk ].NeighborCellId ;
	// 				//V_i term - P point
	// 				ShapeFunction =  var->LSQ_Cx[kk][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 0 ] 
	// 								+var->LSQ_Cy[kk][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 1 ] ;

	// 				s.Add_Entries( row, P, (-Ad_dPN)*ShapeFunction*(-1.0) ) ;

	// 				//V_j term - N point
	// 				if ( m->cell[ i ].face[ kk ]->type == POWER or m->cell[ i ].face[ kk ]->type == GROUND ){

	// 					//calculate gradient using value in the surface.
	// 					electrode_voltage = SineVoltage( m->cell[ i ].face[ kk ]->type, config, var ) ;
	// 					Source += (-1.0)*(-Ad_dPN)*ShapeFunction*electrode_voltage ;

	// 				} else if( m->cell[ i ].face[ kk ]->type == NEUMANN ) {

	// 					s.Add_Entries( row, P, -Ad_dPN*ShapeFunction*(1.0) ) ;

	// 				} else {
	// 					cout<<"ERROR, No Dielectric in domain boundary @ ∇P."<<endl;exit(1) ;
	// 				}
	// 			}//End P point boundary face.
		
	// 			/*--- Non-Orthogonal term 
	// 				For ∇N gradient term @ N point, dVar = V_j - V_i		
	// 				for ( int kk = 0 ; kk < iCell ; kk++ ) {
	// 			*/		
	// 			if ( Cell_j->type == POWER or Cell_j->type == GROUND ){
	// 				//No gradient
	// 			} else {
	// 				for ( int kk = 0 ; kk < jCell ; kk++ ) {

	// 					col = m->Cell[ j ][ kk ].NeighborGlobalCellId ;
	// 					jj  = m->Cell[ j ][ kk ].NeighborCellId ;

	// 					ShapeFunction =  var->LSQ_Cx[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[0] 
	// 							  		+var->LSQ_Cy[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[1] ;
	// 					//V_i term
	// 					s.Add_Entries( row, N, Ad_dPN*ShapeFunction*(-1.0) ) ;

	// 					//V_j term
	// 					if ( Cell_j->type == m->cell[ jj ].type ){

	// 						s.Add_Entries( row, col, Ad_dPN*ShapeFunction*(1.0) ) ;

	// 					} else {

	// 						if ( m->cell[ jj ].type == POWER or m->cell[ jj ].type == GROUND ){

	// 							electrode_voltage = SineVoltage( m->cell[ jj ].type, config, var ) ;
	// 							Source += (-1.0)*Ad_dPN*ShapeFunction*electrode_voltage ;

	// 						} else {

	// 							BC_Value = m->Cell[ j ][ kk ].dNPf*m->Cell[ j ][ kk ].dPPf * m->Cell[ j ][ kk ].SurfaceCharge
	// 						 	 	  / ( var->Eps[ jj ]*m->Cell[ j ][ kk ].dPPf + var->Eps[ j ]*m->Cell[ j ][ kk ].dNPf ) ;
	// 						 	Source += (-1.0)*Ad_dPN*ShapeFunction*BC_Value ;
	// 						 	/*
	// 								BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
	// 						 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );
	// 						 	*/
	// 						 	CN = ( var->Eps[ jj]*m->Cell[ j ][ kk ].dPPf ) / ( var->Eps[ jj ]*m->Cell[ j ][ kk ].dPPf + var->Eps[ j ]*m->Cell[ j ][ k ].dNPf );
	// 						 	CP = ( var->Eps[ j ]*m->Cell[ j ][ kk ].dNPf ) / ( var->Eps[ jj ]*m->Cell[ j ][ kk ].dPPf + var->Eps[ j ]*m->Cell[ j ][ k ].dNPf );
	// 						 	s.Add_Entries( row, col, Ad_dPN*ShapeFunction*CN ) ;
	// 						 	s.Add_Entries( row,   N, Ad_dPN*ShapeFunction*CP ) ;
	// 						}
	// 					}//NCell discontinue.
	// 				}//End N point, interior.

	// 				for ( int kk = jFace ; kk < jCell ; kk++ ) {

	// 					//col = m->Cell[ j ][ kk ].NeighborGlobalCellId ;
	// 					//jj  = m->Cell[ j ][ kk ].NeighborCellId ;

	// 					ShapeFunction =  var->LSQ_Cx[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[0] 
	// 							  		+var->LSQ_Cy[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[1] ;
	// 					//V_i term
	// 					s.Add_Entries( row, N, Ad_dPN*ShapeFunction*(-1.0) ) ;

	// 					if ( m->cell[ j ].face[ kk ]->type == GROUND or m->cell[ j ].face[ kk ]->type == POWER ){

	// 						electrode_voltage = SineVoltage( m->cell[ j ].face[ kk ]->type, config, var ) ;
	// 						Source += (-1.0)*Ad_dPN*ShapeFunction*electrode_voltage ;

	// 					} else if( m->cell[ j ].face[ kk ]->type == NEUMANN ) {

	// 						s.Add_Entries( row, N, Ad_dPN*ShapeFunction*(1.0) ) ;

	// 					} else {
	// 						cout<<"ERROR, No Dielectric in domain boundary @ ∇N"<<endl; exit(1) ;
	// 					}
	// 				}//End N-point boundary
	// 			}//End ∇N.	

	// 		}//Loop over neighbor cells

	// 	/*--------------------------------------------------------------*/
	// 		for( int k = iCell ; k < iFace ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( m->cell[ i ].face[ k ]->type == NEUMANN ){

	// 			} else if ( m->cell[ i ].face[ k ]->type == GROUND or m->cell[ i ].face[ k ]->type == POWER ){

	// 				electrode_voltage = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;
	// 				//if ( m->cell[ i ].face[ k ]->type == POWER ) var->Volt = electrode_voltage ;

	// 				s.Add_Entries( row, row, -Ad_dPN ) ;
	// 				Source += -(electrode_voltage)*Ad_dPN ;

	// 			}else if ( m->cell[ i ].face[ k ]->type == DIELECTRIC ) {

	// 				cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
	// 				exit(1) ;

	// 			}else{
	// 				cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 				cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 				cout<<"error-\" solver_poisson.cpp-3\""<<endl;
	// 			}
	// 		}//Loop over boundary face cells
	// 	/*--------------------------------------------------------------*/
	// 	}
	// 	Source += -var->NetQ[ i ]*(m->cell[ i ].volume) ;
	// 	s.push_source ( row, Source) ;
	// }//Cell Loop
}
void CPoisson::Bulid_A_B_1st( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int iFace=0, iCell=0, j=0 ;
	// int nCell = m->local_cell_number ;//+ m->ghost_cell_number_level_1 ; //Include overloap cell
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0 ;

	// for( int i = 0 ; i < nCell ; i++ ) {

	// 	iFace 	 = Cell_i->face_number ;
	// 	iCell 	 = Cell_i->cell_number ;

	// 	/*--- Reset  ---*/
	// 	row 	 = m->cell[ i ].id ;
	// 	col[ 0 ] = row ;
	// 	ncol 	 = 1 ;
	// 	Source 	 = 0.0 ;
	// 	for( int k = 0 ; k < 5 ; k++ ) C[ k ] = 0.0 ;

	// 	/*--- Loop over electrode cells ---*/
	// 	if ( Cell_i->type == POWER or Cell_i->type == GROUND){

	// 		C[ 0 ] = 1.0 ;
	// 		Source += (SineVoltage( Cell_i->type, config, var )) ;

	// 	} else {//PLASMA or DIELECTRIC

	// 		/*--- Loop over bulk faces ---*/
	// 		for( int k = 0 ; k < iCell ; k++ ){

	// 			j = Cell_i->cell[ k ].local_id ;

	// 			Cell_j = plasma.get_cell( j );

	// 			if ( Cell_j->type == PLASMA or Cell_j->type == DIELECTRIC ){

	// 				HarmonicMean = 	  var->Eps[ i ]*var->Eps[ j ]
	// 							   /( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
	// 							   +  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;
	// 							   //cout<<var->Eps[ i ]<<endl;
					
	// 				Ad_dPN = HarmonicMean * (m->PFM_CELL[ i ][ k ].dArea) ;

	// 				col[ ncol ] = m->PFM_CELL[ i ][ k ].NeighborGlobalCellId ;
	// 				C[ 0 ] 	  += -Ad_dPN ;
	// 				C[ncol]    =  Ad_dPN ;
	// 				ncol ++ ;
	// 				//rhs[j][i] += (-1.0)*val_e[j][i] * permitt_ij*dx[i+1] 
	// 				// /(  permitt_Rplus*dx[i] + permitt_ij*dx[i+1] ) * (domain->Ae[ j ][ i ]);
	// 				Source += (-1.0)*(m->PFM_CELL[ i ][ k ].SurfaceCharge)*var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) 
	// 				 		/ ( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) + var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) * (m->PFM_CELL[ i ][ k ].dArea) ;//;*m->cell[ i ].volume ;//
	// 				//Source += (-1.0)*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].SurfaceCharge/m->PFM_CELL[ i ][ k ].dPPf/( var->Eps[ j ]/m->PFM_CELL[ i ][ k ].dNPf + var->Eps[ i ]/m->PFM_CELL[ i ][ k ].dPPf )*m->cell[ i ].volume  ;
	// 			} else {/*--- For discontuity face ---*/

	// 				Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dPPf)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 				if ( Cell_j->type == GROUND or Cell_j->type == POWER ){

	// 					electrode_voltage = SineVoltage(  Cell_j->type, config, var ) ;
	// 					C[ 0 ] += -Ad_dPN ;
	// 					Source += -(electrode_voltage)*Ad_dPN ;

	// 				}else if ( m->cell[ i ].face[ k ]->type == NEUMANN ) {

	// 				}else{
	// 					cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 					cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 					cout<<"error-\" solver_poisson.cpp-1\""<<endl;
	// 				}
	// 			}//end discontuity face
	// 		}//end bulk face

	// 		/*--- Loop over boundary faces ---*/
	// 		for( int k = iCell ; k < iFace ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( m->cell[ i ].face[ k ]->type == NEUMANN ){

	// 			} else if ( m->cell[ i ].face[ k ]->type == GROUND or m->cell[ i ].face[ k ]->type == POWER ){

	// 				electrode_voltage = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;
	// 				if ( m->cell[ i ].face[ k ]->type == POWER ) var->Volt = electrode_voltage ;
	// 				C[ 0 ] += -Ad_dPN ;
	// 				Source += -(electrode_voltage)*Ad_dPN ;

	// 			}else if ( m->cell[ i ].face[ k ]->type == DIELECTRIC ) {

	// 				cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
	// 				exit(1) ;

	// 			}else{
	// 				cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 				cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 				cout<<"error-\" solver_poisson.cpp-3\""<<endl;
	// 			}
	// 		}
	// 	}
	// 	Source += -var->NetQ[ i ]*(m->cell[ i ].volume) ;

	// 	//cout<<"S: "<<Source<<endl;
	// 	//C[0] = 1.E- ;
	// 	// NormalizeCoeff[ i ] = C[0] ;
	// 	// for( int k = 0 ; k < 5 ; k++ ) C[ k ] /= NormalizeCoeff[ i ] ;
	// 	// Source /= NormalizeCoeff[ i ] ;

	// 	s.push_matrix( row, ncol, col, C ) ;
	// 	s.push_source ( row, Source) ;
	// 	//cout<<C[ 0 ] <<endl ; cout<<endl;
	// }//Cell Loop
}
void CPoisson::Bulid_A_B_2nd( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// /*--- ---*/

	// /*--- ---*/
	// int iFace=0, iCell=0, j=0 ;
	// int nCell = m->local_cell_number ;//+ m->ghost_cell_number_level_1 ; //Include overloap cell
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0, P=0.0, N=0.0 ;

	// for( int i = 0 ; i < nCell ; i++ ) {

	// 	iFace 	 = Cell_i->face_number ;
	// 	iCell 	 = Cell_i->cell_number ;

	// 	/*--- Reset  ---*/
	// 	Source 	 = 0.0 ;
	// 	row 	 = m->cell[ i ].id ;
	// 	/*--- Loop over electrode cells ---*/
	// 	if ( Cell_i->type == POWER or Cell_i->type == GROUND){

	// 		Source += (SineVoltage( Cell_i->type, config, var )) ;

	// 	} else {

	// 		/*--- Loop over bulk faces ---*/
	// 		for( int k = 0 ; k < iCell ; k++ ){

	// 			j = Cell_i->cell[ k ].local_id ;

	// 			if ( Cell_j->type == PLASMA or Cell_j->type == DIELECTRIC ){

	// 				HarmonicMean = 	  var->Eps[ i ]*var->Eps[ j ]
	// 							   /( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf)
	// 							   +  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;
	// 							   //cout<<var->Eps[ i ]<<endl;
					
	// 				Ad_dPN = HarmonicMean * (m->PFM_CELL[ i ][ k ].dArea) ;

	// 				P = -Ad_dPN ;
	// 				N =  Ad_dPN ;

	// 				GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 				GVarP[ 1 ] = var->EField[ 1 ][ i ] ;

	// 				GVarN[ 0 ] = var->EField[ 0 ][ j ] ;
	// 				GVarN[ 1 ] = var->EField[ 1 ][ j ] ;	

	// 				Source += (-1.0)*( DotProduct( GVarN, m->PFM_CELL[ i ][ k ].NNP )*N 
	// 								 + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP )*P ) ;

	// 				// Source += (-1.0)*m->PFM_CELL[ i ][ k ].SurfaceCharge* var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf 
	// 				// 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf )* m->PFM_CELL[ i ][ k ].dArea ;
	// 				Source += (-1.0)*(m->PFM_CELL[ i ][ k ].SurfaceCharge)*var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) 
	// 				 		/ ( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) + var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) * (m->PFM_CELL[ i ][ k ].dArea) ;

	// 			} else {/*--- For discontuity face ---*/

	// 				Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 				if ( Cell_j->type == DIELECTRIC ){

	// 					cout<<"POISSON: dummy if-"<<endl ;
	// 					exit(1) ;

	// 				} else if ( Cell_j->type == GROUND or Cell_j->type == POWER ){


	// 					electrode_voltage = SineVoltage( Cell_j->type, config, var ) ;

	// 					if ( Cell_j->type == POWER ) var->Volt = electrode_voltage ;

	// 					P += -Ad_dPN ;
	// 					Source += -(electrode_voltage)*Ad_dPN ;


	// 					GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 					GVarP[ 1 ] = var->EField[ 1 ][ i ] ;
	// 					Source += (-1.0)*( 0.0 + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)*P ) ;

	// 				}else if ( m->cell[ i ].face[ k ]->type == NEUMANN ) {

	// 				}else{
	// 					cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 					cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 					cout<<"error-\" solver_poisson.cpp-1 \""<<endl;
	// 				}


	// 			}//end discontuity face
	// 		}//end bulk face

	// 		/*--- Loop over boundary faces ---*/
	// 		for( int k = iCell ; k < iFace ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( m->cell[ i ].face[ k ]->type == NEUMANN ){

	// 			} else if ( m->cell[ i ].face[ k ]->type == GROUND or m->cell[ i ].face[ k ]->type == POWER ){

	// 				electrode_voltage = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;
	// 				if ( m->cell[ i ].face[ k ]->type == POWER ) var->Volt = electrode_voltage ;

	// 				P = -Ad_dPN ;
	// 				Source += -(electrode_voltage)*Ad_dPN ;

	// 				GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 				GVarP[ 1 ] = var->EField[ 1 ][ i ] ;
	// 				Source += (-1.0)*( 0.0 + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)*P ) ;

	// 			}else if ( m->cell[ i ].face[ k ]->type == DIELECTRIC ) {
	// 				cout<<"Boundary face can not be DIELECTRIC, pls check w/ K.-L. Chen"<<endl;
	// 				exit(1) ;

	// 			}else{
	// 				cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 				cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 				cout<<"error-\" solver_poisson.cpp -4\""<<endl;
	// 			}
	// 		}
	// 	}
		
	// 	Source += -var->NetQ[ i ]*(m->cell[ i ].volume) ;

	// 	// Source /= NormalizeCoeff[ i ] ;
	// 	s.push_source ( row, Source) ;
	// }//Cell Loop
}
void CPoisson::Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int nCell = m->local_cell_number ;
	// int iFace=0, iCell=0, j=0, NeighborCellIndex=0 ;
	// double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, L=0.0, R=0.0;

	// Cell *Cell_i ;

	// for ( int i = 0 ; i < nCell ; i++ ) {

	// 	Cell_i = plasma.get_cell( i ) ;

	// 	Gx = 0.0 ; Gy = 0.0 ;
	// 	iCell = Cell_i->cell_number ;  
	// 	iFace = Cell_i->face_number ; 
		
	// 	if ( Cell_i->type == POWER ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else if ( Cell_i->type == GROUND ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else {

	// 			/*--- Loop over neighbor "faces" ---*/
	// 		for ( int k = 0 ; k < iCell ; k++ ) {

	// 			j = Cell_i->cell[ k ].local_id ;

	// 			if ( Cell_i->type != Cell_j->type ) {//For discontinued face

	// 				if (  Cell_i->type == POWER or Cell_i->type == GROUND or  Cell_j->type == POWER or Cell_j->type == GROUND ){

	// 					BC_Value = SineVoltage( Cell_i->type, config, var ) ;
	// 					dVar = BC_Value - var->Phi[ i ] ;

	// 				} else {

	// 					BC_Value = m->PFM_CELL[ i ][ k ].dNPf*m->PFM_CELL[ i ][ k ].dPPf * m->PFM_CELL[ i ][ k ].SurfaceCharge
	// 					 		 / ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) ;

	// 					BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
	// 					 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );

									
	// 					dVar = BC_Value - var->Phi[ i ] ;

	// 				}

	// 			} else {

	// 				dVar = var->Phi[ j ] - var->Phi[ i ] ;

	// 			}
	// 			Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	// 	    	Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;
	// 		}

	// 		/*--- Loop over boundary faces ---*/
	// 		for ( int k = iCell ; k < iFace ; k++ ) {

	// 			if ( m->cell[ i ].face[ k ]->type == NEUMANN )	{

	// 				GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 				GVarP[ 1 ] = var->EField[ 1 ][ i ] ;
	// 				BC_Value = var->Phi[ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)  ;
	// 					//cout<<BC_Value - var->Phi[ i ] <<endl;
	// 			} else if( m->cell[ i ].face[ k ]->type != PLASMA ) {

	// 				BC_Value = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;
	// 				//cout<<BC_Value<<endl;

	// 			} else {

	// 				cout<<"ERR@Poisson 303"<<endl ;

	// 			}
				
	// 			dVar = BC_Value - var->Phi[ i ] ;
	// 			//cout<<var->Phi[ i ]<<endl;
	// 			Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	// 	   		Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

	// 		}
	// 		var->EField[ 0 ][ i ] = Gx ;
	// 		var->EField[ 1 ][ i ] = Gy ;

	// 	}//For Dielectric or Plasma
	// }//Loop over all cells
}
void CPoisson::Calculate_Gradient_Neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int nCell = m->local_cell_number ;
	// int iFace=0, iCell=0, j=0, NeighborCellIndex=0 ;
	// double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, L=0.0, R=0.0;
	// Cell *Cell_i, *Cell_j ;

	// for ( int i = 0 ; i < nCell ; i++ ) {

	// 	Cell_i = plasma.get_cell( i ) ;

	// 	Gx = 0.0 ; Gy = 0.0 ;
	// 	iCell = Cell_i->cell_number ;  
	// 	iFace = Cell_i->face_number ; 
		
	// 	if ( Cell_i->type == POWER or  Cell_i->type == GROUND ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else {//could be dielectric or plasma

	// 		/*--- Loop over neighbor "cells" ---*/
	// 		for ( int k = 0 ; k < iCell ; k++ ) {


	// 			j = Cell_i->cell[ k ].local_id ;
	// 			Cell_j = plasma.get_cell( j ) ;

	// 			if ( Cell_i->type != Cell_j->type ) {//For discontinued face

	// 				//If neighboring cell is "electrode", the boundary value should calculate by SineVoltage function. 
	// 				if ( Cell_j->type == POWER or Cell_j->type == GROUND ) {

	// 					BC_Value = SineVoltage( Cell_j->type, config, var ) ;
	// 					dVar = BC_Value - var->Phi[ i ] ;

	// 				//If neighboring cell is "dielectric", the boundary value should calculate by formula. 
	// 				} else {

	// 					BC_Value = m->PFM_CELL[ i ][ k ].dNPf*m->PFM_CELL[ i ][ k ].dPPf * m->PFM_CELL[ i ][ k ].SurfaceCharge
	// 					 		 / ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) ;

	// 					BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
	// 					 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );

									
	// 					dVar = BC_Value - var->Phi[ i ] ;
	// 				}

	// 			} else {

	// 				dVar = var->Phi[ j ] - var->Phi[ i ] ;

	// 			}
	// 			Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	// 	    	Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;
	// 		}
	// 		//cout<<"Loop over boundary faces"<<endl;


	// 		/*--- Loop over boundary faces ---*/
	// 		for ( int k = iCell ; k < iFace ; k++ ) {

	// 			if ( Cell_i->face[ k ]->type == NEUMANN )	{

	// 				BC_Value = var->Phi[ i ] ;

	// 			} else if( Cell_i->face[ k ]->type != PLASMA ) {

	// 				BC_Value = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;

	// 			} else {

	// 				cout<<"ERR@Poisson 303"<<endl ;

	// 			}
				
	// 			dVar = BC_Value - var->Phi[ i ] ;
	// 			//cout<<var->Phi[ i ]<<endl;
	// 			Gx = Gx + m->LSQ[ i ].Cx[ k ]*dVar ;
	// 	   	Gy = Gy + m->LSQ[ i ].Cy[ k ]*dVar ;

	// 		}
	// 		var->EField[ 0 ][ i ] = Gx ;
	// 		var->EField[ 1 ][ i ] = Gy ;

	// 	}//For Dielectric or Plasma
	// }//Loop over all cells
}
void CPoisson::Calculate_Gradient_Neumann2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int nCell = m->local_cell_number ;
	int iFace=0, iCell=0, j=0, NeighborCellIndex=0 ;
	double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, L=0.0, R=0.0;

	Cell *Cell_i, *Cell_j ;

	for ( int i = 0 ; i < nCell ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		Gx = 0.0 ; Gy = 0.0 ;
		iCell = Cell_i->cell_number ;  
		iFace = Cell_i->face_number ; 
		
		if (Cell_i->type == POWER or Cell_i->type == GROUND ){

			var->EField[ 0 ][ i ] = 0.0 ;
			var->EField[ 1 ][ i ] = 0.0 ;

		} else {//could be dielectric or plasma

			/*--- Loop over neighbor "cells" ---*/
			for ( int k = 0 ; k < iCell ; k++ ) {

				j = Cell_i->cell[ k ]->local_id ;

				Cell_j = plasma.get_cell( j ) ;

				if (Cell_i->type != Cell_j->type ) {//For discontinued face

					//If neighboring cell is "electrode", the boundary value should calculate by SineVoltage function. 
					if ( Cell_j->type == POWER or Cell_j->type == GROUND ) {

						BC_Value = SineVoltage( Cell_j->type, config, var ) ;
						dVar = BC_Value - var->Phi[ i ] ;

					//If neighboring cell is "dielectric", the boundary value should calculate by formula. 
					} else {

						BC_Value = m->PFM_CELL[ i ][ k ].dNPf*m->PFM_CELL[ i ][ k ].dPPf * m->PFM_CELL[ i ][ k ].SurfaceCharge
						 		 / ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) ;

						BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
						 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );

									
						dVar = BC_Value - var->Phi[ i ] ;
					}

				} else {

					dVar = var->Phi[ j ] - var->Phi[ i ] ;

				}
				Gx = Gx + var->LSQ_Cx[ k ][ i ]*dVar ;
		    Gy = Gy + var->LSQ_Cy[ k ][ i ]*dVar ;
			}
			//cout<<"Loop over boundary faces"<<endl;


			/*--- Loop over boundary faces ---*/
			for ( int k = iCell ; k < iFace ; k++ ) {

				if (Cell_i->face[ k ]->type == NEUMANN ) {

					BC_Value = var->Phi[ i ] ;

				} else if(Cell_i->face[ k ]->type != PLASMA ) {

					BC_Value = SineVoltage( Cell_i->face[ k ]->type, config, var ) ;

				} else {

					cout<<"ERR@Poisson 303"<<endl ;

				}
				
				dVar = BC_Value - var->Phi[ i ] ;
				//cout<<var->Phi[ i ]<<endl;
				Gx = Gx + var->LSQ_Cx[ k ][ i ]*dVar ;
		    Gy = Gy + var->LSQ_Cy[ k ][ i ]*dVar ;

			}
			var->EField[ 0 ][ i ] = Gx ;
			var->EField[ 1 ][ i ] = Gy ;

		}//For Dielectric or Plasma
	}//Loop over all cells
}
void CPoisson::Calculate_Gradient_GG( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int nCell = m->local_cell_number ;
	// int iFace=0, iCell=0, j=0, NeighborCellIndex=0 ;
	// double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, L=0.0, R=0.0;
	// double Vf = 0.0, dL=0.0, dR=0.0 ; ;
	// for ( int i = 0 ; i < nCell ; i++ ) {

	// 	Gx = 0.0 ; Gy = 0.0 ;
	// 	iCell = Cell_i->cell_number ;  
	// 	iFace = Cell_i->face_number ; 
		
	// 	if ( Cell_i->type == POWER or  Cell_i->type == GROUND ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else {

	// 		//cout<<i<<endl;

	// 			/*--- Loop over neighbor "faces" ---*/
	// 		for ( int k = 0 ; k < iCell ; k++ ) {

	// 			j = Cell_i->cell[ k ].local_id ;

	// 			dL = m->PFM_CELL[ i ][ k ].dNPf / m->PFM_CELL[ i ][ k ].dDist ;
	// 			dR = m->PFM_CELL[ i ][ k ].dPPf / m->PFM_CELL[ i ][ k ].dDist ;

	// 			if ( Cell_i->type != Cell_j->type ) {//For discontinued face

	// 				//cout<<"For discontinued face"<<endl;

	// 				if ( Cell_j->type == POWER or Cell_j->type == GROUND ){

	// 					Vf = SineVoltage( Cell_j->type, config, var ) ;

	// 				} else {

	// 					Vf = m->PFM_CELL[ i ][ k ].dNPf*m->PFM_CELL[ i ][ k ].dPPf * m->PFM_CELL[ i ][ k ].SurfaceCharge
	// 					 		 / ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) ;
	// 					Vf +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
	// 					 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );

	// 				}

	// 			} else {

	// 				Vf = dL*var->Phi[ i ] + dR*var->Phi[ j ] ; 

	// 				//dVar = var->Phi[ j ] - var->Phi[ i ] ;

	// 			}
	// 			Gx = Gx + Vf*m->PFM_CELL[ i ][ k ].Af[ 0 ] ;
	// 	    	Gy = Gy + Vf*m->PFM_CELL[ i ][ k ].Af[ 1 ] ;
	// 		}
	// 		//cout<<"Loop over boundary faces"<<endl;


	// 		/*--- Loop over boundary faces ---*/
	// 		for ( int k = iCell ; k < iFace ; k++ ) {

	// 			if ( m->cell[ i ].face[ k ]->type == NEUMANN )	{

	// 				Vf = var->Phi[ i ] ;

	// 			} else if( m->cell[ i ].face[ k ]->type != PLASMA ) {

	// 				Vf = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;

	// 			} else {

	// 				cout<<"ERR@Poisson 303"<<endl ;

	// 			}

	// 			Gx = Gx + Vf*m->PFM_CELL[ i ][ k ].Af[ 0 ] ;
	// 	    	Gy = Gy + Vf*m->PFM_CELL[ i ][ k ].Af[ 1 ] ;

	// 		}
	// 		var->EField[ 0 ][ i ] = Gx/m->cell[ i ].volume ;
	// 		var->EField[ 1 ][ i ] = Gy/m->cell[ i ].volume ;

	// 	}//For Dielectric or Plasma
	// }//Loop over all cells
}
void CPoisson::Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	int nCell = m->local_cell_number  ;
	for ( int i = 0 ; i < nCell ; i++ ) {	
		var->EField[ 0 ][ i ] = 0.0 ;
		var->EField[ 1 ][ i ] = 0.0 ;
	}//Loop over all cells
	var->EField[ 0 ]=var->EField[ 0 ];
	var->EField[ 1 ]=var->EField[ 1 ];
}
void CPoisson::Calculate_NetCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
	int nCell = m->local_cell_number ;

	Cell *Cell_i ;

	var->NetQ.zero() ;

	for( int i = 0 ; i < nCell ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;


		if ( Cell_i->type == PLASMA ) {

			for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
				if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {
					var->NetQ[ i ] += ( var->Qe*config->Species[ jSpecies ].Charge*var->U0[ jSpecies ][ i ] ) ;
				}
			}//end jSpecies
		}//End Plasma

	}
	var->NetQ = var->NetQ ;
}
void CPoisson::Calculate_NetCharge_2N_minus_N( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// /*
	// 	∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	// */
	// int nCell = m->local_cell_number ;
	// double tmp=0.0 ;
	// for( int i = 0 ; i < nCell ; i++ ) {

	// 	var->NetQ[ i ] = 0.0 ;

	// 	if ( Cell_i->type == PLASMA ){

	// 		for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
	// 			if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {
	// 				tmp = var->Qe*config->Species[ jSpecies ].Charge*( 2.0*var->U0[ jSpecies ][ i ] - var->PreU0[ jSpecies ][ i ] + var->Dt*var->DD_Convection[jSpecies][i]*m->cell[ i ].volume ) ;
	// 				//tmp = var->Qe*config->Species[ jSpecies ].Charge*var->U0[ jSpecies ][ i ] ;
	// 				var->NetQ[ i ] += tmp ;
	// 			}
	// 		}//end jSpecies
	// 		//cout<<var->NetQ[ i ]<<endl;
	// 	}//End Plasma

	// }
	// //exit(1);
	// var->NetQ = var->NetQ ;
}
void CPoisson::CalculateEffectivePermitt( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
	double eps=0.0 ;

	Cell *Cell_i ;

	for( int i = 0 ; i < m->local_cell_number ; i++ ) {
		
		eps = 0.0 ;

		Cell_i  = plasma.get_cell( i ) ;

		if( Cell_i->type == PLASMA ){

			for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
				if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {
						eps += fabs( config->Species[ jSpecies ].Charge ) * (var->Mobi[jSpecies][ i ]) * (var->U0[ jSpecies ][ i ] ) ;
				}
			}//End jSpecies

		}//End PLASMA
		var->Eps[ i ] = (var->Eps0[ i ] + var->Qe*var->Dt*eps) ;
	}
	var->Eps = var->Eps ;
}
void CPoisson::CalculateEffectivePermittEleOnly( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// /*
	// 	∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	// */
	// double eps=0.0 ;

	// for( int i = 0 ; i < m->local_cell_number ; i++ ) {

	// 	eps = 0.0 ;
	// 	if( m -> cell[i].type == PLASMA ){

	// 		for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
	// 			if ( config->Species[ jSpecies ].Type == ELECTRON ) {
	// 					eps += fabs( config->Species[ jSpecies ].Charge ) * (var->Mobi[jSpecies][ i ]) * (var->U0[ jSpecies ][ i ]) ;
	// 			}
	// 		}//End jSpecies

	// 	}//End PLASMA
	// 	var->Eps[ i ] = (var->Eps0[ i ] + var->Qe*var->Dt*eps) ;
	// 	//cout<<var->Eps[ i ]<<endl;
	// }
	// var->Eps = var->Eps ;
	// //exit(1);
}
void CPoisson::CalculatePermitt( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
	for( int i = 0 ; i < m->local_cell_number ; i++ ) {

		var->Eps[ i ] = var->Eps0[ i ] ;

	}
	var->Eps = var->Eps ;
}
void CPoisson::UpdateElectricalMap( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) 
{
	map< int, CElectrical>::iterator Iter;

	for ( Iter=config->ElectricalMap.begin() ; Iter!=config->ElectricalMap.end() ; ++Iter ) {
    	Iter->second.UpdateVoltage( var->PhysicalTime*var->Ref_t ) ;
    	if ( Iter->first == POWER ) var->Volt = Iter->second.BC_Voltage ;
		//if ( mpi_rank == MASTER_NODE ) cout<<Iter->second.BC_Voltage<<endl;
    	//cout<<Iter->first<<"\t"<<Iter->second.BC_Voltage <<endl;
	}
	//exit(1) ;
}
double CPoisson::SineVoltage( int FaceType, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{


	map< int, CElectrical>::iterator Iter;
	Iter = config->ElectricalMap.find( FaceType ) ;
	return Iter->second.BC_Voltage ;

	// map< int, CElectrical>::iterator Iter;
	// double time=0.0 ;
	// Iter = config->ElectricalMap.find( FaceType ) ;
	// double AC_Amplitude = 0.5*Iter->second.Voltage_p2p ;

	// if ( var->PhysicalTime > Iter->second.Period ) {

	// 	int n_period = (int)( var->PhysicalTime/Iter->second.Period ) ;
	// 	time = var->PhysicalTime - (double)( n_period*Iter->second.Period ) ;

	// } else {

	// 	time = var->PhysicalTime ;
	
	// }

	// double duration = 2.0 * (var->PI) * time / Iter->second.Period ;
	// double voltage = AC_Amplitude*sin(duration) ;
	// if( fabs(voltage) < 1.E-10 ) voltage = 0.0 ;
	// return voltage + Iter->second.BiasVoltage ;
}
void CPoisson::CalculateDispCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int nCell = m->local_cell_number ;
	for( int i = 0 ; i < nCell ; i++ ) {
		 var->DispJD[0][i] = var->Eps0[ i ]*( var->EField[ 0 ][ i ] - var->PreEField[ 0 ][ i ] )/var->Dt ;
		 var->DispJD[1][i] = var->Eps0[ i ]*( var->EField[ 1 ][ i ] - var->PreEField[ 1 ][ i ] )/var->Dt ;
		 var->DispJD[2][i] = var->Eps0[ i ]*( var->EField[ 2 ][ i ] - var->PreEField[ 2 ][ i ] )/var->Dt ;

		var->TotalJD[0][i] = var->DispJD[0][ i ] ;
		var->TotalJD[1][i] = var->DispJD[1][ i ] ;
		var->TotalJD[2][i] = var->DispJD[2][ i ] ;

	}//Cell Loop
}


#include "solver_poisson.hpp"
#include "petscsys.h"
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
			Bulid_A_B_0th( m, config, variable ) ;
		} else {
			Bulid_A_B_Orthogonal2( m, config, variable ) ;
		}
		plasma.get_solution( variable->Phi.data ) ;

		variable->Phi = variable->Phi ;

	/*--- Minus ---*/
		if ( Correction < 0 ){
			Calculate_Gradient_Neumann2( m, config, variable ) ;
			//Calculate_Gradient( m, config, variable ) ;
		} else {
			//Calculate_Gradient_GG( m, config, variable ) ;
			Calculate_Gradient_Neumann2( m, config, variable ) ;
		}

		//CalculateReconstructionGradient( m, config, variable ) ;
		for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
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
	// int jFace=0, jCell=0, j=0, P=0, N=0, Pj=0, Nj=0, kk=0, col=0, jj ;
	// 
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0, ShapeFunction=0.0, CP=0.0, CN=0.0, BC_Value=0.0 ;

	// Cell *Cell_i, *Cell_j ;
	// //Face *Face_i ;

	// // 0 stand for poisson equation.
	// plasma.before_matrix_construction( 0 ) ;
	// plasma.before_source_term_construction( 0 ) ;

	// for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	Cell_i  = plasma.get_cell( i ) ;


	// 	/*--- Loop over electrode cells ---*/
	// 	if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ) {

	// 		plasma.add_entry_in_matrix     ( 0, i, Cell_i->id, 1.0 ) ;
	// 		plasma.add_entry_in_source_term( i, SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ) ;

	// 	} else {

	// 		for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

	// 			/*--- Orthogonal term ---*/
	// 			j = Cell_i->cell[ k ]->data_id ;

	// 			Cell_j = plasma.get_cell( j ) ;

	// 			jCell = Cell_j->cell_number ;
	// 			jFace = Cell_j->face_number ;

	// 			HarmonicMean = var->Eps[ i ]*var->Eps[ j ]
	// 						/( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
	// 						+  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;

	// 			Ad_dPN = HarmonicMean * ( m->PFM_CELL[ i ][ k ].dArea ) ;

	// 			plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN ) ;
	// 			plasma.add_entry_in_matrix( i, Cell_j->id,  Ad_dPN ) ;

	// 			//plasma.add_entry_in_source_term( i, SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ) ;
	// 			//Source += (-1.0)*(m->PFM_CELL[ i ][ k ].SurfaceCharge)*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf 
	// 			// 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) * m->PFM_CELL[ i ][ k ].dArea ;	
	// 		}//Loop over neighbor cells

	// 	/*--------------------------------------------------------------*/
	// 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){

	// 			} else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "GROUND" or plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "POWER" ){

	// 				electrode_voltage = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id) , config, var ) ;
	// 				//s.Add_Entries( row, row, -Ad_dPN ) ;
	// 				plasma.add_entry_in_source_term( i, -Ad_dPN ) ;
	// 				//Source += -(electrode_voltage)*Ad_dPN ;
	// 				plasma.add_entry_in_source_term( i, -(electrode_voltage)*Ad_dPN ) ;

	// 			}else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "DIELECTRIC" ) {

	// 				cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
	// 				exit(1) ;

	// 			}else{
	// 				cout<< plasma.get_face_typename( Cell_i->face[ k ]->data_id) name<<endl;
	// 				cout<< plasma.get_face_typename( Cell_i->face[ k ]->data_id) <<endl;
	// 				cout<<"error-\" solver_poisson.cpp-3\""<<endl;
	// 			}
	// 		}//Loop over boundary face cells
	// 	/*--------------------------------------------------------------*/
	// 		//cout<<"Ad_dPN: "<<Ad_dPN<<endl;
	// 	}
	// 	//Source += -var->NetQ[ i ]*(m->cell[ i ].volume) ;
	// 	plasma.add_entry_in_source_term( i, -var->NetQ[ i ]*Cell_i->volume ) ;
	// 	//s.push_source ( row, Source) ;
	// }//Cell Loop

	// plasma.finish_matrix_construction( 0 ) ;
	// plasma.finish_source_term_construction( 0 ) ;

}
void CPoisson::Bulid_A_B_Orthogonal2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int jFace=0, jCell=0, j=0 ;
	
	double Ad_dPN=0.0, HarmonicMean=0.0, electrode_voltage=0.0, ShapeFunction=0.0 ;

	Cell *Cell_i, *Cell_j ;

	// 0 stand for poisson equation.
	plasma.before_matrix_construction() ;
	plasma.before_source_term_construction() ;
	

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		/*--- Loop over electrode cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ) {

			plasma.add_entry_in_matrix     ( i, Cell_i->id, 1.0 ) ;
			plasma.add_entry_in_source_term( i, SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ) ;

		} else {

			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				/*--- Orthogonal term ---*/
				j = Cell_i->cell[ k ]->data_id ;

				Cell_j = plasma.get_cell( j ) ;

				jCell = Cell_j->cell_number ;
				jFace = Cell_j->face_number ;


				HarmonicMean = var->Eps[ i ]*var->Eps[ j ]
							/( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
							+  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;

				Ad_dPN = HarmonicMean * ( m->PFM_CELL[ i ][ k ].dArea )/var->Eps[ i ] ;

				plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN ) ;
				plasma.add_entry_in_matrix( i, Cell_j->id,  Ad_dPN ) ;

			}//Loop over neighbor cells

		/*--------------------------------------------------------------*/
			for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea)/var->Eps[ i ] ;

				if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id) == "NEUMANN" ){

				} else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id) == "GROUND" or plasma.get_face_typename( Cell_i->face[ k ]->data_id) == "POWER" ){

					electrode_voltage = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id), config, var ) ;

					plasma.add_entry_in_matrix     ( i,  Cell_i->id, -Ad_dPN ) ;
					plasma.add_entry_in_source_term( i, -(electrode_voltage)*Ad_dPN ) ;

				}else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "DIELECTRIC" ) {

					cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
					exit(1) ;

				}else{
					cout<<"error-\" solver_poisson.cpp-3\""<<endl;
				}

			}//Loop over boundary face cells
		/*--------------------------------------------------------------*/
		}
		plasma.add_entry_in_source_term( i, -var->NetQ[ i ]*Cell_i->volume/var->Eps[ i ] ) ;
	//	cout<<-var->NetQ[ i ]*Cell_i->volume/var->Eps[ i ]<<endl;
	}//Cell Loop
//	if(mpi_rank ==1){
//		cout<<"plasma.Mesh.cell_number: "<<plasma.Mesh.cell_number<<endl;
//	}

//	cout<<"rank: "<<mpi_rank<<"  before"<<endl;
	plasma.finish_matrix_construction() ;
	plasma.finish_source_term_construction() ;
//	cout<<"rank: "<<mpi_rank<<"  after"<<endl;
//		MPI_Barrier(MPI_COMM_WORLD);
}
void CPoisson::Bulid_A_B_0th( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int jFace=0, jCell=0, j=0, P=0, N=0, Pj=0, Nj=0, kk=0, col=0, jj ;
	//+ m->ghost_cell_number_level_1 ; //Include overloap cell
	double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0, ShapeFunction=0.0, CP=0.0, CN=0.0, BC_Value=0.0 ;

	// 0 stand for poisson equation.
	plasma.before_matrix_construction() ;
	plasma.before_source_term_construction() ;

	Cell *Cell_i, *Cell_j, *Cell_jj, *Cell_ii ;
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;


		/*--- Loop over electrode cells ---*/
		if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ){

			plasma.add_entry_in_matrix     ( i, Cell_i->id, 1.0 ) ;
			plasma.add_entry_in_source_term( i, SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ) ;


		} else {

			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				/*--- Orthogonal term ---*/
				j = Cell_i->cell[ k ]->data_id ;
				Cell_j = plasma.get_cell( j ) ;

				jCell = Cell_j->cell_number ;
				jFace = Cell_j->face_number ;

				HarmonicMean = var->Eps[ i ]*var->Eps[ j ]
							/( var->Eps[ j ]*(m->PFM_CELL[ i ][ k ].dPPf) 
							+  var->Eps[ i ]*(m->PFM_CELL[ i ][ k ].dNPf) ) ;

				Ad_dPN = HarmonicMean * ( m->PFM_CELL[ i ][ k ].dArea ) ;

				plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN ) ;
				plasma.add_entry_in_matrix( i, Cell_j->id,  Ad_dPN ) ;

				Source += (-1.0)*(m->PFM_CELL[ i ][ k ].SurfaceCharge)*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf 
				 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) * m->PFM_CELL[ i ][ k ].dArea ;				 		

				/*--- Non-Orthogonal term 
					For ∇P gradient term @ P point, dVar = V_jj - V_i		
				*/
				for ( int kk = 0 ; kk < Cell_i->cell_number ; kk++ ) {

					jj = Cell_i->cell[ kk ]->data_id ;
					Cell_jj = plasma.get_cell( jj ) ;

					//V_i term - P point
					ShapeFunction = var->LSQ_Cx[ kk ][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 0 ] 
					 							+ var->LSQ_Cy[ kk ][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 1 ] ; //(∇ dot P'P)

					//s.Add_Entries( row, P, -Ad_dPN*ShapeFunction*(-1.0) ) ;
					 plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN*ShapeFunction*(-1.0) ) ;

					//V_jj term
					if ( plasma.get_cell_typename( Cell_i->data_id ) == plasma.get_cell_typename( Cell_jj->data_id ) ){

						//s.Add_Entries( row, col, -Ad_dPN*ShapeFunction*(1.0) ) ;
						plasma.add_entry_in_matrix( i, Cell_jj->id,  -Ad_dPN*ShapeFunction*(1.0) ) ;

					} else {

						if( plasma.get_cell_typename( Cell_jj->data_id ) == "POWER" or plasma.get_cell_typename( Cell_jj->data_id ) == "GROUND" ){

							//calculate gradient using value in the surface.
							electrode_voltage = SineVoltage( plasma.get_cell_typename( Cell_jj->data_id ), config, var ) ;
							//Source += (-1.0)*-Ad_dPN*ShapeFunction*electrode_voltage ;
							plasma.add_entry_in_source_term( i, (-1.0)*-Ad_dPN*ShapeFunction*electrode_voltage ) ;

						} else {

							//plasma-dielectric interface.
							BC_Value = m->PFM_CELL[ i ][ kk ].dNPf*m->PFM_CELL[ i ][ kk ].dPPf * m->PFM_CELL[ i ][ kk ].SurfaceCharge
								  / ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ kk ].dNPf ) ;

							//Source += (-1.0)*-Ad_dPN*ShapeFunction*BC_Value ;
							plasma.add_entry_in_source_term( i, (-1.0)*-Ad_dPN*ShapeFunction*BC_Value ) ;

							/*
								BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
								/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );
							*/

							CN = ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf ) / ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ kk ].dNPf );
							CP = ( var->Eps[ i  ]*m->PFM_CELL[ i ][ kk ].dNPf ) / ( var->Eps[ jj ]*m->PFM_CELL[ i ][ kk ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ kk ].dNPf );
							//s.Add_Entries( row,   P, (-Ad_dPN)*ShapeFunction*CP ) ;
							//s.Add_Entries( row, col, (-Ad_dPN)*ShapeFunction*CN ) ;
							plasma.add_entry_in_matrix( i, Cell_i->id, (-Ad_dPN)*ShapeFunction*CP ) ;
							plasma.add_entry_in_matrix( i, Cell_j->id, (-Ad_dPN)*ShapeFunction*CN ) ;

						}
					}//End P point, interior
				}//End ∇P.

				/*--- Boundary face @ P-point ---*/
				for ( int kk = Cell_i->face_number ; kk < Cell_i->cell_number ; kk++ ) {

					//V_i term - P point
					ShapeFunction = var->LSQ_Cx[kk][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 0 ] 
												+ var->LSQ_Cy[kk][ i ]*m->PFM_CELL[ i ][ k ].PPP[ 1 ] ;

					//s.Add_Entries( row, P, (-Ad_dPN)*ShapeFunction*(-1.0) ) ;
					plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN*ShapeFunction*(-1.0) ) ;

					//V_j term - N point
					if ( plasma.get_face_typename( Cell_i->face[ kk ]->data_id) == "POWER" or plasma.get_face_typename( Cell_i->face[ kk ]->data_id) == "GROUND" ){

						//calculate gradient using value in the surface.
						electrode_voltage = SineVoltage( plasma.get_face_typename( Cell_i->face[ kk ]->data_id), config, var ) ;
						//Source += (-1.0)*(-Ad_dPN)*ShapeFunction*electrode_voltage ;
						plasma.add_entry_in_source_term( i,  (-1.0)*(-Ad_dPN)*ShapeFunction*electrode_voltage ) ;

					} else if( plasma.get_face_typename( Cell_i->face[ kk ]->data_id) == "NEUMANN" ) {

						//s.Add_Entries( row, P, -Ad_dPN*ShapeFunction*(1.0) ) ;
						plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN*ShapeFunction*( 1.0) ) ;

					} else {
						cout<<"ERROR, No Dielectric in domain boundary @ ∇P."<<endl ; exit(1) ;
					}
				}//End P point boundary face.


				/*--- Non-Orthogonal term 
					For ∇N gradient term @ N point, dVar = V_jj - V_i		
				*/		
				if ( plasma.get_cell_typename( Cell_j->data_id ) == "POWER" or plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" ){
					//No gradient inside electrode
				} else {

					for ( int kk = 0 ; kk < jCell ; kk++ ) {

						jj = Cell_j->cell[ kk ]->data_id ;
						Cell_jj = plasma.get_cell( jj ) ;

						ShapeFunction = var->LSQ_Cx[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[0] 
								  				+ var->LSQ_Cy[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[1] ;
						//V_i term
						//s.Add_Entries( row, N, Ad_dPN*ShapeFunction*(-1.0) ) ;
						plasma.add_entry_in_matrix( i, Cell_j->id, Ad_dPN*ShapeFunction*(-1.0) ) ;

						//V_jj term
						if ( plasma.get_cell_typename( Cell_j->data_id ) == plasma.get_cell_typename( Cell_jj->data_id ) ){

							plasma.add_entry_in_matrix( i, Cell_jj->id, Ad_dPN*ShapeFunction*( 1.0) ) ;
							//s.Add_Entries( row, col, Ad_dPN*ShapeFunction*(1.0) ) ;

						} else {

							if ( plasma.get_cell_typename( Cell_jj->data_id ) == "POWER" or plasma.get_cell_typename( Cell_jj->data_id ) == "GROUND" ){

								electrode_voltage = SineVoltage( plasma.get_cell_typename( Cell_jj->data_id ), config, var ) ;

								//Source += (-1.0)*Ad_dPN*ShapeFunction*electrode_voltage ;
								plasma.add_entry_in_source_term( i, (-1.0)*Ad_dPN*ShapeFunction*electrode_voltage ) ;


							} else {

								BC_Value = m->PFM_CELL[ j ][ kk ].dNPf*m->PFM_CELL[ j ][ kk ].dPPf * m->PFM_CELL[ j ][ kk ].SurfaceCharge
							 	 	  / ( var->Eps[ jj ]*m->PFM_CELL[ j ][ kk ].dPPf + var->Eps[ j ]*m->PFM_CELL[ j ][ kk ].dNPf ) ;
							 	Source += (-1.0)*Ad_dPN*ShapeFunction*BC_Value ;
							 	/*
									BC_Value +=( var->Phi[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Phi[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
							 		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );
							 	*/
							 	CN = ( var->Eps[ jj]*m->PFM_CELL[ j ][ kk ].dPPf ) / ( var->Eps[ jj ]*m->PFM_CELL[ j ][ kk ].dPPf + var->Eps[ j ]*m->PFM_CELL[ j ][ k ].dNPf );
							 	CP = ( var->Eps[ j ]*m->PFM_CELL[ j ][ kk ].dNPf ) / ( var->Eps[ jj ]*m->PFM_CELL[ j ][ kk ].dPPf + var->Eps[ j ]*m->PFM_CELL[ j ][ k ].dNPf );
							 	//s.Add_Entries( row,   N, Ad_dPN*ShapeFunction*CP ) ;
							 	//s.Add_Entries( row, col, Ad_dPN*ShapeFunction*CN ) ;
							 	plasma.add_entry_in_matrix( i, Cell_j ->id, (Ad_dPN)*ShapeFunction*CP ) ;
								plasma.add_entry_in_matrix( i, Cell_jj->id, (Ad_dPN)*ShapeFunction*CN ) ;
							}

						}//plasma.Mesh.cell_number discontinue.
					}//End N point, interior.

					for ( int kk = jFace ; kk < jCell ; kk++ ) {

						//col = m->Cell[ j ][ kk ].NeighborGlobalCellId ;
						jj = Cell_j->face[ kk ]->data_id ;

						ShapeFunction =  var->LSQ_Cx[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[0] 
								  					+var->LSQ_Cy[ kk ][ j ]*m->PFM_CELL[ i ][ k ].NNP[1] ;
						//V_i term
						plasma.add_entry_in_matrix( i, Cell_j->id, Ad_dPN*ShapeFunction*(-1.0) ) ;

						//V_jj term
						if ( plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" or plasma.get_face_typename( Cell_j->face[ k ]->data_id) == "POWER" ){

							electrode_voltage = SineVoltage( plasma.get_face_typename( Cell_j->face[ k ]->data_id), config, var ) ;
							//Source += (-1.0)*Ad_dPN*ShapeFunction*electrode_voltage ;
							plasma.add_entry_in_source_term( i,  (-1.0)*Ad_dPN*ShapeFunction*electrode_voltage ) ;

						} else if( plasma.get_face_typename( Cell_j->face[ k ]->data_id) == "NEUMANN" ) {

							//s.Add_Entries( row, N, Ad_dPN*ShapeFunction*(1.0) ) ;
							plasma.add_entry_in_matrix( i, Cell_j->id, Ad_dPN*ShapeFunction*(1.0) ) ;

						} else {
							cout<<"ERROR, No Dielectric in domain boundary @ ∇N"<<endl; exit(1) ;
						}
					}//End N-point boundary
				}//End ∇N.	
			}//Loop over neighbor cells

		/*--------------------------------------------------------------*/
			for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea)/var->Eps[ i ] ;

				if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ){

				} else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "GROUND" or plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "POWER" ){

					electrode_voltage = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id) , config, var ) ;

					plasma.add_entry_in_matrix( i,  Cell_i->id, -Ad_dPN ) ;

					plasma.add_entry_in_source_term( i, -(electrode_voltage)*Ad_dPN ) ;

				}else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "DIELECTRIC" ) {

					cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
					exit(1) ;

				}else{

					cout<<"error-\" solver_poisson.cpp-3\""<<endl;
				}

			}//Loop over boundary face cells
		/*--------------------------------------------------------------*/
		}
		plasma.add_entry_in_source_term( i, -var->NetQ[ i ]*Cell_i->volume/var->Eps[ i ] ) ;
	}//Cell Loop
}
void CPoisson::Bulid_A_B_1st( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// int j=0 ;
	// //+ m->ghost_cell_number_level_1 ; //Include overloap cell
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0 ;

	// for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {


	// 	/*--- Reset  ---*/
	// 	row 	 = m->cell[ i ].id ;
	// 	col[ 0 ] = row ;
	// 	ncol 	 = 1 ;
	// 	Source 	 = 0.0 ;
	// 	for( int k = 0 ; k < 5 ; k++ ) C[ k ] = 0.0 ;

	// 	/*--- Loop over electrode cells ---*/
	// 	if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND"){

	// 		C[ 0 ] = 1.0 ;
	// 		Source += (SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var )) ;

	// 	} else {//PLASMA or DIELECTRIC

	// 		/*--- Loop over bulk faces ---*/
	// 		for( int k = 0 ; k < Cell_i->cell_number ; k++ ){

	// 			j = Cell_i->cell[ k ].data_id ;

	// 			Cell_j = plasma.get_cell( j );

	// 			if ( plasma.get_cell_typename( Cell_j->data_id ) == "PLASMA" or plasma.get_cell_typename( Cell_j->data_id ) == "DIELECTRIC" ){

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

	// 				if ( plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" or plasma.get_cell_typename( Cell_j->data_id ) == "POWER" ){

	// 					electrode_voltage = SineVoltage(  plasma.get_cell_typename( Cell_j->data_id ), config, var ) ;
	// 					C[ 0 ] += -Ad_dPN ;
	// 					Source += -(electrode_voltage)*Ad_dPN ;

	// 				}else if ( m->cell[ i ].face[ k ]->type == "NEUMANN" ) {

	// 				}else{
	// 					cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 					cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 					cout<<"error-\" solver_poisson.cpp-1\""<<endl;
	// 				}
	// 			}//end discontuity face
	// 		}//end bulk face

	// 		/*--- Loop over boundary faces ---*/
	// 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( m->cell[ i ].face[ k ]->type == "NEUMANN" ){

	// 			} else if ( m->cell[ i ].face[ k ]->type == "GROUND" or m->cell[ i ].face[ k ]->type == "POWER" ){

	// 				electrode_voltage = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;
	// 				if ( m->cell[ i ].face[ k ]->type == "POWER" ) var->Volt = electrode_voltage ;
	// 				C[ 0 ] += -Ad_dPN ;
	// 				Source += -(electrode_voltage)*Ad_dPN ;

	// 			}else if ( m->cell[ i ].face[ k ]->type == "DIELECTRIC" ) {

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
	// int j=0 ;
	// //+ m->ghost_cell_number_level_1 ; //Include overloap cell
	// double Ad_dPN=0.0, HarmonicMean=0.0, Source=0.0, electrode_voltage=0.0, P=0.0, N=0.0 ;

	// for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	/*--- Reset  ---*/
	// 	Source 	 = 0.0 ;
	// 	row 	 = m->cell[ i ].id ;
	// 	/*--- Loop over electrode cells ---*/
	// 	if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND"){

	// 		Source += (SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var )) ;

	// 	} else {

	// 		/*--- Loop over bulk faces ---*/
	// 		for( int k = 0 ; k < Cell_i->cell_number ; k++ ){

	// 			j = Cell_i->cell[ k ].data_id ;

	// 			if ( plasma.get_cell_typename( Cell_j->data_id ) == "PLASMA" or plasma.get_cell_typename( Cell_j->data_id ) == "DIELECTRIC" ){

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

	// 				if ( plasma.get_cell_typename( Cell_j->data_id ) == "DIELECTRIC" ){

	// 					cout<<"POISSON: dummy if-"<<endl ;
	// 					exit(1) ;

	// 				} else if ( plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" or plasma.get_cell_typename( Cell_j->data_id ) == "POWER" ){


	// 					electrode_voltage = SineVoltage( plasma.get_cell_typename( Cell_j->data_id ), config, var ) ;

	// 					if ( plasma.get_cell_typename( Cell_j->data_id ) == "POWER" ) var->Volt = electrode_voltage ;

	// 					P += -Ad_dPN ;
	// 					Source += -(electrode_voltage)*Ad_dPN ;


	// 					GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 					GVarP[ 1 ] = var->EField[ 1 ][ i ] ;
	// 					Source += (-1.0)*( 0.0 + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)*P ) ;

	// 				}else if ( m->cell[ i ].face[ k ]->type == "NEUMANN" ) {

	// 				}else{
	// 					cout<< m->cell[ i ].face[ k ]->Typename<<endl;
	// 					cout<< m->cell[ i ].face[ k ]->type<<endl;
	// 					cout<<"error-\" solver_poisson.cpp-1 \""<<endl;
	// 				}


	// 			}//end discontuity face
	// 		}//end bulk face

	// 		/*--- Loop over boundary faces ---*/
	// 		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 			Ad_dPN = var->Eps[ i ]/(m->PFM_CELL[ i ][ k ].dDist)*(m->PFM_CELL[ i ][ k ].dArea) ;

	// 			if ( m->cell[ i ].face[ k ]->type == "NEUMANN" ){

	// 			} else if ( m->cell[ i ].face[ k ]->type == "GROUND" or m->cell[ i ].face[ k ]->type == "POWER" ){

	// 				electrode_voltage = SineVoltage( m->cell[ i ].face[ k ]->type, config, var ) ;
	// 				if ( m->cell[ i ].face[ k ]->type == "POWER" ) var->Volt = electrode_voltage ;

	// 				P = -Ad_dPN ;
	// 				Source += -(electrode_voltage)*Ad_dPN ;

	// 				GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 				GVarP[ 1 ] = var->EField[ 1 ][ i ] ;
	// 				Source += (-1.0)*( 0.0 + DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)*P ) ;

	// 			}else if ( m->cell[ i ].face[ k ]->type == "DIELECTRIC" ) {
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
	// 
	// int j=0, NeighborCellIndex=0 ;
	// double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, L=0.0, R=0.0;

	// Cell *Cell_i ;

	// for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	Cell_i = plasma.get_cell( i ) ;

	// 	Gx = 0.0 ; Gy = 0.0 ;
		
	// 	if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else if ( plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else {

	// 			/*--- Loop over neighbor "faces" ---*/
	// 		for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

	// 			j = Cell_i->cell[ k ].data_id ;

	// 			if ( plasma.get_cell_typename( Cell_i->data_id ) != plasma.get_cell_typename( Cell_j->data_id ) ) {//For discontinued face

	// 				if (  plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" or  plasma.get_cell_typename( Cell_j->data_id ) == "POWER" or plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" ){

	// 					BC_Value = SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ;
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
	// 		for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 			if ( m->cell[ i ].face[ k ]->type == "NEUMANN" )	{

	// 				GVarP[ 0 ] = var->EField[ 0 ][ i ] ;
	// 				GVarP[ 1 ] = var->EField[ 1 ][ i ] ;
	// 				BC_Value = var->Phi[ i ] ;//+ DotProduct( GVarP, m->PFM_CELL[ i ][ k ].PPP)  ;
	// 					//cout<<BC_Value - var->Phi[ i ] <<endl;
	// 			} else if( m->cell[ i ].face[ k ]->type != "PLASMA" ) {

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
	// 
	// int j=0, NeighborCellIndex=0 ;
	// double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, L=0.0, R=0.0;
	// Cell *Cell_i, *Cell_j ;

	// for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	Cell_i = plasma.get_cell( i ) ;

	// 	Gx = 0.0 ; Gy = 0.0 ;
		
	// 	if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or  plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ){

	// 		var->EField[ 0 ][ i ] = 0.0 ;
	// 		var->EField[ 1 ][ i ] = 0.0 ;

	// 	} else {//could be dielectric or plasma

	// 		/*--- Loop over neighbor "cells" ---*/
	// 		for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {


	// 			j = Cell_i->cell[ k ].data_id ;
	// 			Cell_j = plasma.get_cell( j ) ;

	// 			if ( plasma.get_cell_typename( Cell_i->data_id ) != plasma.get_cell_typename( Cell_j->data_id ) ) {//For discontinued face

	// 				//If neighboring cell is "electrode", the boundary value should calculate by SineVoltage function. 
	// 				if ( plasma.get_cell_typename( Cell_j->data_id ) == "POWER" or plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" ) {

	// 					BC_Value = SineVoltage( plasma.get_cell_typename( Cell_j->data_id ), config, var ) ;
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
	// 		for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 			if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" )	{

	// 				BC_Value = var->Phi[ i ] ;

	// 			} else if( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  != "PLASMA" ) {

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
	
	int j=0, NeighborCellIndex=0 ;
	double dVar=0.0, Gx=0.0, Gy=0.0, Gz=0.0, BC_Value=0.0, L=0.0, R=0.0;

	Cell *Cell_i, *Cell_j ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i = plasma.get_cell(i) ;

		Gx = Gy = Gz = 0.0 ;

		if ( plasma.get_cell_typename( Cell_i->data_id ) == "POWER" or plasma.get_cell_typename( Cell_i->data_id ) == "GROUND" ){

			var->EField[ 0 ][ i ] = 0.0 ;
			var->EField[ 1 ][ i ] = 0.0 ;
			if( plasma.Mesh.ndim == 3 )  var->EField[ 2 ][ i ] = 0.0 ;

		} else {//could be dielectric or plasma

			/*--- Loop over neighbor "cells" ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				j = Cell_i->cell[ k ]->data_id ;
				Cell_j = plasma.get_cell( j ) ;

				if (plasma.get_cell_typename( Cell_i->data_id ) != plasma.get_cell_typename( Cell_j->data_id ) ) {//For discontinued face

					//If neighboring cell is "electrode", the boundary value should calculate by SineVoltage function. 
					if ( plasma.get_cell_typename( Cell_j->data_id ) == "POWER" or plasma.get_cell_typename( Cell_j->data_id ) == "GROUND" ) {

						BC_Value = SineVoltage( plasma.get_cell_typename( Cell_j->data_id ), config, var ) ;
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
			for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				if (plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "NEUMANN" ) {

					BC_Value = var->Phi[ i ] ;

				} else if(plasma.get_face_typename( Cell_i->face[ k ]->data_id)  != "PLASMA" ) {

					BC_Value = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id) , config, var ) ;

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
			if( plasma.Mesh.ndim == 3 ) var->EField[ 2 ][ i ] = Gx ;

		}//For Dielectric or Plasma

		//cout<<"Gy: "<<Gy<<endl;
	}//Loop over all cells
	//cout<<endl;
	//exit(1);
}
void CPoisson::Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var )
{
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {	
		var->EField[ 0 ][ i ] = 0.0 ;
		var->EField[ 1 ][ i ] = 0.0 ;
		if( plasma.Mesh.ndim == 3 ) var->EField[ 2 ][ i ] = 0.0 ;
	}//Loop over all cells
	var->EField[ 0 ]=var->EField[ 0 ];
	var->EField[ 1 ]=var->EField[ 1 ];
	if( plasma.Mesh.ndim == 3 ) var->EField[ 2 ]=var->EField[ 2 ];
}
void CPoisson::Calculate_NetCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
	Cell *Cell_i ;

	var->NetQ.zero() ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;


		if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ) {

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
	// 
	// double tmp=0.0 ;
	// for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	var->NetQ[ i ] = 0.0 ;

	// 	if ( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

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

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		
		eps = 0.0 ;

		Cell_i  = plasma.get_cell( i ) ;

		if( plasma.get_cell_typename( Cell_i->data_id ) == "PLASMA" ){

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

	// for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

	// 	eps = 0.0 ;
	// 	if( m -> cell[i].type == "PLASMA" ){

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
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

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
double CPoisson::SineVoltage( string FaceType, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int ft=0 ;
	if     ( FaceType == "POWER" ) ft = POWER ;
	else if( FaceType == "GROUND") ft = GROUND ;

	map< int, CElectrical>::iterator Iter;
	Iter = config->ElectricalMap.find( ft ) ;
	return Iter->second.BC_Voltage ;
}
void CPoisson::CalculateDispCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
		var->DispJD[0][i] = var->Eps0[ i ]*( var->EField[ 0 ][ i ] - var->PreEField[ 0 ][ i ] )/var->Dt ;
		var->DispJD[1][i] = var->Eps0[ i ]*( var->EField[ 1 ][ i ] - var->PreEField[ 1 ][ i ] )/var->Dt ;
		if( plasma.Mesh.ndim == 3 ) var->DispJD[2][i] = var->Eps0[ i ]*( var->EField[ 2 ][ i ] - var->PreEField[ 2 ][ i ] )/var->Dt ;

		var->TotalJD[0][i] = var->DispJD[0][ i ] ;
		var->TotalJD[1][i] = var->DispJD[1][ i ] ;
		if( plasma.Mesh.ndim == 3 ) var->TotalJD[2][i] = var->DispJD[2][ i ] ;

	}//Cell Loop
}

#include "domain_structure.hpp"
using namespace std ;
CDomain::CDomain()
{
}
void CDomain::link_ultraMPP()
{
	MPP = &plasma ;
}
void CDomain::SetupCellFaceType()
{
	int iCell, jCell, iFace, j ;
	int nCell = MPP->Mesh.cell_number_with_overlapping_cell ;
	local_cell_number = MPP->Mesh.cell_number ;
	Face *Face_i ;
	Cell *Cell_i ;
	for( int i = 0 ; i < nCell ; i++ ) {

		Cell_i  = MPP->get_cell( i ) ;

			 if ( Cell_i->Typename == "PLASMA"  ) Cell_i->type = PLASMA ;
	else if ( Cell_i->Typename == "POWER"   ) Cell_i->type = POWER 	 ;
	else if ( Cell_i->Typename == "POWER_1" ) Cell_i->type = POWER_1 ;
	else if ( Cell_i->Typename == "POWER_2" ) Cell_i->type = POWER_2 ;
	else if ( Cell_i->Typename == "POWER_3" ) Cell_i->type = POWER_3 ;
	else if ( Cell_i->Typename == "POWER_4" ) Cell_i->type = POWER_4 ;
	else if ( Cell_i->Typename == "POWER_5" ) Cell_i->type = POWER_5 ;
	else if ( Cell_i->Typename == "POWER_6" ) Cell_i->type = POWER_6 ;
	else if ( Cell_i->Typename == "POWER_7" ) Cell_i->type = POWER_7 ;
	else if ( Cell_i->Typename == "POWER_8" ) Cell_i->type = POWER_8 ;
	else if ( Cell_i->Typename == "POWER_9" ) Cell_i->type = POWER_9 ;

		else if ( Cell_i->Typename == "GROUND"   ) Cell_i->type = GROUND   ;
		else if ( Cell_i->Typename == "GROUND_1" ) Cell_i->type = GROUND_1 ;
		else if ( Cell_i->Typename == "GROUND_2" ) Cell_i->type = GROUND_2 ;
		else if ( Cell_i->Typename == "GROUND_3" ) Cell_i->type = GROUND_3 ;
		else if ( Cell_i->Typename == "GROUND_4" ) Cell_i->type = GROUND_4 ;
		else if ( Cell_i->Typename == "GROUND_5" ) Cell_i->type = GROUND_5 ;
		else if ( Cell_i->Typename == "GROUND_6" ) Cell_i->type = GROUND_6 ;
		else if ( Cell_i->Typename == "GROUND_7" ) Cell_i->type = GROUND_7 ;
		else if ( Cell_i->Typename == "GROUND_8" ) Cell_i->type = GROUND_8 ;
		else if ( Cell_i->Typename == "GROUND_9" ) Cell_i->type = GROUND_9 ;

		else if ( Cell_i->Typename == "DIELECTRIC"   ) Cell_i->type = DIELECTRIC ;
		else if ( Cell_i->Typename == "DIELECTRIC_1" ) Cell_i->type = DIELECTRIC_1 ;
		else if ( Cell_i->Typename == "DIELECTRIC_2" ) Cell_i->type = DIELECTRIC_2 ;
		else if ( Cell_i->Typename == "DIELECTRIC_3" ) Cell_i->type = DIELECTRIC_3 ;
		else if ( Cell_i->Typename == "DIELECTRIC_4" ) Cell_i->type = DIELECTRIC_4 ;
		else if ( Cell_i->Typename == "DIELECTRIC_5" ) Cell_i->type = DIELECTRIC_5 ;
		else if ( Cell_i->Typename == "DIELECTRIC_6" ) Cell_i->type = DIELECTRIC_6 ;
		else if ( Cell_i->Typename == "DIELECTRIC_7" ) Cell_i->type = DIELECTRIC_7 ;
		else if ( Cell_i->Typename == "DIELECTRIC_8" ) Cell_i->type = DIELECTRIC_8 ;
		else if ( Cell_i->Typename == "DIELECTRIC_9" ) Cell_i->type = DIELECTRIC_9 ;
		else if ( Cell_i->Typename == "EMITTER" 	   ) Cell_i->type = EMITTER ;
		else { cout<<"Material Error"<<endl; exit(1) ;}
		//cout<<"TypeName: "<<Cell_i->Typename<<"\t"<<"Cell type: "<<Cell_i->type<<endl;
	}
	 for ( int i = 0 ; i < nCell ; i++ ) {

	 	Cell_i  = MPP->get_cell( i ) ;

	 	iFace 	 = Cell_i->face_number ;
	 	iCell 	 = Cell_i->cell_number ;

		for ( int j = iCell ; j < iFace ; j++ ) {

			Face_i = Cell_i->face[ j ] ;
			//cout<<"Typename: "<<Face_i->Typename<<endl;
				 if ( Face_i->Typename ==   "POWER" ) Face_i->type = POWER ;	 
			else if ( Face_i->Typename == "POWER_1" ) Face_i->type = POWER_1 ;
			else if ( Face_i->Typename == "POWER_2" ) Face_i->type = POWER_2 ;
			else if ( Face_i->Typename == "POWER_3" ) Face_i->type = POWER_3 ;
			else if ( Face_i->Typename == "POWER_4" ) Face_i->type = POWER_4 ;
			else if ( Face_i->Typename == "POWER_5" ) Face_i->type = POWER_5 ;
			else if ( Face_i->Typename == "POWER_6" ) Face_i->type = POWER_6 ;
			else if ( Face_i->Typename == "POWER_7" ) Face_i->type = POWER_7 ;
			else if ( Face_i->Typename == "POWER_8" ) Face_i->type = POWER_8 ;
			else if ( Face_i->Typename == "POWER_9" ) Face_i->type = POWER_9 ;

			else if ( Face_i->Typename ==   "GROUND" ) Face_i->type = GROUND ;	 
			else if ( Face_i->Typename == "GROUND_1" ) Face_i->type = GROUND_1 ;
			else if ( Face_i->Typename == "GROUND_2" ) Face_i->type = GROUND_2 ;
			else if ( Face_i->Typename == "GROUND_3" ) Face_i->type = GROUND_3 ;
			else if ( Face_i->Typename == "GROUND_4" ) Face_i->type = GROUND_4 ;
			else if ( Face_i->Typename == "GROUND_5" ) Face_i->type = GROUND_5 ;
			else if ( Face_i->Typename == "GROUND_6" ) Face_i->type = GROUND_6 ;
			else if ( Face_i->Typename == "GROUND_7" ) Face_i->type = GROUND_7 ;
			else if ( Face_i->Typename == "GROUND_8" ) Face_i->type = GROUND_8 ;
			else if ( Face_i->Typename == "GROUND_9" ) Face_i->type = GROUND_9 ;

			else if ( Face_i->Typename ==   "DIELECTRIC" ) Face_i->type = DIELECTRIC   ;	 
			else if ( Face_i->Typename == "DIELECTRIC_1" ) Face_i->type = DIELECTRIC_1 ;
			else if ( Face_i->Typename == "DIELECTRIC_2" ) Face_i->type = DIELECTRIC_2 ;
			else if ( Face_i->Typename == "DIELECTRIC_3" ) Face_i->type = DIELECTRIC_3 ;
			else if ( Face_i->Typename == "DIELECTRIC_4" ) Face_i->type = DIELECTRIC_4 ;
			else if ( Face_i->Typename == "DIELECTRIC_5" ) Face_i->type = DIELECTRIC_5 ;
			else if ( Face_i->Typename == "DIELECTRIC_6" ) Face_i->type = DIELECTRIC_6 ;
			else if ( Face_i->Typename == "DIELECTRIC_7" ) Face_i->type = DIELECTRIC_7 ;
			else if ( Face_i->Typename == "DIELECTRIC_8" ) Face_i->type = DIELECTRIC_8 ;
			else if ( Face_i->Typename == "DIELECTRIC_9" ) Face_i->type = DIELECTRIC_9 ;


			else if ( Face_i->Typename == "NEUMANN"  ) Face_i->type = NEUMANN ;
			else if ( Face_i->Typename == "OPEN" 	   ) Face_i->type = NEUMANN ;
			else if ( Face_i->Typename == "SYMMETRY" ) Face_i->type = NEUMANN ;
			else if ( Face_i->Typename == "EMITTER"  ) Face_i->type = EMITTER ;
			else { cout<<"Boundary Face Error"<<endl; exit(1) ; }
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
void CDomain::BulidCellStructure()
{
	int iCell, iFace ;
	int nCell = MPP->Mesh.cell_number ;

	nDim = MPP->Mesh.ndim ;

	Cell *Cell_i ;
	//cout<<"AAA"<<endl;
	PFM_CELL = new CCell *[ nCell ] ;
	//cout<<"BBB"<<endl;

	for( int i = 0 ; i < nCell ; i++ ) {

		Cell_i  = MPP->get_cell( i ) ;

		iCell 	 = Cell_i->cell_number ;
		iFace 	 = Cell_i->face_number ;
		//cout<<"CCC"<<endl;
		//cout<<"i:"<<i<<", iFace: "<<iFace<<endl;
		PFM_CELL[ i ] = new CCell [ iFace ] ;
		//cout<<"DDD"<<endl;

		//cout<<"EEE"<<endl;
		for( int k = 0 ; k < iFace ; k++ ) PFM_CELL[ i ][ k ].Init() ;
		//cout<<"FFF"<<endl;

		/*--- For interior cell ---*/
		for( int k = 0 ; k < iCell ; k++ ){

			PFM_CELL[ i ][ k ].NeighborCellId 			= Cell_i->cell[ k ]->local_id ;
			PFM_CELL[ i ][ k ].NeighborGlobalCellId = Cell_i->cell[ k ]->id ;

			if ( Cell_i->type == Cell_i->cell[ k ]->type ){

				/*--- PN Vector ---*/
				PN[ 0 ] = Cell_i->cell[ k ]->r[0] - Cell_i->r[0] ;
				PN[ 1 ] = Cell_i->cell[ k ]->r[1] - Cell_i->r[1] ;

				/*--- Pf Vector ---*/
				Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
				Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

				/*--- Nf Vector ---*/
				Nf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->cell[ k ]->r[0] ;
				Nf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->cell[ k ]->r[1] ;

				/*--- Face Normal & Tangent Vector ---*/
				PFM_CELL[ i ][ k ].nf[ 0 ] = Cell_i->nA[ k ][ 0 ] ;
				PFM_CELL[ i ][ k ].nf[ 1 ] = Cell_i->nA[ k ][ 1 ] ;
				PFM_CELL[ i ][ k ].nf[ 2 ] = Cell_i->nA[ k ][ 2 ] ;

				PFM_CELL[ i ][ k ].mf[ 0 ] = -Cell_i->nA[ k ][ 1 ] ;
				PFM_CELL[ i ][ k ].mf[ 1 ] =  Cell_i->nA[ k ][ 0 ] ;

				/*--- P'f Vector ---*/
				PPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[0] ;
				PPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[1] ;
			    PFM_CELL[ i ][ k ].dPPf = sqrt( DotProduct( PPf, PPf ) ) ; 

				/*--- N'f Vector ---*/
				NPf[ 0 ] = DotProduct( Nf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[0];
				NPf[ 1 ] = DotProduct( Nf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[1];
			  PFM_CELL[ i ][ k ].dNPf = sqrt( DotProduct( NPf, NPf ) ) ; 

				/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
				PFM_CELL[ i ][ k ].PPP[ 0 ] = Pf[ 0 ] - PPf[ 0 ] ;
				PFM_CELL[ i ][ k ].PPP[ 1 ] = Pf[ 1 ] - PPf[ 1 ] ;
				// if ( fabs( Cell[ i ][ k ].PPP[ 0 ] ) < ZERO ) Cell[ i ][ k ].PPP[ 0 ] = 0.0 ;
				// if ( fabs( Cell[ i ][ k ].PPP[ 1 ] ) < ZERO ) Cell[ i ][ k ].PPP[ 1 ] = 0.0 ;

				/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
				PFM_CELL[ i ][ k ].NNP[ 0 ] = Nf[ 0 ] - NPf[ 0 ] ;
				PFM_CELL[ i ][ k ].NNP[ 1 ] = Nf[ 1 ] - NPf[ 1 ] ;
				// if ( fabs( Cell[ i ][ k ].NNP[ 0 ] ) < ZERO ) Cell[ i ][ k ].NNP[ 0 ] = 0.0 ;
				// if ( fabs( Cell[ i ][ k ].NNP[ 1 ] ) < ZERO ) Cell[ i ][ k ].NNP[ 1 ] = 0.0 ;

				/*--- Af Vector ---*/
				PFM_CELL[ i ][ k ].Af[ 0 ] = Cell_i->A[ k ][ 0 ] ;
				PFM_CELL[ i ][ k ].Af[ 1 ] = Cell_i->A[ k ][ 1 ] ;

				/*--- dDist = |P'N'|, dArea = face area ---*/
				PFM_CELL[ i ][ k ].dDist = PFM_CELL[ i ][ k ].dPPf + PFM_CELL[ i ][ k ].dNPf ;
				PFM_CELL[ i ][ k ].dArea = Cell_i->face[ k ]->dA ;

				PFM_CELL[ i ][ k ].SurfaceCharge = 0.0 ;

			} else {//Discontinue face

				/*--- PN Vector ---*/
				PN[ 0 ] = Cell_i->cell[ k ]->r[0] - Cell_i->r[0] ;
				PN[ 1 ] = Cell_i->cell[ k ]->r[1] - Cell_i->r[1] ;

				/*--- Pf Vector ---*/
				Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
				Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

				/*--- Nf Vector ---*/
				Nf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->cell[ k ]->r[0] ;
				Nf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->cell[ k ]->r[1] ;

				/*--- Face Normal & Tangent Vector ---*/
				PFM_CELL[ i ][ k ].nf[ 0 ] = Cell_i->nA[ k ][ 0 ] ;
				PFM_CELL[ i ][ k ].nf[ 1 ] = Cell_i->nA[ k ][ 1 ] ;
				PFM_CELL[ i ][ k ].nf[ 2 ] = Cell_i->nA[ k ][ 2 ] ;

				PFM_CELL[ i ][ k ].mf[ 0 ] = -Cell_i->nA[ k ][ 1 ] ;
				PFM_CELL[ i ][ k ].mf[ 1 ] =  Cell_i->nA[ k ][ 0 ] ;

				/*--- P'f Vector ---*/
				PPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[0] ;
				PPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[1] ;
				PFM_CELL[ i ][ k ].dPPf = sqrt( DotProduct( PPf, PPf ) ) ; 

				/*--- N'f Vector ---*/
				NPf[ 0 ] = DotProduct( Nf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[0];
				NPf[ 1 ] = DotProduct( Nf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[1];
				PFM_CELL[ i ][ k ].dNPf = sqrt( DotProduct( NPf, NPf ) ) ; 

				/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
				PFM_CELL[ i ][ k ].PPP[ 0 ] = Pf[ 0 ] - PPf[ 0 ] ;
				PFM_CELL[ i ][ k ].PPP[ 1 ] = Pf[ 1 ] - PPf[ 1 ] ;
				// if ( fabs( PFM_CELL[ i ][ k ].PPP[ 0 ] ) < ZERO ) PFM_CELL[ i ][ k ].PPP[ 0 ] = 0.0 ;
				// if ( fabs( PFM_CELL[ i ][ k ].PPP[ 1 ] ) < ZERO ) PFM_CELL[ i ][ k ].PPP[ 1 ] = 0.0 ;

				/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
				PFM_CELL[ i ][ k ].NNP[ 0 ] = Nf[ 0 ] - NPf[ 0 ] ;
				PFM_CELL[ i ][ k ].NNP[ 1 ] = Nf[ 1 ] - NPf[ 1 ] ;
				// if ( fabs( PFM_CELL[ i ][ k ].NNP[ 0 ] ) < ZERO ) PFM_CELL[ i ][ k ].NNP[ 0 ] = 0.0 ;
				// if ( fabs( PFM_CELL[ i ][ k ].NNP[ 1 ] ) < ZERO ) PFM_CELL[ i ][ k ].NNP[ 1 ] = 0.0 ;

				/*--- Af Vector ---*/
				PFM_CELL[ i ][ k ].Af[ 0 ] = Cell_i->A[ k ][ 0 ] ;
				PFM_CELL[ i ][ k ].Af[ 1 ] = Cell_i->A[ k ][ 1 ] ;

				/*--- dDist = |P'N'|, dArea = face area ---*/
				PFM_CELL[ i ][ k ].dDist = PFM_CELL[ i ][ k ].dPPf ;//+ PFM_CELL[ i ][ k ].dNPf ;
				PFM_CELL[ i ][ k ].dArea = Cell_i->face[ k ]->dA ;
				PFM_CELL[ i ][ k ].SurfaceCharge = 0.0 ;
			}
		}//End Interior cell

		/*--- For boundary face ---*/
		for( int k = iCell ; k < iFace ; k++ ) {

			PFM_CELL[ i ][ k ].NeighborCellId 			= -999 ;
			PFM_CELL[ i ][ k ].NeighborGlobalCellId = -999 ;

			/*--- PN Vector ---*/
			//PN[ 0 ] = Cell_i->cell[ k ]->r[0] - Cell_i->r[0] ;
			//PN[ 1 ] = Cell_i->cell[ k ]->r[1] - Cell_i->r[1] ;

			/*--- Pf Vector ---*/
			Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
			Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

			/*--- Nf Vector ---*/
			//Nf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->cell[ k ]->r[0] ;
			//Nf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->cell[ k ]->r[1] ;

			/*--- Face Normal & Tangent Vector ---*/
			PFM_CELL[ i ][ k ].nf[ 0 ] = Cell_i->nA[ k ][ 0 ] ;
			PFM_CELL[ i ][ k ].nf[ 1 ] = Cell_i->nA[ k ][ 1 ] ;
			PFM_CELL[ i ][ k ].nf[ 2 ] = Cell_i->nA[ k ][ 2 ] ;

			PFM_CELL[ i ][ k ].mf[ 0 ] = -Cell_i->nA[ k ][ 1 ] ;
			PFM_CELL[ i ][ k ].mf[ 1 ] =  Cell_i->nA[ k ][ 0 ] ;

			/*--- P'f  Normal Vector ---*/
			PPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[0] ;
			PPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[1] ;
		  PFM_CELL[ i ][ k ].dPPf = sqrt( DotProduct( PPf, PPf ) ) ; 

			/*--- N'f  Normal Vector ---*/
			//NPf[ 0 ] = DotProduct( Nf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[0];
			//NPf[ 1 ] = DotProduct( Nf, PFM_CELL[ i ][ k ].nf )*PFM_CELL[ i ][ k ].nf[1];
		   	//PFM_CELL[ i ][ k ].dNPf = sqrt( DotProduct( NPf, NPf ) ) ; 

			/*--- PP'  Normal Vector ( PP' = Pf - P'f ) ---*/
			PFM_CELL[ i ][ k ].PPP[ 0 ] = Pf[ 0 ] - PPf[ 0 ] ;
			PFM_CELL[ i ][ k ].PPP[ 1 ] = Pf[ 1 ] - PPf[ 1 ] ;
			// if ( fabs( PFM_CELL[ i ][ k ].PPP[ 0 ] ) < ZERO ) PFM_CELL[ i ][ k ].PPP[ 0 ] = 0.0 ;
			// if ( fabs( PFM_CELL[ i ][ k ].PPP[ 1 ] ) < ZERO ) PFM_CELL[ i ][ k ].PPP[ 1 ] = 0.0 ;
			/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
			//PFM_CELL[ i ][ k ].NNP[ 0 ] = Nf[ 0 ] - NPf[ 0 ] ;
			//PFM_CELL[ i ][ k ].NNP[ 1 ] = Nf[ 1 ] - NPf[ 1 ] ;

			/*--- Af Vector ---*/
			PFM_CELL[ i ][ k ].Af[ 0 ] = Cell_i->A[ k ][ 0 ] ;
			PFM_CELL[ i ][ k ].Af[ 1 ] = Cell_i->A[ k ][ 1 ] ;

			/*--- dDist = |P'N'|, dArea = face area ---*/
			PFM_CELL[ i ][ k ].dDist = PFM_CELL[ i ][ k ].dPPf ; //+ PFM_CELL[ i ][ k ].dNPf ;
			PFM_CELL[ i ][ k ].dArea = Cell_i->face[ k ]->dA ;

			PFM_CELL[ i ][ k ].SurfaceCharge = 0.0 ;
		}//End boundary face

	}//Cell Loop
	//cout<<"FUCK"<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	//cout<<"FUCK0"<<endl;
 	//Calculate_LSQ_Coeff() ;
}

void CDomain::Calculate_LSQ_Coeff()
{
 
	int nCell = MPP->Mesh.cell_number_with_overlapping_cell ;
	int iFace=0, iCell=0 ;
	double dx=0.0, dy=0.0, a11=0.0, a12=0.0, a21=0.0, a22=0.0, det=0.0 ;
	double ia11=0.0, ia12=0.0, ia21=0.0, ia22=0.0 ;
	LSQ = new CLSQ [ nCell ] ;
	Cell *Cell_i ;
	for ( int i = 0 ; i < nCell ; i++ ) {
		LSQ[ i ].Init() ;
		Cell_i = MPP->get_cell(i) ;

		iCell = Cell_i->cell_number ;
		iFace = Cell_i->face_number ;
		/*--- Reset Matrix ---*/
		a11 = 0.0 ; a12 = 0.0 ;
		a21 = 0.0 ; a22 = 0.0 ;

		/*--- Loop over neighbor "cells" ---*/
		for ( int k = 0 ; k < iCell ; k++ ) {

			if ( Cell_i->type != Cell_i->cell[ k ]->type ) {//For discontinued face

				 Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
				 Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
				fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
				fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
				dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
				dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;
				//dx = Cell_i->face[ k ]->r[0]  - Cell_i->r[0] ;
				//dy = Cell_i->face[ k ]->r[1]  - Cell_i->r[1] ;
			} else {

				dx = Cell_i->cell[ k ]->r[0]  - Cell_i->r[0] ;
				dy = Cell_i->cell[ k ]->r[1]  - Cell_i->r[1] ;
			}

			a11 = a11 + LSQ_Weight( dx, dy )*dx*dx ; 
			a12 = a12 + LSQ_Weight( dx, dy )*dx*dy ;
			a21 = a21 + LSQ_Weight( dx, dy )*dx*dy ; 
			a22 = a22 + LSQ_Weight( dx, dy )*dy*dy ;

		}

		/*--- Loop over domain boundary "faces" ---*/
		for ( int k = iCell ; k < iFace ; k++ ) {

			Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
			Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
			fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
			fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
			dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
			dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

			//dx = Cell_i->face[ k ]->r[0]  - Cell_i->r[0] ;
			//dy = Cell_i->face[ k ]->r[1]  - Cell_i->r[1] ;

			a11 = a11 + LSQ_Weight( dx, dy )*dx*dx ; 	
			a12 = a12 + LSQ_Weight( dx, dy )*dx*dy ;
			a21 = a21 + LSQ_Weight( dx, dy )*dx*dy ;	
			a22 = a22 + LSQ_Weight( dx, dy )*dy*dy ;
		}

		/*--- Cal. LSQ det. ---*/
		det = a11*a22 - a12*a21 ;
		//if( fabs(det) < 1.E-14 ) cout<<"LSQ det Singular"<<endl;

		/*--- Cal. Inverse Matrix ---*/
		ia11 =  a22/det ;
		ia12 = -a21/det ;
		ia21 = -a12/det ;
		ia22 =  a11/det ;

		/*--- Cal. LSQ Coefficient for "cells" ---*/
		for ( int k = 0 ; k < iCell ; k++ ) {

			if ( Cell_i->type != Cell_i->cell[ k ]->type ) {

				Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
				Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
				fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
				fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
				dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
				dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

				//dx =  Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
				//dy =  Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

			} else {

				dx = Cell_i->cell[ k ]->r[0]  - Cell_i->r[0] ;
				dy = Cell_i->cell[ k ]->r[1]  - Cell_i->r[1] ;
			}
			LSQ[ i ].Cx[ k ] = LSQ_Weight( dx, dy )*ia11*dx + LSQ_Weight( dx, dy )*ia12*dy ;
			LSQ[ i ].Cy[ k ] = LSQ_Weight( dx, dy )*ia21*dx + LSQ_Weight( dx, dy )*ia22*dy ;
		}
		//if(mpi_rank==0) cout<<"---------------------------------------"<<endl;
		/*--- Cal. LSQ Coefficient for domain boundary "faces" ---*/
		for ( int k = iCell ; k < iFace ; k++ ) {

			Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
			Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
			fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
			fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
			dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
			dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

			//dx =  Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
			//dy =  Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

			LSQ[ i ].Cx[ k ] = LSQ_Weight( dx, dy )*ia11*dx + LSQ_Weight( dx, dy )*ia12*dy ;
			LSQ[ i ].Cy[ k ] = LSQ_Weight( dx, dy )*ia21*dx + LSQ_Weight( dx, dy )*ia22*dy ;
		}//End boundaty face
	}//End cell loop
	MPI_Barrier(MPI_COMM_WORLD);
}

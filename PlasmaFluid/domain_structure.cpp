#include "domain_structure.hpp"
using namespace std ;
CDomain::CDomain()
{

}
void CDomain::BulidCellStructure()
{
	nDim = plasma.Mesh.ndim ;

	Cell *Cell_i, *Cell_j ;

	PFM_CELL = new CCell *[ plasma.Mesh.cell_number ] ;

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell_i  = plasma.get_cell( i ) ;

		PFM_CELL[ i ] = new CCell [ Cell_i->face_number ] ;

		for( int k = 0 ; k < Cell_i->face_number ; k++ ) PFM_CELL[ i ][ k ].Init() ;

		/*--- For interior cell ---*/
		for( int k = 0 ; k < Cell_i->cell_number ; k++ ){

			PFM_CELL[ i ][ k ].NeighborCellId 			= Cell_i->cell[ k ]->local_id ;
			PFM_CELL[ i ][ k ].NeighborGlobalCellId = Cell_i->cell[ k ]->id ;

			if ( plasma.get_cell_typename( Cell_i->data_id ) ==  plasma.get_cell_typename( Cell_i->cell[ k ]->data_id ) ){

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
		for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

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
	MPI_Barrier(MPI_COMM_WORLD);
}

void CDomain::Calculate_LSQ_Coeff()
{
 
	// double dx=0.0, dy=0.0, a11=0.0, a12=0.0, a21=0.0, a22=0.0, det=0.0 ;
	// double ia11=0.0, ia12=0.0, ia21=0.0, ia22=0.0 ;
	// LSQ = new CLSQ [ plasma.Mesh.cell_number ] ;
	// Cell *Cell_i ;


	// for ( int i = 0 ; i < plasma.Mesh.cell_number_with_overlapping_cell ; i++ ) {

	// 	LSQ[ i ].Init() ;

	// 	Cell_i = plasma.get_cell(i) ;

	// 	/*--- Reset Matrix ---*/
	// 	a11 = 0.0 ; a12 = 0.0 ;
	// 	a21 = 0.0 ; a22 = 0.0 ;

	// 	/*--- Loop over neighbor "cells" ---*/
	// 	for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

	// 		if ( Cell_i->type != Cell_i->cell[ k ]->type ) {//For discontinued face

	// 			 Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
	// 			 Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
	// 			fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
	// 			fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
	// 			dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
	// 			dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;
	// 			//dx = Cell_i->face[ k ]->r[0]  - Cell_i->r[0] ;
	// 			//dy = Cell_i->face[ k ]->r[1]  - Cell_i->r[1] ;
	// 		} else {

	// 			dx = Cell_i->cell[ k ]->r[0]  - Cell_i->r[0] ;
	// 			dy = Cell_i->cell[ k ]->r[1]  - Cell_i->r[1] ;
	// 		}

	// 		a11 = a11 + LSQ_Weight( dx, dy )*dx*dx ; 
	// 		a12 = a12 + LSQ_Weight( dx, dy )*dx*dy ;
	// 		a21 = a21 + LSQ_Weight( dx, dy )*dx*dy ; 
	// 		a22 = a22 + LSQ_Weight( dx, dy )*dy*dy ;

	// 	}

	// 	/*--- Loop over domain boundary "faces" ---*/
	// 	for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 		Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
	// 		Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
	// 		fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
	// 		fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
	// 		dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
	// 		dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

	// 		//dx = Cell_i->face[ k ]->r[0]  - Cell_i->r[0] ;
	// 		//dy = Cell_i->face[ k ]->r[1]  - Cell_i->r[1] ;

	// 		a11 = a11 + LSQ_Weight( dx, dy )*dx*dx ; 	
	// 		a12 = a12 + LSQ_Weight( dx, dy )*dx*dy ;
	// 		a21 = a21 + LSQ_Weight( dx, dy )*dx*dy ;	
	// 		a22 = a22 + LSQ_Weight( dx, dy )*dy*dy ;
	// 	}

	// 	/*--- Cal. LSQ det. ---*/
	// 	det = a11*a22 - a12*a21 ;
	// 	//if( fabs(det) < 1.E-14 ) cout<<"LSQ det Singular"<<endl;

	// 	/*--- Cal. Inverse Matrix ---*/
	// 	ia11 =  a22/det ;
	// 	ia12 = -a21/det ;
	// 	ia21 = -a12/det ;
	// 	ia22 =  a11/det ;

	// 	/*--- Cal. LSQ Coefficient for "cells" ---*/
	// 	for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

	// 		if ( Cell_i->type != Cell_i->cell[ k ]->type ) {

	// 			Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
	// 			Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
	// 			fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
	// 			fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
	// 			dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
	// 			dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

	// 			//dx =  Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
	// 			//dy =  Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

	// 		} else {

	// 			dx = Cell_i->cell[ k ]->r[0]  - Cell_i->r[0] ;
	// 			dy = Cell_i->cell[ k ]->r[1]  - Cell_i->r[1] ;
	// 		}
	// 		LSQ[ i ].Cx[ k ] = LSQ_Weight( dx, dy )*ia11*dx + LSQ_Weight( dx, dy )*ia12*dy ;
	// 		LSQ[ i ].Cy[ k ] = LSQ_Weight( dx, dy )*ia21*dx + LSQ_Weight( dx, dy )*ia22*dy ;
	// 	}
	// 	//if(mpi_rank==0) cout<<"---------------------------------------"<<endl;
	// 	/*--- Cal. LSQ Coefficient for domain boundary "faces" ---*/
	// 	for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

	// 		Pf[ 0 ] = Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
	// 		Pf[ 1 ] = Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;
	// 		fPf[ 0 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[0] ;
	// 		fPf[ 1 ] = DotProduct( Pf, PFM_CELL[ i ][ k ].mf )*PFM_CELL[ i ][ k ].mf[1] ;
	// 		dx = ( -fPf[ 0 ] + Cell_i->face[ k ]->r[0] )  - Cell_i->r[0] ;
	// 		dy = ( -fPf[ 1 ] + Cell_i->face[ k ]->r[1] )  - Cell_i->r[1] ;

	// 		//dx =  Cell_i->face[ k ]->r[0] - Cell_i->r[0] ;
	// 		//dy =  Cell_i->face[ k ]->r[1] - Cell_i->r[1] ;

	// 		LSQ[ i ].Cx[ k ] = LSQ_Weight( dx, dy )*ia11*dx + LSQ_Weight( dx, dy )*ia12*dy ;
	// 		LSQ[ i ].Cy[ k ] = LSQ_Weight( dx, dy )*ia21*dx + LSQ_Weight( dx, dy )*ia22*dy ;
	// 	}//End boundaty face
	// }//End cell loop
	// MPI_Barrier(MPI_COMM_WORLD);
}
void CDomain::Init()
{
	/* Face */
	type_typename[  POWER  ] = "POWER" ;
	typename_type[ "POWER" ] = 	POWER  ;

	type_typename[  GROUND  ]	=	"GROUND" ;
	typename_type[ "GROUND" ]	=	 GROUND  ;

	type_typename[  SYMMETRIC  ]	=	"SYMMETRIC" ;
	typename_type[ "SYMMETRIC" ]	=	 SYMMETRIC  ;


	/* Cell */
	type_typename[  SOLID_POWER  ]	=	"SOLID_POWER" ;
	typename_type[ "SOLID_POWER" ]	=	 SOLID_POWER  ;

	type_typename[  SOLID_GROUND  ]	=	"SOLID_GROUND" ;
	typename_type[ "SOLID_GROUND" ]	=	 SOLID_GROUND  ;


	type_typename[  PLASMA  ]	=	"PLASMA" ;
	typename_type[ "PLASMA" ]	=	 PLASMA  ;

	type_typename[  DIELECTRIC  ]	=	"DIELECTRIC";
	typename_type[ "DIELECTRIC" ]	=	 DIELECTRIC ;

	UltraMPPExtractFaceCellTag() ;

	/* Face */
	for ( auto mpp = MPP_face_tag.cbegin(); mpp != MPP_face_tag.cend(); ++mpp ) {
		for ( auto pfm = typename_type.cbegin(); pfm != typename_type.cend(); ++pfm ){
			if ( mpp->first == pfm-> first ) {
				face_type[ mpp->second ] = pfm->second ;
			}
		}//PFM map
	}//MPP map

	/* Cell */
	for ( auto mpp = MPP_cell_tag.cbegin(); mpp != MPP_cell_tag.cend(); ++mpp ) {
		for ( auto pfm = typename_type.cbegin(); pfm != typename_type.cend(); ++pfm ){
			if ( mpp->first == pfm-> first ) {
				cell_type[ mpp->second ] = pfm->second ;
			}
		}//PFM map
	}//MPP map

	#define mpp_pfm_check false
	#if( mpp_pfm_check == true )
		for ( auto mpp_pfm = face_type.cbegin(); mpp_pfm != face_type.cend(); ++mpp_pfm ){
			cout<< mpp_pfm->first << "\t" << mpp_pfm->second << endl;
		}
		for ( auto mpp_pfm = cell_type.cbegin(); mpp_pfm != cell_type.cend(); ++mpp_pfm ){
			cout<< mpp_pfm->first << "\t" << mpp_pfm->second << endl;
		}
	#endif
}
void CDomain::UltraMPPExtractFaceCellTag()
{
	/* Extract the boundary & cell setting from input json file. */
    json &json_bc_setting    = *(plasma.get_json_input_parameter("boundary_setting")) ;
    json &json_cell_setting  = *(plasma.get_json_input_parameter("volume_setting"  )) ;

    /* Face */
		for ( int ibc = 0; ibc < json_bc_setting["name"].size(); ibc++){
			face_parameter[ json_bc_setting["name"][ ibc ] ] = json_bc_setting["values"][ ibc ] ;
			MPP_face_tag  [ json_bc_setting["name"][ ibc ] ] = plasma.set_bc_mapping(  json_bc_setting["name"][ ibc ],  json_bc_setting["boundary_type"][ ibc ]  );
		}

		/* Cell */
	  for ( int icc = 0; icc < json_cell_setting["name"].size(); icc++ ) {
	  	cell_parameter[ json_cell_setting["name"][ icc ] ] = json_cell_setting["permittivity"][ icc ] ;
	  	MPP_cell_tag  [ json_cell_setting["name"][ icc ] ] = plasma.get_cell_type_mapping( json_cell_setting["name"][ icc ] ) ;
	  }

	  /* Debug */
		#define face_cell_check false
    #if( face_cell_check == true )
	    for ( int ibc = 0; ibc < json_bc_setting["name"].size(); ibc++) {
	      cout << json_bc_setting  ["name"][ ibc ] <<"\t" << json_bc_setting  ["boundary_type"][ ibc ] << "\t" << json_bc_setting  ["values"][ ibc ]<<endl ;
	    }
	    for ( int ibc = 0; ibc < json_cell_setting["name"].size(); ibc++){
	      cout << json_cell_setting["name"][ ibc ] <<"\t" << json_cell_setting["volume_type"  ][ ibc ] << "\t" << json_cell_setting["permittivity"][ ibc ] <<endl ;
	    }
    #endif
}
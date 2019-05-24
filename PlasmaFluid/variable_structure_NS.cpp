
#include "variable_structure_NS.hpp"
using namespace std ;
CVariableNS::CVariableNS()
{
}
void CVariableNS::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	U0.initial( m, CELL, "Rho"  ) ;
	U1.initial( m, CELL, "Rho_U") ;
	U2.initial( m, CELL, "Rho_V") ;
	U3.initial( m, CELL, "Rho_W") ;
	U4.initial( m, CELL, "RHO_E") ;
	Dt.initial( m, CELL, "DT" 	) ;

	PreU0.initial( m, CELL, "Rho"  ) ;
	PreU1.initial( m, CELL, "Rho_U") ;
	PreU2.initial( m, CELL, "Rho_W") ;
	PreU3.initial( m, CELL, "Rho_W") ;
	PreU4.initial( m, CELL, "RHO_E") ;
	MachNumber.initial( m, CELL, "M") ;
}
void CVariableNS::InitialCondition( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, int iSpecies )
{
	int nCell = m->local_cell_number ; 
	double Mass = config->Species[iSpecies].Mass_Kg ;
	double Rho_inf 	= 1.0, 
		   M_inf 	= 1.65, 
		   u_inf 	= M_inf, 
		   v_inf = 0.0, 
		   P_inf = 1.0/config->Species[ iSpecies ].Gamma, 
		   a_inf = sqrt(config->Species[ iSpecies ].Gamma*P_inf/Rho_inf), 
		   H_inf = a_inf*a_inf/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*( u_inf*u_inf + v_inf*v_inf), 
		   Entropy_Inf = pow(Rho_inf,config->Species[ iSpecies ].Gamma)/P_inf ;
		   
	for ( int i = 0 ; i < nCell ; i ++ )
	{
		if(m->cell[i].type == PLASMA ){

			U0[i]	=	Rho_inf ;
			U1[i]	=	Rho_inf*u_inf ;
			U2[i]	=	0.0;
			U3[i]	=	0.0;
			U4[i]	=	P_inf/(config->Species[ iSpecies ].Gamma-1.0)+0.5*Rho_inf*(u_inf*u_inf+v_inf*v_inf) ;
		}
	}
	U0 = U0 ;
	U1 = U1 ;
	U2 = U2 ;
	U3 = U3 ;
	U4 = U4 ;
	UpdateSolution( m, config ) ;
}
void CVariableNS::UpdateSolution( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	int nCell = m->local_cell_number ;
	for ( int i = 0 ; i < nCell ; i ++ )
	{
		if(m->cell[i].type == PLASMA ){

			PreU0[i]	=	U0[i] ;
			PreU1[i]	=	U1[i] ;
			PreU2[i]	=	U2[i] ;
			PreU3[i]	=	U3[i] ;
			PreU4[i]	=	U4[i] ;
		}
	}
	PreU0 = PreU0 ;
	PreU1 = PreU1 ;
	PreU2 = PreU2 ;
	PreU3 = PreU3 ;
	PreU4 = PreU4 ;

}
void CVariableNS::Calculate_LSQ_Coeff_Scalar( boost::shared_ptr<CDomain> &m )
{

	for ( int k = 0 ; k < 6 ; k++ ) {
		LSQ_Cx[ k ].initial( m, CELL, "LSQ_Cx" ) ;			
		LSQ_Cy[ k ].initial( m, CELL, "LSQ_Cy" ) ;			
		LSQ_Cz[ k ].initial( m, CELL, "LSQ_Cz" ) ;			
	} 
	int nCell = m->local_cell_number ;

	int iFace=0, iCell=0 ;
	double dx=0.0, dy=0.0, a11=0.0, a12=0.0, a21=0.0, a22=0.0, det=0.0 ;
	double ia11=0.0, ia12=0.0, ia21=0.0, ia22=0.0 ;

	for ( int i = 0 ; i < nCell ; i++ ) {

		
		iCell = m->cell[ i ].cell_number ;
		iFace = m->cell[ i ].face_number ;

		/*--- Reset Matrix ---*/
		a11 = 0.0 ; a12 = 0.0 ;
		a21 = 0.0 ; a22 = 0.0 ;

		/*--- Loop over neighbor "cells" ---*/
		for ( int k = 0 ; k < iCell ; k++ ) {

			if ( m->cell[ i ].type != m->cell[ i ].cell[ k ]->type ) {//For discontinued face

				Pf[ 0 ] = m->cell[ i ].face[ k ]->x - m->cell[ i ].x ;
				Pf[ 1 ] = m->cell[ i ].face[ k ]->y - m->cell[ i ].y ;
				//cout<<Cell[ i ][ k ].nf[0]<<endl;
				//cout<<Cell[ i ][ k ].nf[1]<<endl;
				fPf[ 0 ] = DotProduct( Pf, m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[0] ;
				fPf[ 1 ] = DotProduct( Pf, m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[1] ;
				dx = ( -fPf[ 0 ] + m->cell[ i ].face[ k ]->x )  - m->cell[ i ].x ;
				dy = ( -fPf[ 1 ] + m->cell[ i ].face[ k ]->y )  - m->cell[ i ].y ;
				//dx = cell[ i ].face[ k ]->x  - cell[ i ].x ;
				//dy = cell[ i ].face[ k ]->y  - cell[ i ].y ;
			} else {

				dx = m->cell[ i ].cell[ k ]->x  - m->cell[ i ].x ;
				dy = m->cell[ i ].cell[ k ]->y  - m->cell[ i ].y ;
			}

			a11 = a11 + dx*dx ; 
			a12 = a12 + dx*dy ;
			a21 = a21 + dx*dy ; 
			a22 = a22 + dy*dy ;

		}

		/*--- Loop over domain boundary "faces" ---*/
		for ( int k = iCell ; k < iFace ; k++ ) {

			Pf[ 0 ] = m->cell[ i ].face[ k ]->x - m->cell[ i ].x ;
			Pf[ 1 ] = m->cell[ i ].face[ k ]->y - m->cell[ i ].y ;
			fPf[ 0 ] = DotProduct( Pf, m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[0] ;
			fPf[ 1 ] = DotProduct( Pf, m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[1] ;
			dx = ( -fPf[ 0 ] + m->cell[ i ].face[ k ]->x )  - m->cell[ i ].x ;
			dy = ( -fPf[ 1 ] + m->cell[ i ].face[ k ]->y )  - m->cell[ i ].y ;

			//dx = cell[ i ].face[ k ]->x  - cell[ i ].x ;
			//dy = cell[ i ].face[ k ]->y  - cell[ i ].y ;

			a11 = a11 + dx*dx ; 	
			a12 = a12 + dx*dy ;
			a21 = a21 + dx*dy ;	
			a22 = a22 + dy*dy ;
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

			if ( m->cell[ i ].type != m->cell[ i ].cell[ k ]->type ) {

				Pf[ 0 ] = m->cell[ i ].face[ k ]->x - m->cell[ i ].x ;
				Pf[ 1 ] = m->cell[ i ].face[ k ]->y - m->cell[ i ].y ;
				fPf[ 0 ] = DotProduct( Pf,  m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[0] ;
				fPf[ 1 ] = DotProduct( Pf,  m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[1] ;
				dx = ( -fPf[ 0 ] + m->cell[ i ].face[ k ]->x )  - m->cell[ i ].x ;
				dy = ( -fPf[ 1 ] + m->cell[ i ].face[ k ]->y )  - m->cell[ i ].y ;

				//dx =  cell[ i ].face[ k ]->x - cell[ i ].x ;
				//dy =  cell[ i ].face[ k ]->y - cell[ i ].y ;

			} else {

				dx = m->cell[ i ].cell[ k ]->x  - m->cell[ i ].x ;
				dy = m->cell[ i ].cell[ k ]->y  - m->cell[ i ].y ;
			}
			LSQ_Cx[ k ][ i ] = ia11*dx + ia12*dy ;
			LSQ_Cy[ k ][ i ] = ia21*dx + ia22*dy ;
		}
		//if(mpi_id==0) cout<<"---------------------------------------"<<endl;
		/*--- Cal. LSQ Coefficient for domain boundary "faces" ---*/
		for ( int k = iCell ; k < iFace ; k++ ) {

			Pf[ 0 ] = m->cell[ i ].face[ k ]->x - m->cell[ i ].x ;
			Pf[ 1 ] = m->cell[ i ].face[ k ]->y - m->cell[ i ].y ;
			fPf[ 0 ] = DotProduct( Pf, m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[0] ;
			fPf[ 1 ] = DotProduct( Pf, m->Cell[ i ][ k ].mf )*m->Cell[ i ][ k ].mf[1] ;
			dx = ( -fPf[ 0 ] + m->cell[ i ].face[ k ]->x )  - m->cell[ i ].x ;
			dy = ( -fPf[ 1 ] + m->cell[ i ].face[ k ]->y )  - m->cell[ i ].y ;

			//dx =  cell[ i ].face[ k ]->x - cell[ i ].x ;
			//dy =  cell[ i ].face[ k ]->y - cell[ i ].y ;
			LSQ_Cx[ k ][ i ] = ia11*dx + ia12*dy ;
			LSQ_Cy[ k ][ i ] = ia21*dx + ia22*dy ;
		}//End boundaty face
	}//End cell loop
	// MPI_Barrier(MPI_COMM_WORLD); exit(1) ;
	for( int k = 0 ;  k < 6 ; k++ ){
		LSQ_Cx[ k ] = LSQ_Cx[ k ] ;
		LSQ_Cy[ k ] = LSQ_Cy[ k ] ;
		LSQ_Cz[ k ] = LSQ_Cz[ k ] ;
	}
}

/**
 * @file domain_structure.cpp
 * @author Kuan-Lin Chen (Postdoc., Aerothermal & Plasma Physics Laboratory (APPL), NCTU)
 * @date 15 May 2019
 * @brief In this file, we extract the ultraMPP (ultra-fast Massive Parallel Platform, Plasma T.I.) connectivity.
 * @version ultraMPP 1.9.3
 */
#pragma once
#include "PFM.hpp"
#include "ultraMPP.h"

using namespace std;

/*----------------------------------------------------------

      n2
      o------------o
     .  \         .
    .    \   N   .
   .  P   f     .
  .		     \   .
 .          \ .       
o------------o      P: Present cell
					N: Neighbor cell 

----------------------------------------------------------*/
class CLSQ
{
   public:
   	double Cx[ 7 ], Cy[ 7 ] ;
   	void Init(){
   		for (int i=0 ; i < 7 ; i++ ){
   			Cx[ i ] = 9.E99 ;
   			Cy[ i ] = 9.E99 ;
   		}
   	}
};
class CCell {
	public:
	int NeighborCellId ;
	int NeighborGlobalCellId ;

	double  //PN[2], /*!< \brief Vector of PN. */
		    //Pf[2], /*!< \brief Vector of Pf. */
		    //Nf[2], /*!< \brief Vector of Nf. */
		    nf[3], /*!< \brief Unit normal vector of surface f. */
		    mf[3], /*!< \brief Unit Tangent vector of surface f. */
		    Af[3], /*!< \brief Surface area vector. */
		   PPP[3], /*!< \brief Vector of PP'. */
		   NNP[3], /*!< \brief Vector of NN'. */
		     dPPf, /*!< \brief Distance between P' and f. */
		   	 dNPf, /*!< \brief Distance between N' and f. */
			dDist, /*!< \brief Distance between P' and N'. */
			dArea; /*!< \brief Face area. */
	double SurfaceCharge ;
	void Init()
	{
		NeighborCellId = 0 ;
		NeighborGlobalCellId = 0 ;
		 dPPf = 0.0 ;
		 dNPf = 0.0 ;
		dDist = 0.0 ;
		dArea = 0.0 ;
		for ( int i = 0 ; i < 3 ; i++ ){
			 nf[ i ] = 0.0 ;
			 mf[ i ] = 0.0 ;
			 Af[ i ] = 0.0 ;
			PPP[ i ] = 0.0 ; 
			NNP[ i ] = 0.0 ; 
		}
		//for ( int i = 0 ; i < 7 ; i++ ) 
		SurfaceCharge = 0.0 ; 
	}

};
class CDomain
{
	public:
		CDomain();

		CCell **PFM_CELL ;

		void BulidCellStructure() ;
		
		/*--- Least-Square Granient Reconstruction ---*/
		CLSQ *LSQ ;

		void Calculate_LSQ_Coeff() ;
		
		int local_cell_number ;

		int test ;
		double PN[ 3 ], Pf[ 3 ], Nf[ 3 ], PPf[ 3 ], NPf[ 3 ], fPf[ 3 ], mf[ 3 ] ;
		double LSQ_Weight( double dx, double dy ) {

			//double GradientWeithtPower = 1.0 ; //Inverse distance.
			//double Dis = 0.0 ;
			double w = 0.0 ;

			//Dis = sqrt(dx*dx + dy*dy) ;
       		//w 	= 1.0 / pow(Dis, GradientWeithtPower ) ;
       		w=1.0;
			return w*w ;
		} ;
		double DotProduct(double *A, double *B ){
			//cout<<"A"<<endl;
			return A[0]*B[0] + A[1]*B[1] ;
		};
		void Init() ;
		void UltraMPPExtractFaceCellTag() ;

};
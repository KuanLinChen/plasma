
#include "solver_poisson_explicit.hpp"

using namespace std ;
CEPoisson::CEPoisson()
{
}
void CEPoisson::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config )
{
	Correction = config->Equation[ POISSON ].Correction ;
	if ( mpi_id == 0 ){
		cout<<"Creat POISSON"<<endl;
		cout<<"Correction: "<<Correction<<endl;
	}
	its=0 ;


	Res[ 0 ].initial ( m, CELL, "res. potential" ) ;
	Res[ 1 ].initial ( m, CELL, "res. Ex" ) ;
	Res[ 2 ].initial ( m, CELL, "res. Ey" ) ;
	Source[ 0 ].initial ( m, CELL, "ρ0" ) ;
	Source[ 1 ].initial ( m, CELL, "ρ1" ) ;
	Source[ 2 ].initial ( m, CELL, "ρ2" ) ;
	u_exact.initial ( m, CELL, "u_exact" ) ;
	p_exact.initial ( m, CELL, "p_exact" ) ;

}
void CEPoisson::Solve( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	exact( m, config, var ) ;
	int nCell = m->local_cell_number ;
	double Max_Res_Phi=0.0 ;	

	for ( int iter=0 ; iter < 1000 ; iter++ ){

		Max_Res_Phi=0.0 ;	
		CalculateSourceTerm( m, config, var ) ;
		ComputeResidue( m, config, var ) ;

		for( int i = 0 ; i < nCell ; i++ ) {
			if( m->cell[ i ].type == PLASMA ){
				if( fabs(Res[ 0 ][i]) > Max_Res_Phi ) Max_Res_Phi = fabs(Res[ 0 ][i] ) ;
				//Max_Res_Phi += fabs( (var->dPrePhi[ i ] -var->dPhi[ i ])/var->dPrePhi[ i ] ) ;
			}
		}
		if ( Max_Res_Phi < 1.E-8 ){
			cout<<iter<<endl;
		 break;
		}
		if(iter == 999 )cout<<"Max"<<endl;
	}

	Plot( m, config, var ) ;
}
void CEPoisson::exact( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	int nCell = m->local_cell_number ;
	for ( int i = 0 ; i < nCell ; i++ ) {
    	u_exact[i] = sin( var->PI*m->cell[i].x ) ;
    	p_exact[i] = var->PI*cos(var->PI*m->cell[i].x) ;
	}
}
void CEPoisson::CalculateSourceTerm( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	int nCell = m->local_cell_number ;

	double Lr =1.0/2.0/var->PI ;
	double nu = 1.0 ;
	double Tr = Lr*Lr / nu ;

	Source[0].zero() ;
	Source[1].zero() ;

	for( int i = 0 ; i < nCell ; i++ ) {
		//not complete
		if( m->cell[ i ].type == PLASMA ){

			Source[ 0 ][ i ] = nu*pow( var->PI, 2.0 )*sin( var->PI*m->cell[i].x ) ;
			Source[ 1 ][ i ] = -var->dPreEx[i]/Tr ;

		} else {
		}

	}
}
void CEPoisson::ComputeResidue( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	int j=0, nCell = m->local_cell_number, iFace=0, iCell=0 ;

	double Lr =1.0/2.0/var->PI ;
	double nu = 1.0 ;
	double Tr = Lr*Lr / nu ;
	double flux0 = 0.0, flux1 = 0.0 ;
	double uR=0.0, uL=0.0, pR=0.0, pL=0.0 ;
	double hmin = fabs( m->cell[0].x-m->cell[1].x ) ;
	double dtau = 0.8* hmin/(nu/Lr)  ;
	//cout<<dtau<<endl;exit(1);

	Res[0].zero() ;
	Res[1].zero() ;

	for( int i = 0 ; i < nCell ; i++ ) {

		iFace 	 = m->cell[ i ].face_number ;
		iCell 	 = m->cell[ i ].cell_number ;

		for ( int k = 0 ; k < iCell ; k++ ) {

			j 	  = m->Cell[ i ][ k ].NeighborCellId ;


	   //      if ( m->Cell[ i ][ k ].nf[ 0 ] > 0.0 ){

		  //       uL = var->dPrePhi[ i ] ;
		  //       pL = var-> dPreEx[ i ] ;

		  //       uR = var->dPrePhi[ j ] ;
		  //       pR = var-> dPreEx[ j ] ;

		  //      	flux0 = 0.5*( ( nu*pR + nu*pL )    + nu/Lr*( uR - uL ) ) ;
    //     		flux1 = 0.5*( (    uR +    uL )/Tr + nu/Lr*( pR - pL ) ) ;
	   //      	Res[ 0 ][ i ] += flux0*m->Cell[ i ][ k ].dArea* m->Cell[ i ][ k ].nf[ 0 ] ;
				// Res[ 1 ][ i ] += flux1*m->Cell[ i ][ k ].dArea* m->Cell[ i ][ k ].nf[ 0 ] ;

	   //      } else if ( m->Cell[ i ][ k ].nf[ 0 ] < 0.0){

		  //       uL = var->dPrePhi[ j ] ;
		  //       pL = var-> dPreEx[ j ] ;

		  //       uR = var->dPrePhi[ i ] ;
		  //       pR = var-> dPreEx[ i ] ;

		  //      	flux0 = 0.5*( ( nu*pR + nu*pL )    + nu/Lr*( uR - uL ) ) ;
    //     		flux1 = 0.5*( (    uR +    uL )/Tr + nu/Lr*( pR - pL ) ) ;
	   //      	Res[ 0 ][ i ] += flux0*m->Cell[ i ][ k ].dArea* m->Cell[ i ][ k ].nf[ 0 ] ;
				// Res[ 1 ][ i ] += flux1*m->Cell[ i ][ k ].dArea* m->Cell[ i ][ k ].nf[ 0 ] ;
	   //      }



			uL = var->dPrePhi[ i ] ;
			pL = var-> dPreEx[ i ] ;

			uR = var->dPrePhi[ j ] ;
			pR = var-> dPreEx[ j ] ;

			flux0 = 0.5*( ( nu*pR + nu*pL )    )* m->Cell[ i ][ k ].nf[ 0 ] + 0.5*(nu/Lr*( uR - uL ) ) ;
			flux1 = 0.5*( (    uR +    uL )/Tr )* m->Cell[ i ][ k ].nf[ 0 ] + 0.5*(nu/Lr*( pR - pL ) ) ;

			Res[ 0 ][ i ] += flux0*m->Cell[ i ][ k ].dArea ;
			Res[ 1 ][ i ] += flux1*m->Cell[ i ][ k ].dArea ;
				
		}//Loop over neighbor cells

		/*--------------------------------------------------------------*/
		for( int k = iCell ; k < iFace ; k++ ) {

			if ( m->cell[ i ].face[ k ]->type == NEUMANN ) {

				uL = var->dPrePhi[ i ] ;
		        pL = var-> dPreEx[ i ] ;

				flux0 = ( nu*pL )    ;
	        	flux1 = (    uL )/Tr;

	        	Res[ 0 ][ i ] += flux0*m->Cell[ i ][ k ].dArea * m->Cell[ i ][ k ].nf[0] ;
				Res[ 1 ][ i ] += flux1*m->Cell[ i ][ k ].dArea * m->Cell[ i ][ k ].nf[0] ;

			}else{

	        uL = var->dPrePhi[ i ] ;
	        pL = var-> dPreEx[ i ] ;

	        uR = 0 ;
	        pR = var-> dPreEx[ i ] ;//neumann

			//flux0 = 0.5*( ( nu*pR + nu*pL )    + nu/Lr*( uR - uL ) )* m->Cell[ i ][ k ].nf[ 0 ] ;
        	//flux1 = 0.5*( (    uR +    uL )/Tr + nu/Lr*( pR - pL ) )* m->Cell[ i ][ k ].nf[ 0 ] ;

			flux0 = ( nu*pL )* m->Cell[ i ][ k ].nf[ 0 ] ;
        	flux1 = 	0.0  * m->Cell[ i ][ k ].nf[ 0 ];

        	Res[ 0 ][ i ] += flux0*m->Cell[ i ][ k ].dArea ;//* m->Cell[ i ][ k ].nf[0] ;
			Res[ 1 ][ i ] += flux1*m->Cell[ i ][ k ].dArea ;//* m->Cell[ i ][ k ].nf[0] ;
			}

		}

		Res[ 0 ][ i ] += Source[ 0 ][ i ]*m->cell[i].volume ;
		Res[ 1 ][ i ] += Source[ 1 ][ i ]*m->cell[i].volume ;

		var->dPhi[ i ] = var->dPrePhi[ i ] + (dtau/m->cell[i].volume)*Res[ 0 ][ i ] ;
        var-> dEx[ i ] = var-> dPreEx[ i ] + (dtau/m->cell[i].volume)*Res[ 1 ][ i ] ;
	}//Cell Loop

	for( int i = 0 ; i < nCell ; i++ ) {
		var->dPrePhi[ i ] = var->dPhi[ i ] ;
		var-> dPreEx[ i ] = var-> dEx[ i ] ;
	}
}
void CEPoisson::Plot( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	IO Output ;
	Output.set_output_type ( config->OutputFormat ) ;
	Output.set_filename( "Poisson.dat" ) ;
	
	Output <<  u_exact ;
	Output <<  p_exact ;
	Output <<  var->dPhi ;
	Output <<  var->dEx ;
	Output <<  Res[0] ;
	Output <<  Res[1] ;
	Output <<  Source[0] ;
	Output <<  Source[1] ;


	Output<<endl;
}
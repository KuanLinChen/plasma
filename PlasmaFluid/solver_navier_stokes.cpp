
#include "solver_navier_stokes.hpp"

using namespace std ;
CNavierStokes::CNavierStokes()
{
}
void CNavierStokes::Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, int index )
{
	iSpecies = index ;
	SpeciesType = config->Species[ index ].Type ;

	/*--- PETSc Solver ---*/	
		Correction = config->Equation[ SpeciesType ].Correction ;
		if ( mpi_id == 0 ){
			cout<<"Creat "<<config->Species[iSpecies].Name<<" navier stokes solver, index: "<<index<<", charge: "<<config->Species[ index ].Charge<<", Speciec Type: "<<SpeciesType<<endl ;
			cout<<"Correction: "<<Correction<<endl;
			cout<<"Gamma: "<<config->Species[ iSpecies ].Gamma<<endl;
		} 
		if ( m->cylindrical_y ) {
			Omaga = 1.0 ;
		} else {
			Omaga = 0.0 ;
		}
		Mass = config->Species[iSpecies].Mass_Kg ;

	/*--- Allocate Residue ---*/
		Res = new double*[ 5 ] ;/* N, Ux, Uy, Uz, E */
		for ( int iEqn = 0 ; iEqn < 5 ; iEqn ++ ) {
			Res[iEqn] = new double [ m->local_cell_number ] ;
		}
		buffer = new double*[ 5 ] ;
		for ( int iEqn = 0 ; iEqn < 5 ; iEqn ++ ) buffer[iEqn] = new double [ mpi_size ] ;
}
void CNavierStokes::Solve( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &variable, boost::shared_ptr<CPost> &post  )
{
	cout.precision(5);
	for ( int iter = 0 ; iter < 10000000  ; iter ++ ) {	

		variable->UpdateSolution( m, config ) ;
		CalculateLocalTimeStep( m, config, variable ) ;
		ComputeFlux_HLL( m, config, variable ) ;
		ContinuityIntegral( m, config, variable ) ;
		MomentumIntegral( m, config, variable ) ;
		EnergyIntegral( m, config, variable ) ;

		for (int i=0 ; i < 5 ; i++ ){
			MPI_Allgather( &MaxERR[i], 1, MPI_DOUBLE, buffer[i], 1, MPI_DOUBLE, MPI_COMM_WORLD ) ;
		}

		for ( int iEqn = 0 ; iEqn < 5 ; iEqn ++ ) {
			for ( int i = 0 ; i < mpi_size ; i++ ) {
				if ( buffer[iEqn][ i ] > MaxERR[iEqn] ) MaxERR[iEqn] = buffer[iEqn][ i ] ;
			}
		}
		if( mpi_id==MASTER_NODE and iter %100 == 0 ) {
			for ( int iEqn = 0 ; iEqn < 5 ; iEqn ++ ){
				cout<< scientific <<MaxERR[iEqn]<<setw(22);
			}cout<<endl;
		}
		if( iter %2000 == 0 ){
			post->OutputFlow_NS( m, config, variable ) ;
			if( mpi_id==MASTER_NODE ) cout<<"PLOT"<<endl;
		} 
		if( MaxERR[0] < 1.E-12 ){
			cout<<MaxERR[0]<<endl;
			break ;
		} 
	}
}
void CNavierStokes::ComputeFlux_HLL( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var )
{
	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;

	double Rho=0.0, RhoL=0.0, RhoR=0.0 ;
	double u  =0.0, uL  =0.0, uR  =0.0 ; double v=0.0, vL=0.0, vR=0.0 ; double w=0.0, wL=0.0, wR=0.0 ; double  un=0.0, unL=0.0,  unR=0.0 ;
	double p  =0.0, pL  =0.0, pR  =0.0 ; double a=0.0, aL=0.0, aR=0.0 ; double H=0.0, HL=0.0, HR=0.0 ;

	double F0 =0.0, F1 =0.0, F2 =0.0, F3 =0.0, F4 =0.0 ; /* iFlux @ interface */
	double F0L=0.0, F1L=0.0, F2L=0.0, F3L=0.0, F4L=0.0 ; /* iFlux @ left cell */
	double F0R=0.0, F1R=0.0, F2R=0.0, F3R=0.0, F4R=0.0 ; /* iFlux @ Reft cell */
	double  SL=0.0, SR=0.0, SLm=0.0, SRp=0.0 ;
	double LogicalSwitch = 0.0, En = 0.0, GradUvel=0.0, GradU4=0.0, uma=0.0, upa=0.0 ;
	double RT=0.0, RiemannPlus=0.0, RiemannMinus=0.0, Vn=0.0, Entropy=0.0 ;
	/*------------------*/
	double Rho_inf 	= 1.0, 
		   M_inf 	= 1.65, 
		   u_inf 	= M_inf, 
		   v_inf = 0.0, P_inf = 1.0/config->Species[ iSpecies ].Gamma, 
		   un_inf = 0.0,
		   a_inf = sqrt(config->Species[ iSpecies ].Gamma*P_inf/Rho_inf), 
		   H_inf = a_inf*a_inf/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*( u_inf*u_inf + v_inf*v_inf), 
		   Entropy_Inf = pow(Rho_inf,config->Species[ iSpecies ].Gamma)/P_inf ;

	for( int i = 0 ; i < nCell ; i++ ) {

		iFace 	 = m->cell[ i ].face_number ;
		iCell 	 = m->cell[ i ].cell_number ;

		for( int k = 0 ; k < 5 ; k++ ) Res[ k ][ i ] = 0.0 ;

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*m->cell[ i ].x ; 

		/*--- Loop over PLASMA cells ---*/
		if ( m->cell[ i ].type == PLASMA ){

			/*--- Loop over bulk faces ---*/
			for ( int k = 0 ; k < iCell ; k++ ){

				j = m->Cell[ i ][ k ].NeighborCellId ;

				if ( m->cell[ j ].type == PLASMA ){

					RhoL = var->PreU0[ i ] ;
					RhoR = var->PreU0[ j ] ;

					uL 	 = var->PreU1[ i ]/RhoL ;
					uR 	 = var->PreU1[ j ]/RhoR ;

					vL 	 = var->PreU2[ i ]/RhoL ;
					vR 	 = var->PreU2[ j ]/RhoR ;

					wL 	 = var->PreU3[ i ]/RhoL ;
					wR 	 = var->PreU3[ j ]/RhoR ;

					unL  = uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] ;// + wL*m->Cell[ i ][ k ].nf[ 2 ] ;
					unR  = uR*m->Cell[ i ][ k ].nf[ 0 ] + vR*m->Cell[ i ][ k ].nf[ 1 ] ;// + wR*m->Cell[ i ][ k ].nf[ 2 ] ;

					pL 	= (config->Species[ iSpecies ].Gamma-1.0)*(  var->PreU4[ i ] - 0.5*RhoL*(uL*uL + vL*vL) ) ;
					pR  = (config->Species[ iSpecies ].Gamma-1.0)*(  var->PreU4[ j ] - 0.5*RhoR*(uR*uR + vR*vR) ) ;

					aL  = sqrt(config->Species[ iSpecies ].Gamma*pL/RhoL) ;
					aR  = sqrt(config->Species[ iSpecies ].Gamma*pR/RhoR) ;

					HL  = aL*aL/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*(uL*uL+vL*vL) ;
					HR  = aR*aR/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*(uR*uR+vR*vR) ;


					/*---  Davis' wave estimate---*/
					RT 	= sqrt( RhoR/RhoL ) ;
					Rho = RT*RhoL ;
    				un 	= ( unL+RT*unR)/( 1.0 + RT ) ;
     				H 	= ( HL+RT*HR  )/( 1.0 + RT ) ;
     				a = sqrt( (config->Species[ iSpecies ].Gamma-1.0)*( H - 0.5*un*un) ) ;
     				
					uma = un - a ;
					upa = un + a ;
					SL 	= min( unL - aL, uma ) ;
					SR 	= max( unR + aR, upa ) ;

					// SL 	= min( unL - aL, unR - aR ) ;
					// SR 	= max( unL + aL, unR + aR ) ;
					SLm = min( SL, 0.0 ) ;
					SRp = max( SR, 0.0 ) ;


					/*---  Left Physical Flux  ---*/
					F0L = RhoL*unL ;
					F1L = RhoL*unL * uL + pL*m->Cell[ i ][ k ].nf[ 0 ] ;
					F2L = RhoL*unL * vL + pL*m->Cell[ i ][ k ].nf[ 1 ] ;
					F3L = RhoL*unL * wL + pL*m->Cell[ i ][ k ].nf[ 2 ] ;
					F4L = RhoL*unL * HL ;


					/*---  Right Physical Flux  ---*/
					F0R = RhoR*unR ;
					F1R = RhoR*unR * uR + pR*m->Cell[ i ][ k ].nf[ 0 ] ;
					F2R = RhoR*unR * vR + pR*m->Cell[ i ][ k ].nf[ 1 ] ;
					F3R = RhoR*unR * wR + pR*m->Cell[ i ][ k ].nf[ 2 ] ;
					F4R = RhoR*unR * HR ;

					/*---  Numerical Flux, Note: outward flux  ---*/
				    F0 = ( SRp*F0L - SLm*F0R + (SLm*SRp)*(var->PreU0[ j ]-var->PreU0[ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F1 = ( SRp*F1L - SLm*F1R + (SLm*SRp)*(var->PreU1[ j ]-var->PreU1[ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F2 = ( SRp*F2L - SLm*F2R + (SLm*SRp)*(var->PreU2[ j ]-var->PreU2[ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F3 = ( SRp*F3L - SLm*F3R + (SLm*SRp)*(var->PreU3[ j ]-var->PreU3[ i ]) )/(SRp-SLm);//+1.e-19) ;
				    F4 = ( SRp*F4L - SLm*F4R + (SLm*SRp)*(var->PreU4[ j ]-var->PreU4[ i ]) )/(SRp-SLm);//+1.e-19) ;

				    /*---  Adding to residue  ---*/
				    Res[ 0 ][ i ] += LogicalSwitch*F0*m->Cell[ i ][ k ].dArea ;
				    Res[ 1 ][ i ] += LogicalSwitch*F1*m->Cell[ i ][ k ].dArea ;
				    Res[ 2 ][ i ] += LogicalSwitch*F2*m->Cell[ i ][ k ].dArea ;
				    Res[ 3 ][ i ] += LogicalSwitch*F3*m->Cell[ i ][ k ].dArea ;
				    Res[ 4 ][ i ] += LogicalSwitch*F4*m->Cell[ i ][ k ].dArea ;

	 			} else {/*--- For discontuity face ---*/

	 				RhoL = var->PreU0[ i ] ;
					uL 	 = var->PreU1[ i ]/RhoL ;
					vL 	 = var->PreU2[ i ]/RhoL ;
					wL 	 = var->PreU3[ i ]/RhoL ;
					unL  = uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] ;//+ wL*m->Cell[ i ][ k ].nf[ 2 ] ;
					pL 	= (config->Species[ iSpecies ].Gamma-1.0)*(  var->PreU4[ i ] - 0.5*RhoL*(uL*uL + vL*vL) ) ;

					F0 = 0.0 ;
					F1 = pL*m->Cell[ i ][ k ].nf[ 0 ] ;
					F2 = pL*m->Cell[ i ][ k ].nf[ 1 ] ;
					F3 = pL*m->Cell[ i ][ k ].nf[ 2 ] ;
					F4 = 0.0 ;

				    Res[ 0 ][ i ] += LogicalSwitch*F0*m->Cell[ i ][ k ].dArea ;
				    Res[ 1 ][ i ] += LogicalSwitch*F1*m->Cell[ i ][ k ].dArea ;
				    Res[ 2 ][ i ] += LogicalSwitch*F2*m->Cell[ i ][ k ].dArea ;
				    Res[ 3 ][ i ] += LogicalSwitch*F3*m->Cell[ i ][ k ].dArea ;
				    Res[ 4 ][ i ] += LogicalSwitch*F4*m->Cell[ i ][ k ].dArea ;
	 			}
	 		}//End bulk face

	 		
			/*--- Loop over boundary faces ---*/
	 		for( int k = iCell ; k < iFace ; k++ ) {

	 			if ( m->cell[ i ].face[ k ]->type == INLET or m->cell[ i ].face[ k ]->type == OUTLET ) {

	 				RhoL= var->PreU0[ i ] ;
					uL 	= var->PreU1[ i ]/RhoL ;
					vL 	= var->PreU2[ i ]/RhoL ;
					un 	= uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] ;
					pL 	= (config->Species[ iSpecies ].Gamma-1.0)*(  var->PreU4[ j ] - 0.5*Rho*(uL*uL + vL*vL) ) ;
					aL 	= sqrt(config->Species[ iSpecies ].Gamma*pL/RhoL) ;


	 				un_inf = u_inf*m->Cell[ i ][ k ].nf[ 0 ] + v_inf*m->Cell[ i ][ k ].nf[ 1 ] ;

					if ( un_inf > -a_inf) {
				        /*--- Subsonic inflow or outflow ---*/
				    	RiemannPlus = un 	 + 2.0 * aL 	/(config->Species[ iSpecies ].Gamma-1.0);
					} else {
				        /*--- Supersonic inflow ---*/
				        RiemannPlus = un_inf + 2.0 * a_inf	/(config->Species[ iSpecies ].Gamma-1.0);
					}

					/*--- Check whether (u.n-c) is greater or less than zero ---*/
					if (un_inf > a_inf) {
						/*--- Supersonic outflow ---*/
						RiemannMinus = un - 2.0*aL/(config->Species[ iSpecies ].Gamma-1.0);
					} else {
						/*--- Subsonic outflow ---*/
						RiemannMinus = un_inf - 2.0*a_inf/(config->Species[ iSpecies ].Gamma-1.0);
					}

				      /*--- Compute a new value for the local normal velocity and speed of
				         sound from the Riemann invariants. ---*/

				      Vn = 0.5 * (RiemannPlus + RiemannMinus);
				      a = 0.25 * (RiemannPlus - RiemannMinus)*(config->Species[ iSpecies ].Gamma-1.0);

				      /*--- Construct the primitive variable state at the boundary for
				         computing the flux for the weak boundary condition. The values
				         that we choose to construct the solution (boundary or freestream)
				         depend on whether we are at an inflow or outflow. At an outflow, we
				         choose boundary information (at most one characteristic is incoming),
				         while at an inflow, we choose infinity values (at most one
				         characteristic is outgoing). ---*/

				      if ( un_inf > 0.0)   {

				        /*--- Outflow conditions ---*/
				        u = uL + (Vn - un )*m->Cell[ i ][ k ].nf[ 0 ] ;
				        v = vL + (Vn - un )*m->Cell[ i ][ k ].nf[ 1 ] ;
				        Entropy= pow(RhoL,config->Species[ iSpecies ].Gamma)/pL ;

				      } else  {

				        /*--- Inflow conditions ---*/
				        u = u_inf + (Vn - un_inf )*m->Cell[ i ][ k ].nf[ 0 ] ;
				        v = v_inf + (Vn - un_inf )*m->Cell[ i ][ k ].nf[ 1 ] ;
				        Entropy = Entropy_Inf;

				      }

				      /*--- Recompute the primitive variables. ---*/
				      RhoL = pow( Entropy*a*a/config->Species[ iSpecies ].Gamma, 1.0/(config->Species[ iSpecies ].Gamma-1.0));
				      pL = RhoL*a*a/config->Species[ iSpecies ].Gamma ;

				      HL   = a*a/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*(u*u+v*v) ;
				      unL = u*m->Cell[ i ][ k ].nf[ 0 ] + v*m->Cell[ i ][ k ].nf[ 1 ] ;

					// RhoL = Rho_inf ;
					// uL   = u_inf ;
			  		//vL   = 0.0 ;
			  		//wL	 = 0.0 ;
			  		//pL   = P_inf ;
			  		//unL  = uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] ;//+ wL*m->Cell[ i ][ k ].nf[ 2 ] ;
			  		//aL   = sqrt(config->Species[ iSpecies ].Gamma*pL/RhoL) ;
			  		//HL   = aL*aL/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*(uL*uL+vL*vL) ;

					F0 = RhoL*unL ;
					F1 = RhoL*unL * uL + pL*m->Cell[ i ][ k ].nf[ 0 ] ;
					F2 = RhoL*unL * vL + pL*m->Cell[ i ][ k ].nf[ 1 ] ;
					F3 = RhoL*unL * wL + pL*m->Cell[ i ][ k ].nf[ 2 ] ;
					F4 = RhoL*unL * HL ;

    				Res[ 0 ][ i ] += LogicalSwitch*F0*m->Cell[ i ][ k ].dArea ;
				    Res[ 1 ][ i ] += LogicalSwitch*F1*m->Cell[ i ][ k ].dArea ;
				    Res[ 2 ][ i ] += LogicalSwitch*F2*m->Cell[ i ][ k ].dArea ;
				    Res[ 3 ][ i ] += LogicalSwitch*F3*m->Cell[ i ][ k ].dArea ;
				    Res[ 4 ][ i ] += LogicalSwitch*F4*m->Cell[ i ][ k ].dArea ;

	 			// }else if( m->cell[ i ].face[ k ]->type == OUTLET ){

					// /*--- Left state ---*/
	 			// 	RhoL = var->PreU0[ i ] ;
					// uL 	 = var->PreU1[ i ]/RhoL ;
					// vL 	 = var->PreU2[ i ]/RhoL ;
					// wL 	 = var->PreU3[ i ]/RhoL ;

					// unL  = uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] ;//+ wL*m->Cell[ i ][ k ].nf[ 2 ] ;

					// pL 	= ( config->Species[ iSpecies ].Gamma-1.0 )*(  var->PreU4[ i ] - 0.5*RhoL*(uL*uL + vL*vL) ) ;

					// aL  = sqrt(config->Species[ iSpecies ].Gamma*pL/RhoL) ;

					// HL  = aL*aL/(config->Species[ iSpecies ].Gamma-1.0) + 0.5*(uL*uL+vL*vL) ;
					

					// F0 = RhoL*unL ;
					// F1 = RhoL*unL * uL + pL*m->Cell[ i ][ k ].nf[ 0 ] ;
					// F2 = RhoL*unL * vL + pL*m->Cell[ i ][ k ].nf[ 1 ] ;
					// F3 = RhoL*unL * wL + pL*m->Cell[ i ][ k ].nf[ 2 ] ;
					// F4 = RhoL*unL * HL ;

					// /*---  Adding to residue  ---*/
				 // 	Res[ 0 ][ i ] += LogicalSwitch*F0*m->Cell[ i ][ k ].dArea ;
				 // 	Res[ 1 ][ i ] += LogicalSwitch*F1*m->Cell[ i ][ k ].dArea ;
				 // 	Res[ 2 ][ i ] += LogicalSwitch*F2*m->Cell[ i ][ k ].dArea ;
				 // 	Res[ 3 ][ i ] += LogicalSwitch*F3*m->Cell[ i ][ k ].dArea ;
				 // 	Res[ 4 ][ i ] += LogicalSwitch*F4*m->Cell[ i ][ k ].dArea ;

	 			} else if ( m->cell[ i ].face[ k ]->type == EULER_WALL ) {

	 				RhoL = var->PreU0[ i ] ;
					uL 	 = var->PreU1[ i ]/RhoL ;
					vL 	 = var->PreU2[ i ]/RhoL ;
					wL 	 = var->PreU3[ i ]/RhoL ;
					unL  = uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] ;// + wL*m->Cell[ i ][ k ].nf[ 2 ] ;
					pL 	= (config->Species[ iSpecies ].Gamma-1.0)*(  var->PreU4[ i ] - 0.5*RhoL*(uL*uL + vL*vL) ) ;

					F0 = 0.0 ;
					F1 = pL*m->Cell[ i ][ k ].nf[ 0 ] ;
					F2 = pL*m->Cell[ i ][ k ].nf[ 1 ] ;
					F3 = pL*m->Cell[ i ][ k ].nf[ 2 ] ;
					F4 = 0.0 ;

				    Res[ 0 ][ i ] += LogicalSwitch*F0*m->Cell[ i ][ k ].dArea ;
				    Res[ 1 ][ i ] += LogicalSwitch*F1*m->Cell[ i ][ k ].dArea ;
				    Res[ 2 ][ i ] += LogicalSwitch*F2*m->Cell[ i ][ k ].dArea ;
				    Res[ 3 ][ i ] += LogicalSwitch*F3*m->Cell[ i ][ k ].dArea ;
				    Res[ 4 ][ i ] += LogicalSwitch*F4*m->Cell[ i ][ k ].dArea ;
	 			}
	 		}

	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		Res[ 0 ][ i ] = 0.0 ;
			Res[ 1 ][ i ] = 0.0 ;
			Res[ 2 ][ i ] = 0.0 ;
			Res[ 3 ][ i ] = 0.0 ;
			Res[ 4 ][ i ] = 0.0 ;
	 	}
	 	//cout<<"R2: "<<Res[ 2 ][ i ]<<endl;
	}//Cell Loop
	//exit(1) ;
}
void CNavierStokes::ContinuityIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var  )
{
	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double LogicalSwitch = 0.0, SourceSink=0.0, Source=0.0 ;
	MaxERR[ 0 ] = 0.0 ;
	for( int i = 0 ; i < nCell ; i++ ) {

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*m->cell[ i ].x ; 

		/*--- Loop over PLASMA cells ---*/
		if ( m->cell[ i ].type == PLASMA ){

	 		/*--- Source/Sink term ---*/
	 		//SourceSink = (double)*( var->ReactionRatePoint + i  ) ;
	 		//Source = (-1.0)*SourceSink*m->cell[ i ].volume ;
	 		//Res[ 0 ][ i ] += LogicalSwitch*Source ;
	 		Res[ 0 ][ i ]  = Res[ 0 ][ i ]/m->cell[ i ].volume/LogicalSwitch ;
	 		//cout<<LogicalSwitch<<endl;
	 		if ( fabs( Res[ 0 ][ i ] ) > MaxERR[0] ) MaxERR[0] = fabs( Res[ 0 ][ i ] ) ;
			/*--- Integrate ---*/
			var->U0[ i ]  = var->PreU0[ i ] - var->Dt[i] * Res[ 0 ][ i ] ;

	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		var->U0[ i ] = 0.0 ;
	 	}
	}//Cell Loop
// 	var->U0 = var->U0 ;
}
void CNavierStokes::MomentumIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var  )
{

	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double U=0.0, V=0.0, W=0.0, VTOT2=0.0, ReducedMass=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0 ;
	double LogicalSwitch=0.0 ;
	double RhoL=0.0, uL=0.0, vL=0.0, wL=0.0, pL=0.0 ;
	MaxERR[1]=0.0;
	MaxERR[2]=0.0;
	MaxERR[3]=0.0;
	for( int i = 0 ; i < nCell ; i++ ) {

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*m->cell[ i ].x ; 

		/*--- Loop over PLASMA cells ---*/
		if ( m->cell[ i ].type == PLASMA ){

			/*--- Y-Axis-symmetric ---*/
			RhoL = var->PreU0[ i ] ;
			uL 	 = var->PreU1[ i ]/RhoL ;
			vL 	 = var->PreU2[ i ]/RhoL ;
			wL 	 = var->PreU3[ i ]/RhoL ;
			//unL  = uL*m->Cell[ i ][ k ].nf[ 0 ] + vL*m->Cell[ i ][ k ].nf[ 1 ] + wL*m->Cell[ i ][ k ].nf[ 2 ] ;
			pL 	= (config->Species[ iSpecies ].Gamma-1.0)*(  var->U4[ i ] - 0.5*RhoL*(uL*uL + vL*vL) ) ;

	        Res[ 1 ][ i ] += (-1.0) * Omaga * pL * m->cell[ i ].volume ;

	 		/*--- Integrate ---*/
	 		Res[ 1 ][ i ]  = Res[ 1 ][ i ]/m->cell[ i ].volume/LogicalSwitch ;
	 		Res[ 2 ][ i ]  = Res[ 2 ][ i ]/m->cell[ i ].volume/LogicalSwitch ;
	 		Res[ 3 ][ i ]  = Res[ 3 ][ i ]/m->cell[ i ].volume/LogicalSwitch ;
	 		if ( fabs( Res[ 1 ][ i ] ) > MaxERR[ 1 ] ) MaxERR[ 1 ] = fabs( Res[ 1 ][ i ] ) ;
	 		if ( fabs( Res[ 2 ][ i ] ) > MaxERR[ 2 ] ) MaxERR[ 2 ] = fabs( Res[ 2 ][ i ] ) ;
	 		if ( fabs( Res[ 3 ][ i ] ) > MaxERR[ 3 ] ) MaxERR[ 3 ] = fabs( Res[ 3 ][ i ] ) ;
			var->U1[ i ]  = var->PreU1[ i ] - var->Dt[i]* Res[ 1 ][ i ] ;
			var->U2[ i ]  = var->PreU2[ i ] - var->Dt[i]* Res[ 2 ][ i ] ;
			var->U3[ i ]  = var->PreU3[ i ] - var->Dt[i]* Res[ 3 ][ i ] ;


	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		var->U1[ i ] = 0.0 ;
	 		var->U2[ i ] = 0.0 ;
	 		var->U3[ i ] = 0.0 ;
	 	}
	}//Cell Loop

	// var->U1 = var->U1 ;
	// var->U2 = var->U2 ;
	// var->U3 = var->U3 ;
}
void CNavierStokes::EnergyIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var  )
{

	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double U=0.0, V=0.0, W=0.0, VTOT2=0.0, ReducedMass=0.0, JHeating=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0, F4=0.0, CapitalE=0.0 ;
	double RhoR=0.0, RhoL=0.0, uL=0.0, uR=0.0, vL=0.0, vR=0.0, wL=0.0, wR=0.0, unL=0.0, unR=0.0, pL=0.0, pR=0.0, aL=0.0, aR=0.0, TiL=0.0, TiR=0.0, TeL=0.0, TeR=0.0, eL=0.0, eR=0.0 ;
	MaxERR[4] = 0.0 ;
	double LogicalSwitch=0.0, Ad_dPN=0.0 ;

	for( int i = 0 ; i < nCell ; i++ ) {

		LogicalSwitch = ( 1.0 - Omaga ) + Omaga*m->cell[ i ].x ; 

		iCell 	 = m->cell[ i ].cell_number ;

		/*--- Loop over PLASMA cells ---*/
		if ( m->cell[ i ].type == PLASMA ){

	 		/*--- Integrate ---*/
	 		Res[ 4 ][ i ]  = Res[ 4 ][ i ]/m->cell[ i ].volume/LogicalSwitch ;
			var->U4[ i ]  = var->PreU4[ i ] - var->Dt[i] * Res[ 4 ][ i ] ;

	 		if ( fabs( Res[ 4 ][ i ] ) > MaxERR[ 4 ] ) MaxERR[ 4 ] = fabs( Res[ 4 ][ i ] ) ;
	 	/*--- Loop over SOLID cells ---*/
	 	} else {
	 		var->U4[ i ] = 0.0 ;
	 	}
	}//Cell Loop
	//exit(1) ;
	//var->U4 = var->U4 ;
}
void CNavierStokes::CalculateLocalTimeStep( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var  )
{

	int iFace=0, iCell=0, j=0 ;
	int nCell = m->local_cell_number ;
	double U=0.0, V=0.0, W=0.0, P=0.0, VTOT2=0.0 ;
	double U0=0.0, U1=0.0, U2=0.0, U3=0.0, U4=0.0, F4=0.0, CapitalE=0.0 ;
	double Rho=0.0, Sound=0.0, Max_Dt=0.0, CFL=0.1, Length=0.0 ;

	for ( int i = 0 ; i < nCell ; i++ ) {

		Rho = var->PreU0[ i ] ;
		U = var->PreU1[ i ]/var->PreU0[ i ] ;
		V = var->PreU2[ i ]/var->PreU0[ i ] ;
		VTOT2 = sqrt( pow(U, 2.0 ) + pow(V, 2.0 ) ) ;

		Length = sqrt( m->cell[ i ].volume ) ;

		P = (config->Species[ iSpecies ].Gamma-1.0)*(  var->PreU4[ i ] - 0.5*Rho*(U*U + V*V) ) ;

		Sound = sqrt( config->Species[ iSpecies ].Gamma*P/Rho ) ;

		var->Dt[ i ] = 	CFL*Length/( fabs(VTOT2)+Sound )  ;
		//SoundSpeed[ i ] = Sound ;
		var->MachNumber[ i ] = VTOT2/Sound ;
		//if ( var->Dt[ i ] > Max_Dt ) Max_Dt = var->Dt[ i ] ;
	}
	var->Dt = var->Dt ;
}

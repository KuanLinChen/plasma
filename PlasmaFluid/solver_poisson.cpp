
#include "solver_poisson.hpp"
#include "petscsys.h"
using namespace std ;
CPoisson::CPoisson()
{
}
void CPoisson::Init( boost::shared_ptr<CConfig> &config )
{
	if ( mpi_rank == 0 ){
		cout<<"Creat POISSON"<<endl;
		cout<<"Correction: "<<config->Equation[ POISSON ].Correction<<endl;
	}
	its=0 ;
}
void CPoisson::SOLVE( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	UpdateElectricalMap( config, var ) ;
	/*--- compute the transport coefficient (permittivity). ---*/
	switch ( config->Equation[ POISSON ].Equation ) {
		case 0: 	UltraMPPCOmputePermitt                ( config, var ) ; break ;
		case 1:		UltraMPPCOmputeEffectivePermitt       ( config, var ) ; break ;
		case 2:		UltraMPPCOmputeEffectivePermittEleOnly( config, var ) ; break ;
		default: 
		if( mpi_rank == 0 ) cout << "No such equation option, Pls contact K.-L. Chen " << endl ; exit(1) ; break ;
	}

	/*--- Calculate the net charge density for poisson's source term ---*/
	UltraMPPComputeNetCharge( config, var ) ;

	/*--- set the cell paramtert ---*/
	plasma.set_cell_property_parameter( var->eps_eff ) ;

	/*--- Matrix A ---*/
	plasma.before_matrix_construction() ;
	plasma.add_laplacian_matrix_form_op();
	plasma.finish_matrix_construction() ;

	/* face value assiang */
	double voltage= SineVoltage( "POWER", config, var ) ;

	//cout<<"voltage : "<<voltage<<endl;

	plasma.set_bc_value( MPP_face_tag["POWER" ], voltage, var->Potential.face ) ;
	plasma.set_bc_value( MPP_face_tag["GROUND"],     0.0, var->Potential.face ) ;
	plasma.set_bc_value( MPP_face_tag["NEUMANN"],     0.0, var->Potential.face ) ;

	UltraMPPCOmputeSurfaceCharge( config, var ) ;

	plasma.syn_parallel_cell_data( var->Potential.tag_previous, var->Potential.tag_current);

	/* Source B */
	plasma.before_source_term_construction();
	plasma.add_laplacian_source_term_op(var->ChargeDen, var->Potential.face ) ;
	plasma.finish_source_term_construction();

	/* SOLVE */
	plasma.get_solution( var->Potential.current );

	its = plasma.get_iteration_number() ;

	plasma.get_gradient(&var->Potential);

	/* inverse gradient for Electric field*/
	for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) 
	{
		var->Ex[ i ] = - var->Potential.gradient[0][i] ;
		var->Ey[ i ] = - var->Potential.gradient[1][i] ;
	}
	plasma.syn_parallel_cell_data( var->VarTag["Ex"] );
	plasma.syn_parallel_cell_data( var->VarTag["Ey"] );

	if( plasma.Mesh.ndim == 3 ) {	
		for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {
			var->Ez[ i ] = - var->Potential.gradient[2][i] ;
		}
		plasma.syn_parallel_cell_data( var->VarTag["Ez"] );
	} 

	//Tmp code.
	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ){
		var->EField[ 0 ][ i ] = var->Ex[ i ] ; 
		var->EField[ 1 ][ i ] = var->Ey[ i ] ; 
	} 
	var->EField[0] = var->EField[0] ; 
	var->EField[1] = var->EField[1] ;

	UltraMPPComputeDispCurrentDensity( var ) ;

}
void CPoisson::UltraMPPComputeNetCharge( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*Only the 'plasma' region. */
  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;
    var->ChargeDen[ i ] = 0.0 ;

    if ( cell_type[ cell->type ] == PLASMA ) {

			for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

				if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {

					var->ChargeDen[ i ] += ( var->Qe*config->Species[ jSpecies ].Charge*var->U0[ jSpecies ][ i ] ) ;

				}//if charged species.

			}//species loop.

    }//if plasma region.
    //cout<<var->ChargeDen[ i ]<<endl;
    //var->ChargeDen[ i ] = var->ChargeDen[ i ]*(0.0) ;
  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["ChargeDen"] );
}
void CPoisson::UltraMPPCOmputePermitt( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {
		var->eps_eff[i] = var->eps[i]  ;
  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["effective_permittivity"] );
}
void CPoisson::UltraMPPCOmputeEffectivePermitt( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
	double tmp=0.0 ;

  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;
    tmp = 0.0 ;
    if ( cell_type[ cell->type ] == PLASMA ) {

			for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
				if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {
						tmp += fabs( config->Species[ jSpecies ].Charge ) * (var->Mobi[jSpecies][ i ]) * (var->U0[ jSpecies ][ i ] ) ;
				
				}//if charged species.
			}//species loop.

    }//if plasma region.
		var->eps_eff[i] = var->eps[i] + var->Qe*var->Dt*tmp  ;

  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["effective_permittivity"] );
}
void CPoisson::UltraMPPCOmputeEffectivePermittEleOnly( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*
		∇.(ε∇Φ) = -ρ = -e( Ni-Ne )
	*/
	double tmp=0.0 ;

  for ( int i=0 ; i < plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;
    tmp = 0.0 ;

    if ( cell_type[ cell->type ] == PLASMA ) {

			for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {
				if ( config->Species[ jSpecies ].Type == ELECTRON ) {
						tmp += fabs( config->Species[ jSpecies ].Charge ) * (var->Mobi[jSpecies][ i ]) * (var->U0[ jSpecies ][ i ] ) ;
				
				}//if charged species.
			}//species loop.

    }//if plasma region.
		var->eps_eff[i] = var->eps[i] + var->Qe*var->Dt*tmp  ;

	}//cell loop.
	plasma.syn_parallel_cell_data( var->VarTag["effective_permittivity"] );
}
void CPoisson::UltraMPPCOmputeSurfaceCharge( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	double tmp=0.0 ;
  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;

    tmp = 0.0 ;

    if ( cell_type[ cell->type ] == PLASMA ) {

			for ( int k = 0 ; k < cell->cell_number ; k++ ){
				Cell *cell2 = plasma.get_cell( cell->cell[ k ]->local_id ) ; 
				if ( cell_type[ cell2->type ] == DIELECTRIC ) {

						for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

							if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {

								var->Potential.face[ cell->face[k]->data_id ] += var->Dt*var->Qe*config->Species[jSpecies].Charge
								*fabs( var->U1[ jSpecies ][ i ]*cell->nA[ k ][ 0 ]*(-1.0) 
								+      var->U2[ jSpecies ][ i ]*cell->nA[ k ][ 1 ]*(-1.0) ) ;
							}

						}//end jspecies
				}//discontiuity face
			}//cell2 loop
		}//enf plasma cell
  }//cell loop.
}
void CPoisson::UltraMPPComputeDispCurrentDensity( boost::shared_ptr<CVariable> &var )
{
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		var->DispJD[0][i] = var->eps[ i ]*( var->Ex[ i ] - var->PreEx[ i ] )/var->Dt ;
		var->DispJD[1][i] = var->eps[ i ]*( var->Ey[ i ] - var->PreEy[ i ] )/var->Dt ;

		if( plasma.Mesh.ndim == 3 ) var->DispJD[2][i] = var->eps[ i ]*( var->Ez[ i ] - var->PreEz[ i ] )/var->Dt ;

		var->TotalJD[0][i] = var->DispJD[0][ i ] ;
		var->TotalJD[1][i] = var->DispJD[1][ i ] ;
		if( plasma.Mesh.ndim == 3 ) var->TotalJD[2][i] = var->DispJD[2][ i ] ;
	}//Cell Loop
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
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
		case 0: 	UltraMPPComputePermitt                ( config, var ) ; break ;
		case 1:		UltraMPPComputeEffectivePermitt       ( config, var ) ; break ;
		case 2:		UltraMPPComputeEffectivePermittEleOnly( config, var ) ; break ;
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

	/* face value assigning */
	double voltage= SineVoltage( "POWER", config, var ) ;

	plasma.set_bc_value( MPP_face_tag["POWER" ], voltage, var->Potential.face ) ;
	plasma.set_bc_value( MPP_face_tag["GROUND"],     0.0, var->Potential.face ) ;
	plasma.set_bc_value( MPP_face_tag["NEUMANN"],    0.0, var->Potential.face ) ;

	UltraMPPComputeSurfaceCharge( config, var ) ;

	/* Source B */
	plasma.before_source_term_construction();
	plasma.add_laplacian_source_term_op( var->ChargeDen, var->Potential.face ) ;
	plasma.finish_source_term_construction();

	/* SOLVE */
	plasma.get_solution( var->Potential.current );
	plasma.syn_parallel_cell_data( var->Potential.tag_current );

	its = plasma.get_iteration_number() ;

	plasma.get_gradient(&var->Potential);

	/* inverse gradient for Electric field*/
	for ( int i=0 ; i < plasma.Mesh.cell_number ; i++ ) 
	{
		var->Ex[ i ] = - var->Potential.gradient[0][i] ;
		var->Ey[ i ] = - var->Potential.gradient[1][i] ;
	}
	plasma.syn_parallel_cell_data( var->VarTag["Ex"] );
	plasma.syn_parallel_cell_data( var->VarTag["Ey"] );

	if ( plasma.Mesh.ndim == 3 ) {	
		for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {
			var->Ez[ i ] = - var->Potential.gradient[2][i] ;
		}
		plasma.syn_parallel_cell_data( var->VarTag["Ez"] );
	} 
	var->UltraMPPComputeReducedElectricField();

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
    var->ChargeDen[ i ] = -var->ChargeDen[ i ]/ var->eps_eff[i]  ;
  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["ChargeDen"] );
}
void CPoisson::UltraMPPComputePermitt( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {
		var->eps_eff[i] = var->eps[i]  ;
  }//cell loop.
  plasma.syn_parallel_cell_data( var->VarTag["effective_permittivity"] );
}
void CPoisson::UltraMPPComputeEffectivePermitt( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
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
void CPoisson::UltraMPPComputeEffectivePermittEleOnly( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
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
void CPoisson::UltraMPPComputeSurfaceCharge( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
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
	}
}
double CPoisson::SineVoltage( string FaceType, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	int ft=0 ;
	if     ( FaceType == "POWER"  or FaceType == "SOLID_POWER"  ) ft = POWER ;
	else if( FaceType == "GROUND" or FaceType == "SOLID_GROUND" ) ft = GROUND ;

	map< int, CElectrical>::iterator Iter;
	Iter = config->ElectricalMap.find( ft ) ;
	return Iter->second.BC_Voltage ;
}
// void CPoisson::MatA_SourceB( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
// {

// 	int jFace=0, jCell=0, j=0 ;
// 	double Ad_dPN=0.0, HarmonicMean=0.0, electrode_voltage=0.0, Source=0.0 ;
// 	int cell_index ;
// 	plasma.before_matrix_construction() ;
// 	plasma.before_source_term_construction() ;
	

// 	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

// 		Cell *Cell_i = plasma.get_cell(i) ;

// 		/*--- Loop over electrode cells ---*/
// 		if ( cell_type[ Cell_i->type ] == SOLID_POWER or cell_type[ Cell_i->type ] == SOLID_GROUND ) {

// 			plasma.add_entry_in_matrix     ( i, Cell_i->id, 1.0 ) ;
// 			plasma.add_entry_in_source_term( i, SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ) ;

// 		} else {

// 			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

// 				/*--- Orthogonal term ---*/
// 				j = Cell_i->cell[ k ]->data_id ;

// 				Cell * Cell_j = plasma.get_cell( Cell_i->cell[ k ]->data_id ) ;

// 				cell_index = Cell_i->face_index[k] ;

// 				HarmonicMean = var->eps_eff[ i ]*var->eps_eff[ j ]
// 							/( var->eps_eff[ j ]*( Cell_i->face[k]->dr_c2f[    cell_index] ) 
// 							+  var->eps_eff[ i ]*( Cell_i->face[k]->dr_c2f[1 - cell_index] ) ) ;

// 				Ad_dPN = HarmonicMean * ( Cell_i->face[k]->dA ) ;


// 				plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN ) ;
// 				plasma.add_entry_in_matrix( i, Cell_j->id,  Ad_dPN ) ;

// 				plasma.add_entry_in_source_term( i, Source ) ;

// 			}//Loop over neighbor cells

// 		/*--------------------------------------------------------------*/
// 			for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

// 				Ad_dPN = var->eps_eff[ i ]/(Cell_i->face[k]->dr_c2c)*(Cell_i->face[k]->dA) ;

// 				if ( face_type[ Cell_i->face[k]->type ] == NEUMANN ) {

// 				} else if ( face_type[ Cell_i->face[k]->type ] == GROUND or face_type[ Cell_i->face[k]->type ] == POWER ){

// 					electrode_voltage = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id), config, var ) ;

// 					plasma.add_entry_in_matrix     ( i,  Cell_i->id, -Ad_dPN ) ;
// 					plasma.add_entry_in_source_term( i, -(electrode_voltage)*Ad_dPN ) ;

// 				}else if ( plasma.get_face_typename( Cell_i->face[ k ]->data_id)  == "DIELECTRIC" ) {

// 					cout<<"Boundary face don't have DIELECTRIC, pls check w/ K.-L. Chen-2"<<endl;
// 					exit(1) ;

// 				}else{
// 					cout<<"error-\" solver_poisson.cpp-3\""<<endl;
// 				}

// 			}//Loop over boundary face cells
// 		/*--------------------------------------------------------------*/
// 		}
// 		plasma.add_entry_in_source_term( i, var->ChargeDen[ i ]*Cell_i->volume ) ;
// 	}//Cell Loop
// 	plasma.finish_matrix_construction() ;
// 	plasma.finish_source_term_construction() ;
// }
// void CPoisson::ComputeGradient( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
// {
	
// 	int j=0, NeighborCellIndex=0 ;
// 	double dVar=0.0, Gx=0.0, Gy=0.0, Gz=0.0, BC_Value=0.0, L=0.0, R=0.0;

// 	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

// 		Cell *Cell_i = plasma.get_cell(i) ;

// 		Gx = Gy = Gz = 0.0 ;

// 		if ( cell_type[ Cell_i->type ] == POWER or cell_type[ Cell_i->type ] == GROUND ) {

// 			var->Ex[ i ] = 0.0 ;
// 			var->Ey[ i ] = 0.0 ;
// 			if( plasma.Mesh.ndim == 3 )  var->Ez[ i ] = 0.0 ;

// 		} else {//could be dielectric or plasma

// 			/*--- Loop over neighbor "cells" ---*/
// 			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

// 				j = Cell_i->cell[ k ]->data_id ;
// 				Cell * Cell_j = plasma.get_cell( j ) ;

// 				if ( cell_type[ Cell_i->type ] != cell_type[ Cell_j->type ]  ) {//For discontinued face

// 					//If neighboring cell is "electrode", the boundary value should calculate by SineVoltage function. 
// 					if ( cell_type[ Cell_j->type ] == POWER or cell_type[ Cell_j->type ] == GROUND ) {

// 						BC_Value = SineVoltage( plasma.get_cell_typename( Cell_j->data_id ), config, var ) ;
// 						dVar = BC_Value - var->Potential.current[ i ] ;

// 					//If neighboring cell is "dielectric", the boundary value should calculate by formula. 

// 					} else {

// 						cout<<"error for calculate gradient"<<endl ;
// 						// BC_Value = m->PFM_CELL[ i ][ k ].dNPf*m->PFM_CELL[ i ][ k ].dPPf * m->PFM_CELL[ i ][ k ].SurfaceCharge
// 						//  		 / ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) ;

// 						// BC_Value +=( var->Potential.current[ j ]*var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Potential.current[ i ]*var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf ) 
// 						//  		/ ( var->Eps[ j ]*m->PFM_CELL[ i ][ k ].dPPf + var->Eps[ i ]*m->PFM_CELL[ i ][ k ].dNPf );

									
// 						// dVar = BC_Value - var->Potential.current[ i ] ;
// 					}

// 				} else {

// 					dVar = var->Potential.current[ j ] - var->Potential.current[ i ] ;

// 				}
// 				Gx = Gx + var->LSQ_Cx[ k ][ i ]*dVar ;
// 		    Gy = Gy + var->LSQ_Cy[ k ][ i ]*dVar ;
// 			}
// 			//cout<<"Loop over boundary faces"<<endl;


// 			/*--- Loop over boundary faces ---*/
// 			for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

// 				if ( face_type[ Cell_i->face[ k ]->type ]  == NEUMANN ) {

// 					BC_Value = var->Potential.current[ i ] ;

// 				} else if( face_type[ Cell_i->face[ k ]->type ]  != PLASMA ) {

// 					BC_Value = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id) , config, var ) ;

// 				} else {

// 					cout<<"ERR@Poisson 388"<<endl ;

// 				}
				
// 				dVar = BC_Value - var->Potential.current[ i ] ;
// 				//cout<<var->Potential.current[ i ]<<endl;
// 				Gx = Gx + var->LSQ_Cx[ k ][ i ]*dVar ;
// 		    Gy = Gy + var->LSQ_Cy[ k ][ i ]*dVar ;

// 			}
// 			var->Ex[ i ] = Gx ;
// 			var->Ey[ i ] = Gy ;
// 			if( plasma.Mesh.ndim == 3 ) var->Ez[ i ] = Gz ;

// 		}//For Dielectric or Plasma

// 		var->Ex[ i ] = (-1.0)*var->Ex[ i ] ;
// 		var->Ey[ i ] = (-1.0)*var->Ey[ i ] ;
// 	}//Loop over all cells
// 	plasma.syn_parallel_cell_data( var->VarTag["Ex"] );
// 	plasma.syn_parallel_cell_data( var->VarTag["Ey"] );
// }
// void CPoisson::SOLVE_TEST( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
// {
// 	UpdateElectricalMap( config, var ) ;
	
// 	/*--- compute the transport coefficient (permittivity). ---*/
// 	switch ( config->Equation[ POISSON ].Equation ) {
// 		case 0: 	UltraMPPComputePermitt                ( config, var ) ; break ;
// 		case 1:		UltraMPPComputeEffectivePermitt       ( config, var ) ; break ;
// 		case 2:		UltraMPPComputeEffectivePermittEleOnly( config, var ) ; break ;
// 		default: 
// 		if( mpi_rank == 0 ) cout << "No such equation option, Pls contact K.-L. Chen " << endl ; exit(1) ; break ;
// 	}

// 	/*--- Calculate the net charge density for poisson's source term ---*/
// 	UltraMPPComputeNetCharge( config, var ) ;

// 	/*--- set the cell paramtert ---*/
// 	//plasma.set_cell_property_parameter( var->eps_eff ) ;


// 	MatA_SourceB( config, var ) ;

// 	/* SOLVE */
// 	plasma.get_solution( var->Potential.current );
// 	plasma.syn_parallel_cell_data( var->Potential.tag_current );

// 	its = plasma.get_iteration_number() ;

// 	ComputeGradient( config, var ) ;
// 	UltraMPPComputeDispCurrentDensity( var ) ;

// }

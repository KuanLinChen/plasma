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
	}
	its=0 ;
	//fstream  file ; 
	left .open("left.dat" ) ; 
	right.open("right.dat") ; 
	cout<<"openfile"<<endl;
	//left<<"AA" ;
	//right<<"BB" ;
	//PetscEnd();
	data_count = 0;
	stable_coefficient = 1.0 ;
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

//	cout<<"B"<<endl;PetscEnd();

	
	/* face value assigning */
	double voltage= SineVoltage( "POWER", config, var ) ;
	//cout<<voltage<<endl;
	//PetscEnd();
	plasma.set_bc_value( MPP_face_tag["POWER" ], voltage, var->Potential.face ) ;
	plasma.set_bc_value( MPP_face_tag["GROUND"],     0.0, var->Potential.face ) ;
	plasma.set_bc_value( MPP_face_tag["NEUMANN"],    0.0, var->Potential.face ) ;
	
	if( plasma.Mesh.ndim == 3 ){
		UltraMPPComputeSurfaceCharge( config, var ) ;
	} else {
		UltraMPPComputeSurfaceCharge_2d( config, var ) ;
	}

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
	for ( int i=0 ; i < plasma.Mesh.cell_number ; i++ ) {
		var->Ex[ i ] = - var->Potential.gradient[0][i] ;
		var->Ey[ i ] = - var->Potential.gradient[1][i] ;
		//cout<<var->Potential.gradient[0][i] <<endl;
	}
	plasma.syn_parallel_cell_data( var->VarTag["Ex"] );
	plasma.syn_parallel_cell_data( var->VarTag["Ey"] );

	if( plasma.Mesh.ndim == 3 ) {
		for ( int i=0 ; i < plasma.Mesh.cell_number ; i++ ) {
			var->Ez[ i ] = - var->Potential.gradient[2][i] ;
		}
		plasma.syn_parallel_cell_data( var->VarTag["Ez"] );
	}
	/*Compute real net charge density [C/m^3]*/
	for ( int i=0 ; i < plasma.Mesh.cell_number ; i++ ) {
		var->RealChargeDen[ i ]   = - var->ChargeDen[ i ] * vacuum_permittivity / var->Qe ;
	}

	var->UltraMPPComputeReducedElectricField() ;

	UltraMPPComputeDispCurrentDensity( var ) ;

}
void CPoisson::UltraMPPComputeNetCharge( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	/*Only the 'plasma' region. */
  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;
    var->ChargeDen[ i ] = 0.0 ;

    if ( cell->type == MPP_cell_tag[ "PLASMA" ] ) {

			for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

				if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {

					var->ChargeDen[ i ] += ( var->Qe*config->Species[ jSpecies ].Charge*var->U0[ jSpecies ][ i ] ) ;

				}//if charged species.

			}//species loop.

    }//if plasma region.
    var->ChargeDen[ i ] = -var->ChargeDen[ i ] / var->eps_eff[i]  ;
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

  for ( int i=0 ; i< plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;
    tmp = 0.0 ;
    if ( cell->type == MPP_cell_tag[ "PLASMA" ] ) {

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

  for ( int i=0 ; i< plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;
    tmp = 0.0 ;
    if ( cell->type == MPP_cell_tag[ "PLASMA" ] ) {

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
	double sign_nA = 0.0 ;
  for ( int i=0 ; i<plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;

    tmp = 0.0 ;

    if ( cell->type == MPP_cell_tag[ "PLASMA" ] ) {

			for ( int k = 0 ; k < cell->cell_number ; k++ ){

				sign_nA = pow( -1.0 , cell->face_index[k] ) ;
				
				Cell *cell2 = plasma.get_cell( cell->cell[ k ]->data_id ) ; 

				if ( cell2->type == MPP_cell_tag[ "DIELECTRIC" ] ) {

						for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

							if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {

								var->surface_charge[ cell->face[k]->data_id ] += var->Dt*var->Qe*config->Species[jSpecies].Charge
								*fabs( var->U1[ jSpecies ][ i ]*sign_nA*cell->face[k]->nA[0] +
								   		 var->U2[ jSpecies ][ i ]*sign_nA*cell->face[k]->nA[1] +
								   		 var->U3[ jSpecies ][ i ]*sign_nA*cell->face[k]->nA[2] ) ;
							}
						}//end jspecies

						var->Potential.face[ cell->face[k]->data_id ] = var->surface_charge[ cell->face[k]->data_id ] ;
						
				}//discontiuity face
			}//cell2 loop
		}//enf plasma cell
  }//cell loop.
}
void CPoisson::UltraMPPComputeSurfaceCharge_2d( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	//double surface_charge=0.0 ;
	double sign_nA = 0.0 ;
	data_count = data_count + 1 ;
  for ( int i=0 ; i < plasma.Mesh.cell_number ; i++ ) {

    Cell *cell = plasma.get_cell( i ) ;

    //surface_charge = 0.0 ;

    if ( cell->type == MPP_cell_tag[ "PLASMA" ] ) {
			for ( int k = 0 ; k < cell->cell_number ; k++ ){
				sign_nA = pow( -1.0 , cell->face_index[k] ) ;
				
				Cell *cell2 = plasma.get_cell( cell->cell[ k ]->data_id ) ; 
				if ( cell2->type == MPP_cell_tag[ "DIELECTRIC" ] ) {

						for ( int jSpecies = 0 ; jSpecies < config->TotalSpeciesNum ; jSpecies++ ) {

							if ( config->Species[ jSpecies ].Type == ELECTRON or config->Species[ jSpecies ].Type == ION ) {

								 var->surface_charge[ cell->face[k]->data_id ] += var->Dt*var->Qe*config->Species[jSpecies].Charge
								*fabs( var->U1[ jSpecies ][ i ]*sign_nA*cell->face[k]->nA[0] +
								   		 var->U2[ jSpecies ][ i ]*sign_nA*cell->face[k]->nA[1] ) ;
							}

						}//end jspecies
						//var->Potential.face[ cell->face[k]->data_id ] = var->Potential.face[ cell->face[k]->data_id ]*(-1.0) ;

						var->Potential.face[ cell->face[k]->data_id ] = var->surface_charge[ cell->face[k]->data_id ] ;
						//cout<<"AAAAAAAAAAa"<<endl;
				}//discontiuity face
			}//cell2 loop
		}//enf plasma cell
  }//cell loop.
}
void CPoisson::UltraMPPComputeDispCurrentDensity( boost::shared_ptr<CVariable> &var )
{
	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		var-> DispJD[0][i] = var->eps[ i ]*( var->Ex[ i ] - var->PreEx[ i ] )/var->Dt ;
		var-> DispJD[1][i] = var->eps[ i ]*( var->Ey[ i ] - var->PreEy[ i ] )/var->Dt ;
		var-> DispJD[2][i] = var->eps[ i ]*( var->Ez[ i ] - var->PreEz[ i ] )/var->Dt ;
		var->TotalJD[0][i] = var->DispJD[0][ i ] ;
		var->TotalJD[1][i] = var->DispJD[1][ i ] ;
		var->TotalJD[2][i] = var->DispJD[2][ i ] ;
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
void CPoisson::MatA_SourceB( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{

	int jFace=0, jCell=0, j=0 ;
	double Ad_dPN=0.0, HarmonicMean=0.0, electrode_voltage=0.0, Source=0.0, phi_s = 0.0, SurfaceQ=0.0, tmp=0.0 ;
	double dL, dR ;
	int cell_index ;
	plasma.before_matrix_construction() ;
	plasma.before_source_term_construction() ;
	

	for( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell *Cell_i = plasma.get_cell(i) ;

		/*--- Loop over electrode cells ---*/
		if ( cell_type[ Cell_i->type ] == SOLID_POWER or cell_type[ Cell_i->type ] == SOLID_GROUND ) {

			//plasma.add_entry_in_matrix     ( i, Cell_i->id, 1.0 ) ;
			//plasma.add_entry_in_source_term( i, SineVoltage( plasma.get_cell_typename( Cell_i->data_id ), config, var ) ) ;
			cout<<"err-1"<<endl; PetscEnd() ;

		} else {

			//This part include the plasma and Dielectric.
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				/*--- Orthogonal term ---*/
				j = Cell_i->cell[ k ]->data_id ;

				Cell * Cell_j = plasma.get_cell( Cell_i->cell[ k ]->data_id ) ;

				cell_index = Cell_i->face_index[k] ;
				dL = Cell_i->face[k]->dr_c2f[  cell_index] ;
				dR = Cell_i->face[k]->dr_c2f[1-cell_index] ;

				HarmonicMean = var->eps_eff[ i ]*var->eps_eff[ j ]
							/( var->eps_eff[ j ]*dL +  var->eps_eff[ i ]*dR ) ;

				Ad_dPN = HarmonicMean * ( Cell_i->face[k]->dA )  ;


				plasma.add_entry_in_matrix( i, Cell_i->id, -Ad_dPN ) ;
				plasma.add_entry_in_matrix( i, Cell_j->id,  Ad_dPN ) ;


				//if ( (Cell_i->type == MPP_cell_tag["PLASMA"] and Cell_j->type == MPP_cell_tag["DIELECTRIC"]) or  ) 
				//{
					SurfaceQ = var->Potential.face[ Cell_i->face[k]->data_id ] ;

					tmp = var->eps_eff[ i ]*dR /( var->eps_eff[ j ] * dL + var->eps_eff[ i ]*dR );

					plasma.add_entry_in_source_term( i, -SurfaceQ*tmp*Cell_i->face[k]->dA ) ;

				//}

			}//Loop over neighbor cells

		/*--------------------------------------------------------------*/
			for( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				Ad_dPN = var->eps_eff[ i ]/(Cell_i->face[k]->dr_c2c)*(Cell_i->face[k]->dA)   ;

				if ( Cell_i->face[ k ]->type == MPP_face_tag[ "NEUMANN" ] ) {

				} else if ( Cell_i->face[ k ]->type == MPP_face_tag[ "GROUND" ] or Cell_i->face[ k ]->type == MPP_face_tag[ "POWER" ]  ){

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
		plasma.add_entry_in_source_term( i, var->ChargeDen[ i ]*Cell_i->volume*var->eps_eff[i] ) ;
	}//Cell Loop
	plasma.finish_matrix_construction() ;
	plasma.finish_source_term_construction() ;
}
void CPoisson::ComputeGradient( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	//cout<<"This will be remove in the feture. @ solver_poisson.cpp"<<endl;
	int j=0, NeighborCellIndex=0 ;
	double dVar=0.0, Gx=0.0, Gy=0.0, BC_Value=0.0, dL=0.0, dR=0.0, SurfaceQ, phi_s ;
	int cell_index=0 ;

	for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {

		Cell *Cell_i = plasma.get_cell(i) ;

		Gx = Gy = 0.0 ;

		if ( cell_type[ Cell_i->type ] == POWER or cell_type[ Cell_i->type ] == GROUND ) {

			//var->Ex[ i ] = 0.0 ;
			//var->Ey[ i ] = 0.0 ;
			cout<<"err-2"<<endl;PetscEnd();

		} else {//could be dielectric or plasma

			/*--- Loop over neighbor "cells" ---*/
			for ( int k = 0 ; k < Cell_i->cell_number ; k++ ) {

				j = Cell_i->cell[ k ]->data_id ;
				Cell * Cell_j = plasma.get_cell( j ) ;
				cell_index = Cell_i->face_index[k] ;

				dL = Cell_i->face[k]->dr_c2f[  cell_index] ;
				dR = Cell_i->face[k]->dr_c2f[1-cell_index] ;


				if ( cell_type[ Cell_i->type ] != cell_type[ Cell_j->type ]  ) {//For discontinued face

					if ( Cell_j->type == MPP_cell_tag["DIELECTRIC"] or Cell_j->type == MPP_cell_tag["PLASMA"] ) {

						SurfaceQ = var->Potential.face[ Cell_i->face[k]->data_id ] ;

						phi_s =	  dL*dR/( var->eps_eff[ j ]*dL + var->eps_eff[ i ]*dR ) * SurfaceQ ;

						phi_s += ( var->Potential.current[ j ]*var->eps_eff[ j ]*dL + var->Potential.current[ i ]*var->eps_eff[ i ]*dR ) / ( var->eps_eff[ j ]*dL + var->eps_eff[ i ]*dR );


						// BC_Value = SineVoltage( plasma.get_cell_typename( Cell_j->data_id ), config, var ) ;
						// //dVar = BC_Value - var->Potential.current[ i ] ;
						if ( pow( -1.0, Cell_i->face_index[ k ] ) * Cell_i->face[ k ]->nA[0] > 0.0) {
						 	Gx = Gx + (phi_s - var->Potential.current[ i ] )/dL*0.5 ;
						} else{
						 	Gx = Gx+  (var->Potential.current[ i ] - phi_s )/dL*0.5 ;
						}

					} else {

						cout<<"error for calculate gradient"<<endl ;
						
					}

					//	PetscEnd();

				} else {
					//cout<<"face normal: "<<Cell_i->face[k]->nA[0]<<endl;
					if (  pow( -1.0, Cell_i->face_index[ k ] ) * Cell_i->face[ k ]->nA[0] > 0.0 ){

						Gx = Gx + (var->Potential.current[ j ] - var->Potential.current[ i ])/Cell_i->face[k]->dr_c2c*0.5 ;
						//cout<<"eface: "<<(var->Potential.current[ j ] - var->Potential.current[ i ])/Cell_i->face[k]->dr_c2c*0.5<<endl;
					} else {

						Gx = Gx + (var->Potential.current[ i ] - var->Potential.current[ j ])/Cell_i->face[k]->dr_c2c*0.5 ;
						//cout<<"wface: "<<(var->Potential.current[ i ] - var->Potential.current[ j ])/Cell_i->face[k]->dr_c2c*0.5<<endl;
					}
				}

				Gy = 0.0 ;
			}
			//cout<<"Loop over boundary faces"<<endl;
			//cout<<Gx<<endl;

			/*--- Loop over boundary faces ---*/
			for ( int k = Cell_i->cell_number ; k < Cell_i->face_number ; k++ ) {

				if ( Cell_i->face[ k ]->type == MPP_face_tag[ "NEUMANN" ] ) {


				} else {

					cell_index = Cell_i->face_index[k] ;

					BC_Value = SineVoltage( plasma.get_face_typename( Cell_i->face[ k ]->data_id) , config, var ) ;

					//cout<<plasma.get_face_typename( Cell_i->face[ k ]->data_id)<<"\t"<<BC_Value<<endl;

					if ( pow( -1.0, Cell_i->face_index[ k ] ) * Cell_i->face[ k ]->nA[0] > 0.0 ){

						Gx = Gx + ( BC_Value - var->Potential.current[ i ]  )/Cell_i->face[k]->dr_c2f[cell_index]*0.5 ;
						//cout<<"A1"<<endl;
					} else {

						Gx = Gx + ( var->Potential.current[ i ] - BC_Value )/Cell_i->face[k]->dr_c2f[cell_index]*0.5 ;
						//cout<<"A2"<<endl;
					}

				}
			}

			var->Ex[ i ] = Gx ;
			var->Ey[ i ] = 0.0 ;

		}//For Dielectric or Plasma

		var->Ex[ i ] = (-1.0)*var->Ex[ i ] ;
		//var->Ey[ i ] = (-1.0)*var->Ey[ i ] ;
		//cout<<endl;
	}//Loop over all cells
	//PetscEnd();
	plasma.syn_parallel_cell_data( var->VarTag["Ex"] );
	plasma.syn_parallel_cell_data( var->VarTag["Ey"] );
}
void CPoisson::SOLVE_TEST( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  )
{
	//cout<<"This will be remove in the feture. @ solver_poisson.cpp"<<endl;
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
	UltraMPPComputeSurfaceCharge_2d( config, var ) ;

	MatA_SourceB( config, var ) ;

	/* SOLVE */
	plasma.get_solution( var->Potential.current );
	plasma.syn_parallel_cell_data( var->Potential.tag_current );

	its = plasma.get_iteration_number() ;

	ComputeGradient( config, var ) ;
	var->UltraMPPComputeReducedElectricField();
	UltraMPPComputeDispCurrentDensity( var ) ;

}

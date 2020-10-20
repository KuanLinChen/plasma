#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <mpi.h>
#include "json.hpp"
#include "ultraMPP.h"

using namespace std ;
using nlohmann::json ;

int cpu_size, myid ;

int main( int argc, char * argv[] )
{
	ultraMPP target, source, read ;
	
	//Initial ultraMPP object.
	target.initial( argc,argv, &myid, &cpu_size ) ;
	source.initial( argc,argv, &myid, &cpu_size ) ;
	  read.initial( argc,argv, &myid, &cpu_size ) ;

	source.load_mesh( "source.json" ) ;
	target.load_mesh( "target.json" ) ;
	  read.load_mesh( "read.json"   ) ;

	double atom_mass = 6.6335209E-26 ; //for argon
	double 	*rho_source, *rho_target, *output, *rho_read ;

	int tag_rho_source =	source.set_parallel_cell_data( &rho_source, "Rho" ) ;

	int tag_rho_target =	target.set_parallel_cell_data( &rho_target, "Rho" ) ;
	int tag_rho_output =	target.set_parallel_cell_data( &output    , "N_Ar" ) ;

	int tag_rho_read   =	  read.set_parallel_cell_data( &rho_read  , "N_Ar" ) ;

	source.set_initial_variable( tag_rho_source ) ;
	source.get_initial_value( "test.dat" ) ;
	cout<<"Read test file"<<endl;

	target.syn_parallel_cell_data( tag_rho_target, tag_rho_source ) ;
	cout<<"interpolated to new domain"<<endl;
	Cell *iCell ;

	//Convert rho to number density.
	for ( int i=0 ; i < target.Mesh.cell_number ; i++ ) {
		iCell = target.get_cell(i) ;
		output[ i ] = rho_target[i]/atom_mass ;
	}



	//source 
	source.set_output  ( "source" ) ;
	source.set_output  ( tag_rho_source ) ;
	source.write_output("source" ) ;
	cout<<"write source.dat file"<<endl;

	//target
	target.set_output  ( "target" ) ;
	target.set_output  ( tag_rho_target ) ;
	target.set_output  ( tag_rho_output ) ;
	target.write_output("target" ) ;
	cout<<"write target.dat file"<<endl;

	//read
	read.set_initial_variable( tag_rho_read ) ;
	read.get_initial_value( "target.dat" ) ;
	cout<<"read target.dat file"<<endl;

	read.set_output  ( "read" ) ;
	read.set_output  ( tag_rho_read ) ;
	read.write_output("read" ) ;
	cout<<"write read.dat file"<<endl;

	return 0 ;
}

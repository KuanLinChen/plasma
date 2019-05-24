#include "main.h"
#include <string>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <iostream>

#ifdef USEMPI
	#include <mpi.h>
#endif

using namespace std;

string generate_id ()
{
	int i, mpiId, mpiSize ;
	char buf[37] ;
	string str ;

	str.clear();

	boost::uuids::uuid id = boost::uuids::random_generator()();
	str = boost::uuids::to_string( id ) ;


	#ifdef USEMPI
	MPI_Comm_size ( MPI_COMM_WORLD, &mpiSize ) ;
	MPI_Comm_rank ( MPI_COMM_WORLD, &mpiId) ;
	if ( mpiSize > 0 )
	{
		strcpy( buf, str.c_str() );
		MPI_Bcast ( buf, str.size(), MPI_CHAR, 0, MPI_COMM_WORLD ) ;
		string s( buf ) ;
		str = s ;
	}
	#endif

	return  str ;

}
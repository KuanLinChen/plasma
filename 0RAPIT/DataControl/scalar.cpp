#include <iostream>
#include "scalar.h"
#include "datapool.h"
#include "sys_log.h"
//#include "fvm.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include "uuid_functions.h"
using namespace std;


bool operator< ( int & i , Scalar & s ) { return ( i < s.local_data_number ) ;}
bool operator< ( Scalar & s , int & i ) { return ( i > s.local_data_number ) ;}
bool operator> ( int & i , Scalar & s ) { return ( i > s.local_data_number ) ;}
bool operator> ( Scalar & s , int & i ) { return ( i < s.local_data_number ) ;}
bool operator<= ( int & i , Scalar & s) { return ( i <= s.local_data_number );}
bool operator<= ( Scalar & s , int & i) { return ( i >= s.local_data_number );}
bool operator>= ( int & i , Scalar & s) { return ( i >= s.local_data_number );}
bool operator>= ( Scalar & s , int & i) { return ( i <= s.local_data_number );}
bool operator== ( int & i , Scalar & s) { return ( i == s.local_data_number );}
bool operator== ( Scalar & s , int & i) { return ( i == s.local_data_number );}
bool operator!= ( int & i , Scalar & s) { return ( i == s.local_data_number );}
bool operator!= ( Scalar & s , int & i) { return ( i == s.local_data_number );}


DataManager Scalar::datamanager ;
set<string> Scalar::ID_bank ;

Scalar::Scalar ( boost::shared_ptr<Domain> p, int t )
{
	int i ;

	data_number = 0 ;
	local_data_number = 0 ;
	domain_ptr = NULL ;
	share_memory_flag = false ;
	flag_collecting_face =false ;

	domain_ptr = p ;
	type = t ;

	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;

	if ( type == CELL )
	{
		data_number = p->local_cell_number + p->ghost_cell_number ;
		local_data_number = p->local_cell_number  ;
	} else if ( type == NODE )
	{
		data_number = p->local_node_number + p->ghost_node_number ;
		local_data_number = p->local_node_number;
	} else if ( type == FACE )
	{
		data_number = p->local_face_number + p->ghost_face_number;
		local_data_number = p->local_face_number ;
	}

	data = boost::shared_array <double>( new double [ data_number ] ) ;
	
	Log( p->comm ).MPITagDump( logLEVEL4 ) << "Scalar @" << data.get() << " (" <<  data_number << ")" ;
	for ( i = 0; i < data_number ; i++ )
		data[i] = 0. ; // 1. * i  ;

}

Scalar::Scalar ( int n )
{
	int i ;

	data_number = 0 ;
	local_data_number = 0 ;
	domain_ptr = NULL ;
	share_memory_flag = false ;
	flag_collecting_face =false ;

	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;

	data_number = n ;
	local_data_number = n ;
	data = boost::shared_array<double> (  new double [ data_number ] );
	for ( i = 0; i < data_number; i++ )
		data[i] = 0. ; 
}

/*!
	The empty scalar object
*/

Scalar::Scalar ( )
{
	data_number = 0 ;
	local_data_number = 0 ;
	domain_ptr = NULL ;
	share_memory_flag = false ;
	flag_collecting_face =false ;

	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;
}

Scalar::~Scalar ( )
{
	ID_bank.erase( ID ) ;
}

Scalar::Scalar ( boost::shared_ptr<Domain> p, int t, string n )
{
	int i ;

	data_number = 0 ;
	local_data_number = 0 ;
	domain_ptr = NULL ;
	share_memory_flag = false ;
	flag_collecting_face =false ;

	name = n ;

	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;
	
	domain_ptr = p ;
	type = t ;
	
	if ( type == CELL )
	{
		data_number = p->local_cell_number + p->ghost_cell_number ;
		local_data_number = p->local_cell_number  ;
	} else if ( type == NODE )
	{
		data_number = p->local_node_number + p->ghost_node_number ;
		local_data_number = p->local_node_number;
	} else if ( type == FACE )
	{
		data_number = p->local_face_number + p->ghost_face_number ;
		local_data_number = p->local_face_number ;
	}

	data = boost::shared_array <double>( new double [ data_number ] );
	Log( p->comm ).MPITagDump( logLEVEL4 ) << "Scalar: " << name << " @" << data.get() << " (" <<  data_number << ")" ;

	for ( i = 0; i < data_number; i++ )
		data[i] = 0. ;
}

void Scalar::initial ( boost::shared_ptr<Domain> p, int t, string n )
{
	int i ;

	name = n ;

	domain_ptr = p ;
	type = t ;

	// Skip if not belong to this communicator
	if ( p->comm == MPI_COMM_NULL ) return;

	if ( type == CELL )
	{
		data_number = p->local_cell_number + p->ghost_cell_number ;
		local_data_number = p->local_cell_number  ;
	}
	if ( type == NODE )
	{
		data_number = p->local_node_number + p->ghost_node_number ;
		local_data_number = p->local_node_number;
	}
	if ( type == FACE )
	{
		data_number = p->local_face_number + p->ghost_face_number ;
		local_data_number = p->local_face_number  ;
	}

	data = boost::shared_array<double> ( new double [data_number] ) ;

	Log( p->comm ).MPITagDump( logLEVEL4 ) << "Scalar: " << name << " @" << data.get() << " (" <<  data_number << ")" ;
	
	for ( i = 0; i < data_number; i++ )
		data[i] = 0. ;

}

void Scalar::initial ( boost::shared_ptr<Domain> p, int t )
{
	int i ;

	domain_ptr = p ;
	type = t ;
	
	// Skip if not belong to this communicator
	if ( p->comm == MPI_COMM_NULL ) return;
	
	if ( type == CELL )
	{
		data_number = p->local_cell_number + p->ghost_cell_number ;
		local_data_number = p->local_cell_number  ;
	} else if ( type == NODE )
	{
		data_number = p->local_node_number + p->ghost_node_number ;
		local_data_number = p->local_node_number;
	} else if ( type == FACE )
	{
		data_number = p->local_face_number + p->ghost_face_number;
		local_data_number = p->local_face_number ;
	}

	data = boost::shared_array<double> ( new double [data_number] ) ;

	Log( p->comm ).MPITagDump( logLEVEL4 ) << "Scalar: " << name << " @" << data.get() << " (" <<  data_number << ")" ;
	
	for (int i = 0; i < data_number; i++ )
		data[i] = 0. ; // 1. * i  ;

}


void Scalar::collecting_face ()
{
	int i, buffer;

	flag_collecting_face = true ;
	
	// Skip if not belong to this communicator
	if ( domain_ptr->comm == MPI_COMM_NULL ) return;
	
	if ( type != FACE ) 
	{ 
		Log().TagDump ( logERROR )  << "Not supported for the collecting of faces." ;
		return ;
	}

	buffer = collected_faces.size() ;

	for ( i = 0  ; i < domain_ptr->local_face_number + domain_ptr->ghost_face_number ; i++ ) 
	{
		if ( data[i] > 0.0 )
			collected_faces.push_back ( i ) ;
	}

	// Remove current ID and generate another one
	ID_bank.erase( ID ) ;
	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;
	
	Log( domain_ptr->comm ).MPITagDump ( logLEVEL4 ) << ID <<  " There are originally " << buffer << " face(s) and " << collected_faces.size() - buffer << " face(s) added." ;	
	

}


void Scalar::collecting_face ( double *d )
{
	int i, buffer;

	flag_collecting_face = true ;
	
	// Skip if not belong to this communicator
	if ( domain_ptr->comm == MPI_COMM_NULL ) return;
	
	if ( type != FACE ) 
	{
		Log().TagDump ( logERROR )  << "Not supported for the collecting of faces." ;
		return ;
	}

	buffer = collected_faces.size() ;

	for ( i = 0  ; i < domain_ptr->local_face_number + domain_ptr->ghost_face_number ; i++ ) 
	{
		if ( d[i] > 0.0 )
			collected_faces.push_back ( i ) ;
	}

	// Remove current ID and generate another one
	ID_bank.erase( ID ) ;
	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;	
	
	Log().TagDump ( logLEVEL4 ) << ID << " There are originally " << buffer << " face(s) and " << collected_faces.size() - buffer << " face(s) added." ;

}


void Scalar::collecting_face (  string s  )
{
	int i, buffer;

	flag_collecting_face = true ;
	
	// Skip if not belong to this communicator
	if ( domain_ptr->comm == MPI_COMM_NULL ) return;
	
	if ( type != FACE )
	{
		Log().TagDump ( logERROR )  << "Not supported for the collecting of faces." ;
		return;	
	} 

	buffer = -1000 ;
	for ( i = 0 ; i < domain_ptr->BCFace_type_number ; i++ )
	{
		if (  boost::algorithm::iequals ( domain_ptr->BCFace_typename [i], s ) )
			buffer = i ; 
	}

	if ( buffer == -1000 )
	{
		Log().TagDump ( logERROR )  << "Face type not found: " << s ;
		return ; 
	}

	for ( i = 0  ; i < domain_ptr->BCFace[buffer].size() ; i++ ) 
	{
		collected_faces.push_back ( domain_ptr->BCFace[buffer][i] ) ;
	}

	// Remove current ID and generate another one
	ID_bank.erase( ID ) ;
	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;	
	
	Log().TagDump ( logLEVEL4 ) << ID << " There are originally " << collected_faces.size() - domain_ptr->BCFace[buffer].size() << " face(s) and " << domain_ptr->BCFace[buffer].size() << " face(s) added." ;	
}


void Scalar::reset_collecting_face ( )
{
	int buffer ;
	buffer = collected_faces.size() ;
	collected_faces.clear() ;
	flag_collecting_face = false ;

	// Remove current ID and generate another one
	ID_bank.erase( ID ) ;
	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;	
	Log().TagDump ( logLEVEL4 ) << ID << " There are " << buffer  << " face(s) has been remove to " << collected_faces.size() << " face." ;	
}

void Scalar::zero (  )
{
	int i ;
	for ( i = 0; i < data_number; i++ )
		data[i] = 0. ;
}


void Scalar::share_memory_from ( const Scalar & c )
{
	int i ;
	share_memory_flag	= true ;
	data_number			= c.data_number ;
	local_data_number	= c.local_data_number ;
	domain_ptr			= c.domain_ptr ;
	name				= c.name ;
	type				= c.type ;

	data.reset() ;
	data = c.data ;
}


Scalar & Scalar::operator= ( const Scalar & c ) 
{
	int id; 
	boost::shared_array<double> target, source ;
	boost::shared_ptr<Domain> domain_target, domain_source ;

	target = data;
	source = c.data ;
	domain_target = domain_ptr ;
	domain_source = c.domain_ptr ;


	if ( type == CELL && c.type == CELL )
	{
		datamanager.cell_interpolation ( domain_target.get(), domain_source.get(), target, source ) ;
	} else if (  type == FACE && c.type == FACE )
	{
		if ( c.flag_collecting_face )
		{
			datamanager.face_interpolation ( domain_target.get(), domain_source.get(), c.ID, target, source, c.collected_faces ) ;
		} else
		{
			datamanager.face_interpolation ( domain_target.get(), domain_source.get(), target, source ) ;
		}

	} else  
	{
		//cout << "Interpolation type of " << type << " " << c.type << " is not supported. "  <<  endl;
		Log().TagDump( logERROR ) << "Interpolation type of " << type << " " << c.type << " is not supported. " ;
	}
	
	return *this ;
}

Scalar & Scalar::operator= ( const double*  c ) 
{
	int i ;
	for ( i = 0 ; i < local_data_number ; i++)
		data[i] = c[i] ;
}

void Scalar::remove_relation ( const Scalar & c ) 
{
	int id ;

	boost::shared_ptr<Domain> domain_target, domain_source ;
	domain_target = domain_ptr ;
	domain_source = c.domain_ptr ;

}



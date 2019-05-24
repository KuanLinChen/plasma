#include "io.h"
#include "scalar.h"
#include "sys_log.h"
#include "tecplot_output.h"
#include "vtk_output.h"

IO::IO()
{
	is_new_file = true ;
	domain_ptr = NULL;
	dimension = 2 ;
	output_type = IO_VTK;
	if ( output_type == IO_TECPLOT or output_type == IO_TECPLOT1D )
		output_filename = "output_" + NowTime() + ".dat" ;
	else if ( output_type == IO_VTK )
		output_filename = "output_" + NowTime() + ".vtk" ;
}

IO::~IO()
{

}

IO::IO( string _s )
{
	is_new_file = true ;
	domain_ptr = NULL;
	dimension = 2 ;
	output_type = IO_VTK;
	set_filename ( _s ) ;
}

void IO::set_filename ( string _s )
{
	is_new_file = true ;
	output_filename = _s ;
}

void IO::set_output_type ( int _t )
{
	output_type = _t ;
}

void IO::set_dimension ( int _d )
{
	dimension = _d ;
}

void IO::set_domain ( Domain * _d )
{
	domain_ptr = _d ;
	dimension = _d->dimension ;
}

IO & IO::operator<< ( double * v )
{
	output_data.push_back (v) ;
	return *this ;
}

IO & IO::operator<< ( const string &c )
{
	output_fields.push_back(c) ;
	return *this ;
}

IO & IO::operator<< ( const Scalar &s )
{
	if ( domain_ptr == NULL )
	{
		domain_ptr = s.domain_ptr.get() ;
		dimension = domain_ptr->dimension  ;
	} else if ( domain_ptr != s.domain_ptr.get() )
	{
		Log().TagDump ( logERROR ) << "Inconsistent assigned domain in IO object.";
	}

	output_data.push_back( s.data.get() ) ;
	output_fields.push_back( s.name ) ;
	return *this ;
}

// endl;
IO& IO::operator<< ( ostream & (*f)( ostream& ) )
{
	ostringstream ss;
	string _zone_name ;

	flush() ;
	reset() ;
	return *this ;
}

void IO::reset()
{
	domain_ptr = NULL ;
	is_new_file = true ;
	output_fields.clear() ;
	output_data.clear() ;
	output_filename.clear() ;
}

void IO::flush() 
{
	flush ( "ZONE" ) ;
}

void IO::flush( string zone_name )
{
	int i, field_number;
	int mpi_size, mpi_rank ;
	boost::shared_array<string> output_field_name ;
	boost::shared_array<double *> output_field_data ;

	mpi_size = domain_ptr->comm_size ;
	mpi_rank = domain_ptr->comm_rank ;

	if ( mpi_rank == 0 && (  output_fields.size() !=  output_data.size() ) )
		Log().TagDump( logERROR ) << "Output field captions and data are not matched!" << endl ;

	field_number = output_fields.size() ;
	output_field_name =  boost::shared_array<string> ( new string [field_number] );
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_name[i] = output_fields[i] ;
	};
	output_field_data =  boost::shared_array<double *>  ( new double * [ field_number ] ) ;
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_data[i] = output_data[i] ;
	};

	if ( output_type == IO_TECPLOT )
	{
		tecplot_output ( zone_name, field_number, output_field_name.get() , output_field_data.get() , output_filename , domain_ptr , is_new_file ) ;
	}
	else if ( output_type == IO_VTK )
	{
		vtk_output ( zone_name, field_number, output_field_name.get() , output_field_data.get() , output_filename , domain_ptr , is_new_file ) ;

	} else if ( output_type == IO_TECPLOT1D ){

		tecplot_output1D ( zone_name, field_number, output_field_name.get() , output_field_data.get() , output_filename , domain_ptr , is_new_file ) ;

	}

	if ( is_new_file ) is_new_file = false ;

}


void IO::flush_domain_info() 
{
	flush_domain_info ( domain_ptr->meshfile + '_' + domain_ptr->ID + ".dat" ) ;
}

void IO::flush_domain_info( string info_filename )
{
	int i, field_number;
	int mpi_size, mpi_rank ;
	boost::shared_array<string> output_field_name ;
	boost::shared_array<double *> output_field_data ;
	string zone_name ;

	mpi_size = domain_ptr->comm_size ;
	mpi_rank = domain_ptr->comm_rank ;

	if ( mpi_rank == 0 && (  output_fields.size() !=  output_data.size() ) )
		Log().TagDump( logERROR ) << "Output field captions and data are not matched!" << endl ;

	field_number = 3  ;
	is_new_file = true ;
	zone_name = "Domain Info" ;
	output_field_name =  boost::shared_array<string> ( new string [field_number] );
	output_field_data =  boost::shared_array<double *>  ( new double * [ field_number ] ) ;

	output_field_name[0] = "MPI_ID" ;
	output_field_name[1] = "LOCAL_CELL_ID" ;
	output_field_name[2] = "GLOBAL_CELL_ID" ;

	output_field_data[0] = new double [ domain_ptr->local_cell_number ] ;
	output_field_data[1] = new double [ domain_ptr->local_cell_number ] ;
	output_field_data[2] = new double [ domain_ptr->local_cell_number ] ;
	
	for ( i = 0 ; i < domain_ptr->local_cell_number ; i++ )
	{
		output_field_data[0][i] = 1. * domain_ptr->cell[i].mpi_id ;
		output_field_data[1][i] = 1. * domain_ptr->cell[i].local_id ;
		output_field_data[2][i] = 1. * domain_ptr->cell[i].id ;
	}

	if ( output_type == IO_TECPLOT )
	{
		tecplot_output ( zone_name, field_number, output_field_name.get() , output_field_data.get() , info_filename , domain_ptr , is_new_file ) ;
	}
	else if ( output_type == IO_VTK )
	{
		vtk_output ( zone_name, field_number, output_field_name.get() , output_field_data.get() , info_filename , domain_ptr , is_new_file ) ;
	}
	else if ( output_type == IO_TECPLOT1D )
	{
		cout<<"Error in in.cpp"<<endl; exit(1) ;
	}

	delete [] output_field_data[0]  ;
	delete [] output_field_data[1]  ;
	delete [] output_field_data[2]  ;
}


void IO::root_flush() 
{
	root_flush ( "ZONE" ) ;
}

void IO::root_flush( string zone_name )
{
	int i, field_number;
	int mpi_size, mpi_rank ;
	boost::shared_array<string> output_field_name ;
	boost::shared_array<double *> output_field_data ;

	mpi_size = domain_ptr->comm_size ;
	mpi_rank = domain_ptr->comm_rank ;

	if ( mpi_rank == 0 && (  output_fields.size() !=  output_data.size() ) )
		Log().TagDump( logERROR ) << "Output field captions and data are not matched!" << endl ;

	field_number = output_fields.size() ;
	output_field_name =  boost::shared_array<string> ( new string [field_number] );
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_name[i] = output_fields[i] ;
	};
	output_field_data =  boost::shared_array<double *>  ( new double * [ field_number ] ) ;
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_data[i] = output_data[i] ;
	};

	tecplot_output_root ( zone_name, field_number, output_field_name.get() , output_field_data.get() , output_filename , domain_ptr , is_new_file ) ;
	//vtk_output_root ( zone_name, field_number, output_field_name.get() , output_field_data.get() , output_filename , domain_ptr , is_new_file ) ;

	if ( is_new_file ) is_new_file = false ;

}

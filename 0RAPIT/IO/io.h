#include <fstream>
#include "scalar.h"
#include "domain.h"

using namespace std;

#ifndef IO_H
#define IO_H

#define IO_TECPLOT		0 
#define IO_VTK 			1
#define IO_TECPLOT1D	2 

class IO
{
	public:
		IO () ;
		IO ( string ) ;
		~IO () ;

		Domain *domain_ptr ;
		bool is_new_file ;
		int output_type, dimension;
		string output_filename ;
		vector<string>   output_fields ;
		vector<double *> output_data ;

		IO & operator << ( ostream & (*f)( ostream& ) ) ;
		IO & operator << ( double * ) ;
		IO & operator << ( const string & ) ;
		IO & operator << ( const Scalar & ) ;

		void set_filename ( string ) ;
		void set_domain ( Domain * ) ;
		void set_dimension ( int ) ;
		void set_output_type ( int ) ;
		void reset() ;
		void flush ( string ) ;
		void flush () ;
		void flush_domain_info (  ) ;
		void flush_domain_info ( string ) ;
		void root_flush ( string ) ;
		void root_flush () ;

};

#endif

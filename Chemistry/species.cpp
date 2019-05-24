#include "species.h"
#include "reaction.h"
#include "chemistry.h"
#include "main.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <libconfig.h++>

using namespace std ;
using namespace libconfig;

void const_diffusivity ( Species *s, int N )
{
	int i;
	for ( i = 0 ; i < N ; i++ )
	{
		s->diffusivity[i] = s->diffusivity_table [0] ;
	} 
}

void const_diffusivity_N ( Species *s, int N )
{
	int i;
	for ( i = 0 ; i < N ; i++ )
	{
		s->diffusivity[i] = s->diffusivity_table [0] / s->background_density [i] ;
	} 
}

void const_mobility ( Species *s, int N )
{
	int i;
	for ( i = 0 ; i < N ; i++ )
	{
		s->mobility[i] = s->mobility_table [0] ;
	} 
}

void const_mobility_N ( Species *s, int N )
{
	int i;
	for ( i = 0 ; i < N ; i++ )
	{
		s->mobility[i] = s->mobility_table [0] / s->background_density [i] ;
	} 
}


/*!
	The user defind function returning diffusivity of a given species from the table. The table is read and initialed when in Species::init(). 
*/
void table_diffusivity ( Species *s, int N )
{
	int i;
	double v;
	for ( i = 0 ; i < N ; i++ )
	{
		v = s->diffusivity_dependent_variable[i] ;

		s->diffusivity[i] = s->diffusivity_table [ (int) ( v / s->diffusivity_interval ) ] ;
	} 
}

/*!
	The user defind function returning diffusivity of a given species from the table. The table is read and initialed when in Species::init(). 
*/
void table_diffusivity_N ( Species *s, int N )
{
	int i;
	double v;
	for ( i = 0 ; i < N ; i++ )
	{
		v = s->diffusivity_dependent_variable[i] ;
		
		s->diffusivity[i] = s->diffusivity_table [ (int) ( v / s->diffusivity_interval ) ]  / s->background_density [i]  ;
	} 
}


/*!
	The user defind function returning mobility of a given species from the table. The table is read and initialed when in Species::init(). 
*/
void table_mobility (  Species *s, int N  )
{
	int i;
	double v ;
	for ( i = 0 ; i < N ; i++ )
	{
		v = s->mobility_dependent_variable[i] ;
		s->mobility[i] = s->mobility_table [ (int) ( v / s->mobility_interval ) ] ;
	} 
}

/*!
	The user defind function returning mobility of a given species from the table. The table is read and initialed when in Species::init(). 
*/
void table_mobility_N (  Species *s, int N  )
{
	int i ;
	double v ;
	for ( i = 0 ; i < N ; i++ )
	{

		v =s->mobility_dependent_variable[i] ;

		s->mobility[i] = s->mobility_table [ (int) ( v / s->mobility_interval ) ] / s->background_density[i] ;
	} 
}

void function_mobility (  Species *s, int N  )
{
	int i ;
	double v ;
	double pressure ;
	double EoP ;

	pressure = 0.5; //Torr
	for ( i = 0 ; i < N ; i++ )
	{

		v = s->mobility_dependent_variable[i] ;


		EoP = v / 100. / pressure ; 

		if ( EoP <= 60 )
			s->mobility[i] = 0.1 / pressure * ( 1. - 2.22e-3 * EoP ) ; 
		else
			s->mobility[i] = 8.25e-1 / pressure / sqrt( EoP ) * ( 1. - 86.52 / pow ( EoP , 1.5 ) ) ; 

		//cout << i << " e "<< v << "\tEoP " << EoP << "\tmobility " << s->mobility[i] <<  endl ;
	} 
}


void Species::init ( string filename )
{
	stringstream	ss_buffer ;
	string			line, tmp1, tmp2 ;
	fstream			file ;
	int				i ;

	name = filename.substr( 0, filename.find('.') ) ;
	file.open( filename.c_str(), ios::in ) ;
	while ( getline( file, line ) )
	{
		line = analysor ( line ) ;
		if ( line.size() > 0 && line[0] != '#' )
		{	
			ss_buffer.str(""); 
		   	ss_buffer.clear();
			ss_buffer << line.substr( 0, line.find('=') ) ;
			ss_buffer >> tmp1 ;
			if ( tmp1 == "NAME" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				//cout << line << " "<< official_name <<endl ;
				ss_buffer >> official_name ;
				//cout << official_name << endl ;
			} else if ( tmp1 == "DIFFUSIVITY" )
			{
				f_diffusivity_N  = false ;
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> tmp2 ;
				//if ( diffusivity_init ( tmp2, 0.01 ) != FILE_EXIST  )
				if ( diffusivity_init ( tmp2, 0.01 ) != 0  )
				{
					
					diffusivity_function = &const_diffusivity ; 
				} else 
				{
					diffusivity_function = &table_diffusivity ; 
				} 

			} else if ( tmp1 == "DIFFUSIVITY_N" )
			{
				f_diffusivity_N  = true ;
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> tmp2 ;
				//if ( diffusivity_init ( tmp2, 0.01 ) != FILE_EXIST  )
				if ( diffusivity_init ( tmp2, 0.01 ) != 0 )
				{
					
					diffusivity_function = &const_diffusivity_N ; 
				} else 
				{
					diffusivity_function = &table_diffusivity_N ; 
					//for ( i = 0 ; i < diffusivity_data_number ; i++ )
					//	cout << "DD " << diffusivity_table [i] << endl ;
				} 

 
			}  else if ( tmp1 == "MOBILITY" )
			{
				f_mobility_N  = false ;
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> tmp2 ;
				//if ( mobility_init ( tmp2, 0.01 ) != FILE_EXIST )
				if ( tmp2 == "DriftFunction" )
				{
					mobility_function = &function_mobility ;
				} else if ( mobility_init ( tmp2, 0.01 ) != 0 )
				{
					mobility_function = &const_mobility ; 
					//cout << mobility_table [0] << endl ;
				} else
				{
					mobility_function = &table_mobility  ;
					//for ( i = 0 ; i < mobility_data_number ; i++ )
					//	cout << mobility_table [i] << endl ;
				} 
				

			}  else if ( tmp1 == "MOBILITY_N" )
			{
				f_mobility_N  = true ;
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> tmp2 ;
				//if ( mobility_init ( tmp2, 0.01 ) != FILE_EXIST )
				if ( mobility_init ( tmp2, 0.01 ) != 0 )
				{
					mobility_function = &const_mobility_N ; 
				} else
				{
					mobility_function = &table_mobility_N  ;
					//for ( i = 0 ; i < mobility_data_number ; i++ )
					//	cout << mobility_table [i] << endl ;
				} 
			} else if ( tmp1 == "DEPENDENT_VARIABLE" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> dependent_variable_name  ;
				cout << "DepVar = " << dependent_variable_name << endl ;
			} else if ( tmp1 == "MASS_AMU" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> mass_amu ;
			} else if ( tmp1 == "MASS_KG" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> mass_kg ;
			} else if ( tmp1 == "VISCOSITY" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> viscosity ;
			} else if ( tmp1 == "POLARIZABILITY" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> polarizability ;
			} else if ( tmp1 == "LJ_POTENTIAL" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> LJ_potential ;
			} else if ( tmp1 == "BINARY_DIAMETER" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> binary_diameter ;
			} else if ( tmp1 == "CHARGE" )
			{
				ss_buffer.str(""); 
		   		ss_buffer.clear();
				ss_buffer << line.substr( line.find('=') + 1 ) ;
				ss_buffer >> charge ;

			}

		}
	}
	file.close() ;
}

int Species::mobility_init ( string filename , double interval )
{
	int				i , j ;
	stringstream	ss_buffer ;
	fstream			file ;
	string			line ;
	int				input_table_number  ;
	double			*dummy_dependent_variable ;
	double			*dummy_mobility ;

	input_table_number = 0 ;
	mobility_interval = interval ;
	file.open( filename.c_str(), ios::in ) ;
	if ( !file.is_open() ) 
	{
		mobility_data_number = 1 ;
		mobility_table = new double [1] ;
		mobility_table[0] = atof ( filename.c_str() ) ;

		return 1 ;
	}
	while ( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#' )
		{	
				input_table_number++ ;
		}
	}
	file.close() ;

	//cout<<"table_number= "<< input_table_number << endl ;

	dummy_dependent_variable = new double [ input_table_number ] ;
	dummy_mobility = new double [ input_table_number ] ;

	i = 0 ;
	file.open( filename.c_str(), ios::in ) ;
	while( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#')  
		{
			ss_buffer.str("") ;
			ss_buffer.clear() ;
			ss_buffer << line ;
			ss_buffer >> dummy_dependent_variable[i] >> dummy_mobility[i] ;
			i++ ;
		}
	}
	file.close() ;

	//mobility_data_number  = 1 + int ( ( dummy_dependent_variable[input_table_number-1] - dummy_dependent_variable[0] )/ interval ) ;
	mobility_data_number  = 1 + int ( ( 20. - 0.) / interval ) ;
		
	mobility_table = new double [ mobility_data_number ] ;
	mobility_dependent_variable_table= new double [ mobility_data_number ] ;

	for( i = 0 ; i < mobility_data_number  ; i++ )
	{ 
		//mobility_dependent_variable_table[i] = dummy_dependent_variable[0] + interval * i ;
		mobility_dependent_variable_table[i] = interval * i ;
		
		if ( mobility_dependent_variable_table[i] < dummy_dependent_variable[0] )
			mobility_table[i] = dummy_mobility [0] ;
		else if ( mobility_dependent_variable_table[i] >= dummy_dependent_variable[input_table_number -  1 ] )
			mobility_table[i] = dummy_mobility [ input_table_number -  1 ] ;
		else 
			for( j = 0 ; j < input_table_number-1 ; j++ )
			{
				if (  mobility_dependent_variable_table[i] >= dummy_dependent_variable[j] && mobility_dependent_variable_table[i] < dummy_dependent_variable[j+1] )
					mobility_table[i] = dummy_mobility [j] + ( dummy_mobility[j+1] - dummy_mobility [j] ) / ( dummy_dependent_variable[j+1] - dummy_dependent_variable[j] )  * ( mobility_dependent_variable_table[i] - dummy_dependent_variable[j] ) ;
				//cout << "AAA " <<  dummy_mobility [j]  << endl ;
			}

//		cout << "inin " << i << " " << mobility_dependent_variable_table[i] << " " << mobility_table[i] << endl;
	}


	return 0 ;
}


int Species::diffusivity_init ( string filename , double interval)
{
	int				i , j ;
	stringstream	ss_buffer ;
	fstream			file ;
	string			line ;
	int				input_table_number  ;
	double			*dummy_dependent_variable ;
	double			*dummy_diffusivity ;

	input_table_number = 0 ;
	diffusivity_interval = interval ;
	file.open( filename.c_str(), ios::in ) ;
	if ( !file.is_open() ) 
	{
		diffusivity_data_number = 1 ;
		diffusivity_table = new double [1] ;
		diffusivity_table[0] = atof ( filename.c_str() ) ;

		return 1 ;
	}
	while ( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#' )
		{	
				input_table_number++ ;
		}
	}
	file.close() ;

	//cout<<"table_number= "<< input_table_number << endl ;

	dummy_dependent_variable = new double [ input_table_number ] ;
	dummy_diffusivity = new double [ input_table_number ] ;

	i = 0 ;
	file.open( filename.c_str(), ios::in ) ;
	while( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#')  
		{
			//file_I.seekg(-line.size()-1,ios::cur ) ;
			ss_buffer.str("") ;
			ss_buffer.clear() ;
			ss_buffer << line ;
			ss_buffer >> dummy_dependent_variable[i] >> dummy_diffusivity [i] ;
			//file_I >> dummy[mm].Electron_temper >> dummy[mm].mobility ;
			i++ ;
		}
	}
	file.close() ;

	//diffusivity_data_number  = 1 + int ( ( dummy_dependent_variable[input_table_number-1] - dummy_dependent_variable[0] )/ interval ) ;
	diffusivity_data_number  = 1 + int ( ( 20. - 0.) / interval ) ;
		
	diffusivity_table = new double [ diffusivity_data_number ] ;
	diffusivity_dependent_variable_table= new double [ diffusivity_data_number ] ;

	for( i = 0 ; i < diffusivity_data_number  ; i++ )
	{ 
		//diffusivity_dependent_variable_table[i] = dummy_dependent_variable[0] + interval * i ;
		diffusivity_dependent_variable_table[i] = interval * i ;


		if ( diffusivity_dependent_variable_table[i] < dummy_dependent_variable[0] )
			diffusivity_table[i] = dummy_diffusivity[0] ;
		else if ( diffusivity_dependent_variable_table[i] >= dummy_dependent_variable[input_table_number -  1 ] )
			diffusivity_table[i] = dummy_diffusivity [ input_table_number -  1 ] ;
		else 
			for( j = 0 ; j < input_table_number-1 ; j++ )
			{
				if (  diffusivity_dependent_variable_table[i] >= dummy_dependent_variable[j] && diffusivity_dependent_variable_table[i] < dummy_dependent_variable[j+1] )
					diffusivity_table[i] = dummy_diffusivity [j] + ( dummy_diffusivity [j+1] - dummy_diffusivity [j] ) / ( dummy_dependent_variable[j+1] - dummy_dependent_variable[j] )  * ( diffusivity_dependent_variable_table[i] - dummy_dependent_variable[j] ) ;
			}

	}
	return 0;
}

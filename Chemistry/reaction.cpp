#include <string>
#include <cstdlib>
#include <iostream>
#include "species.h"
#include "reaction.h"
#include "chemistry.h"
#include <sstream>
#include <fstream>
#include "main.h"
#include <libconfig.h++>

using namespace std ;
using namespace libconfig;

double const_reaction_rate_coefficient( Reaction *r , int i )
{	
	return r->rate_coefficient_table [0];
}


double table_reaction_rate_coefficient( Reaction *r , int i )
{
	double v;

	v = r->dependent_variable[i] ;

	r->rate_coefficient = r->rate_coefficient_table [int ( v / r->dependent_variable_interval ) ] ;

	return r->rate_coefficient ;
}

Reaction::Reaction()
{
	reactant_number 			= 0 ;
	product_number 			= 0 ;
	dependent_variable_max	= 20. ;
	dependent_variable_min	= 0. ;
	rate_data_number			= 10001. ;
	threshold_energy			= 0. ;
}

/*! \brief Analysis the chemistry equation.

	This function will analysis a given equation and find the reactants and products in the equation. The coefficient infront of the species will also be read. At this moment, this function can not analysis the balance of atoms in the equation.

*/
void Reaction::equation_init( string line )
{
	int					i, j, counter ;
	stringstream		ss_buffer ;
	string				reactants, products ;
	string				tmp1, tmp2 ;	
	string::iterator		iter;
	size_t				start, end; 
	
	ss_buffer.str(""); 
	ss_buffer.clear();

	ss_buffer << line.substr( 0, line.find("->") ) ;
	ss_buffer >> reactants;
	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << line.substr( line.find("->") +2 ) ; 
	ss_buffer >> products;


	counter = 1 ;
	for( iter = reactants.begin() ; iter != reactants.end() ; iter++ ) 
	{
		if(  *iter  == '+')
		{
			counter++ ;
		}	
	}	

	reactant_number = counter ;

	counter = 1 ;
	for( iter = products.begin() ; iter != products.end() ;iter++ ) 
	{
		if(  *iter  == '+')
		{
			counter++ ;
		}	
	}	
	product_number = counter ;

	reactant_name  = new string [ reactant_number ] ;
	product_name  = new string [ product_number ] ;
	reactant_coefficient = new double [ reactant_number ] ;
	product_coefficient = new double [ product_number ] ;


	ss_buffer.str(""); 
	ss_buffer.clear();
	//ss_buffer << line.substr( 0, line.find('>')-1 ) ;
	ss_buffer << line.substr( 0, line.find("->") ) ;
	ss_buffer >> reactants ;

	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << line.substr( line.find("->")+2 ) ;
	ss_buffer >> products ;

	/*  reactants */
	j = 0 ;
	start = reactants.find_first_not_of('+')	;
	end = 0 ;
	while( start != string::npos )
	{
		end = reactants.find_first_of( '+' , start+1 ) ;
		if( end == string::npos )
			end = reactants.length() ;
				
		tmp1 = reactants.substr ( start , end-start )	;
				
		if ( tmp1.find_first_of("0123456789") == 0 )
		{
			reactant_name[j] = tmp1.substr( 1 , end - start ) ;
			//species_names.push_back ( reaction[i].reactant_name[j] ) ;
			ss_buffer.str(""); 
			ss_buffer.clear();
			ss_buffer << tmp1[0] ;
			ss_buffer >> reactant_coefficient[j] ;
		} else 
		{
			reactant_name[j] = tmp1.substr( 0 , end - start ) ;
			//species_names.push_back ( reaction[i].reactant_name[j] ) ;
			reactant_coefficient[j] = 1. ;
		}
		start = reactants.find_first_not_of( '+' , end+1 ) ;
		j++ ;
	}

	/* products */
	j = 0 ;
	start = products.find_first_not_of('+')	;
	end = 0 ;
	while( start != string::npos )
	{
		end = products.find_first_of( '+' , start+1 ) ;
		if( end== string::npos )
			end = products.length() ;

		tmp1 = products.substr ( start , end-start )	;
				
		if ( tmp1.find_first_of("0123456789") == 0 )
		{
			product_name[j] = tmp1.substr( 1 , end-start ) ;
			//species_names.push_back ( reaction[i].reactant_name[j] ) ;
			ss_buffer.str(""); 
			ss_buffer.clear();
			ss_buffer << tmp1[0] ;
			ss_buffer >> product_coefficient[j] ;
		} else
		{
			product_name[j] = tmp1.substr( 0, end - start ) ;
			//species_names.push_back ( reaction[i].product_name[j] ) ;
			product_coefficient[j] = 1. ;

		}
		start = products.find_first_not_of( '+' , end+1 ) ;	
		j++ ;
	}


	reactant = new double * [ reactant_number ] ; 
	product = new double * [ product_number ] ; 

}



void Reaction::init( string line )
{
	int					i, j, counter ;
	stringstream		ss_buffer ;
	string				reactants, products ;
	string				tmp1, tmp2 ;	
	string::iterator		iter;
	size_t				start, end; 
	
	ss_buffer.str(""); 
	ss_buffer.clear();

	ss_buffer << line.substr( 0, line.find("->") ) ;
	ss_buffer >> reactants;
	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << line.substr( line.find("->") +2 ) ; 
	ss_buffer >> tmp1 ;

	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << tmp1.substr( 0, tmp1.find(';') ) ;
	ss_buffer >> products ;

	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << tmp1.substr( tmp1.find(';')+1 ) ;
	ss_buffer >> tmp2 ;

	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << tmp2.substr( 0 , tmp2.find(':') ) ;
	ss_buffer >> type ;
				
	if ( type == "table" ) 
	{
		ss_buffer.str(""); 
		ss_buffer.clear();
		ss_buffer << tmp2.substr(  tmp2.find(':') + 1 ) ;
		ss_buffer >> rate_table_filename;

		cout << rate_table_filename << endl ;

		rate_table_init  ( 0.01  ) ;
		
		/* Assign rate calculating function */
		rate_function = &table_reaction_rate_coefficient  ;

		/* Need to complete this part since not all the dependent variables are Te. */
		dependent_variable_name = "Te" ;
	}else if ( type == "const" )
	{
		ss_buffer.str(""); 
		ss_buffer.clear();
		ss_buffer << tmp2.substr(  tmp2.find(':') + 1 ) ;
		rate_coefficient_table =  new double [1] ;
		ss_buffer >> rate_coefficient_table [0] ;
		rate_table_filename  = "None" ;

		/* Assign rate calculating function */
		rate_function = &const_reaction_rate_coefficient;

	} else
	{
		/* Assign rate calculating function */
		rate_function = NULL ;
	}
	
	counter = 1 ;
	for( iter = reactants.begin() ; iter != reactants.end() ; iter++ ) 
	{
		if(  *iter  == '+')
		{
			counter++ ;
		}	
	}	

	reactant_number = counter ;

	counter = 1 ;
	for( iter = products.begin() ; iter != products.end() ;iter++ ) 
	{
		if(  *iter  == '+')
		{
			counter++ ;
		}	
	}	
	product_number = counter ;

	reactant_name  = new string [ reactant_number ] ;
	product_name  = new string [ product_number ] ;
	reactant_coefficient = new double [ reactant_number ] ;
	product_coefficient = new double [ product_number ] ;


	ss_buffer.str(""); 
	ss_buffer.clear();
	//ss_buffer << line.substr( 0, line.find('>')-1 ) ;
	ss_buffer << line.substr( 0, line.find("->") ) ;
	ss_buffer >> reactants ;
	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << line.substr( line.find("->")+2 ) ;
	ss_buffer >> tmp1 ;

	ss_buffer.str(""); 
	ss_buffer.clear();
	ss_buffer << tmp1.substr( 0, tmp1.find(';') ) ;
	ss_buffer >> products ;

	/*  reactants */
	j = 0 ;
	start = reactants.find_first_not_of('+')	;
	end = 0 ;
	while( start != string::npos )
	{
		end = reactants.find_first_of( '+' , start+1 ) ;
		if( end == string::npos )
			end = reactants.length() ;
				
		tmp1 = reactants.substr ( start , end-start )	;
				
		if ( tmp1.find_first_of("0123456789") == 0 )
		{
			reactant_name[j] = tmp1.substr( 1 , end - start ) ;
			//species_names.push_back ( reaction[i].reactant_name[j] ) ;
			ss_buffer.str(""); 
			ss_buffer.clear();
			ss_buffer << tmp1[0] ;
			ss_buffer >> reactant_coefficient[j] ;
		} else 
		{
			reactant_name[j] = tmp1.substr( 0 , end - start ) ;
			//species_names.push_back ( reaction[i].reactant_name[j] ) ;
			reactant_coefficient[j] = 1. ;
		}
		start = reactants.find_first_not_of( '+' , end+1 ) ;
		j++ ;
	}

	/* products */
	j = 0 ;
	start = products.find_first_not_of('+')	;
	end = 0 ;
	while( start != string::npos )
	{
		end = products.find_first_of( '+' , start+1 ) ;
		if( end== string::npos )
			end = products.length() ;

		tmp1 = products.substr ( start , end-start )	;
				
		if ( tmp1.find_first_of("0123456789") == 0 )
		{
			product_name[j] = tmp1.substr( 1 , end-start ) ;
			//species_names.push_back ( reaction[i].reactant_name[j] ) ;
			ss_buffer.str(""); 
			ss_buffer.clear();
			ss_buffer << tmp1[0] ;
			ss_buffer >> product_coefficient[j] ;
		} else
		{
			product_name[j] = tmp1.substr( 0, end - start ) ;
			//species_names.push_back ( reaction[i].product_name[j] ) ;
			product_coefficient[j] = 1. ;

		}
		start = products.find_first_not_of( '+' , end+1 ) ;	
		j++ ;
	}


	reactant = new double * [ reactant_number ] ; 
	product = new double * [ product_number ] ; 

}

void Reaction::rate_table_init ( string filename, double interval )
{
	int				i , j ;
	stringstream 	ss_buffer ;
	fstream			file ;
	string			line, tmp1 ;
	double			*dummy_rate, *dummy_dependent_variable ;
	int 				input_table_number  ;

	dependent_variable_table	= new double [ rate_data_number ] ;
	rate_coefficient_table	= new double [ rate_data_number ] ;
	dependent_variable_interval = ( dependent_variable_max - dependent_variable_min ) / (rate_data_number - 1) ;

	rate_table_filename = filename;
	cout << "Reading rate coefficient from " << rate_table_filename <<"\n  Lower bound: " <<dependent_variable_min <<"\n  Upper bound: "  << dependent_variable_max << "\n Iterval: " <<  dependent_variable_interval  << endl; 

	input_table_number = 0 ;
	file.open( rate_table_filename.c_str(), ios::in ) ;
	while ( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#' )
		{	
			input_table_number++ ;
		}
	}
	file.close() ;

	dummy_rate = new double [ input_table_number ] ;
	dummy_dependent_variable = new double [ input_table_number ] ;

	i = 0;
	file.open( rate_table_filename.c_str(), ios::in ) ;
	while( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#' )
		{	
			ss_buffer.str("") ;
			ss_buffer.clear() ;
			ss_buffer << line ;
			ss_buffer >> dummy_dependent_variable[i] >> dummy_rate[i] ;
			i++ ;
		}
	}
	file.close() ;

	for( i = 0 ; i < rate_data_number ; i++ )
	{ 
		dependent_variable_table[i] = dependent_variable_interval * i ;
		
		if(  dependent_variable_table[i] <= dummy_dependent_variable[0]  )
			rate_coefficient_table[i] = dummy_rate [0] ;
		else if (  dependent_variable_table[i] > dummy_dependent_variable[ input_table_number-1 ]  )
			rate_coefficient_table[i] = dummy_rate [ input_table_number-1 ] ;
		else
		{
			for( j = 0 ; j < input_table_number-1 ; j++ )
			{
				if ( dependent_variable_table[i] >= dummy_dependent_variable[j] && dependent_variable_table[i] <dummy_dependent_variable[j+1] )
				{
					rate_coefficient_table [i] = dummy_rate [j] + ( dependent_variable_table[i] - dummy_dependent_variable[j] ) * ( dummy_rate[j+1] - dummy_rate[j] ) / ( dummy_dependent_variable[j+1] - dummy_dependent_variable[j] ) ;
				}
			}
		}
		//cout << dependent_variable_table[i] << " " << rate_coefficient_table [i] << endl ;
	}
}

void Reaction::rate_table_init ( double interval )
{
	int				i , j ;
	stringstream 	ss_buffer ;
	fstream			file ;
	string			line, tmp1 ;
	double			*dummy_rate, *dummy_dependent_variable ;
	int 				input_table_number  ;

	dependent_variable_table	= new double [ rate_data_number ] ;
	rate_coefficient_table	= new double [ rate_data_number ] ;
	dependent_variable_interval = ( dependent_variable_max - dependent_variable_min ) / (rate_data_number - 1) ;

	cout << "Reading rate coefficient from " << rate_table_filename <<"\n  Lower bound: " <<dependent_variable_min <<"\n  Upper bound: "  << dependent_variable_max << "\n Iterval: " <<  dependent_variable_interval  << endl; 

	input_table_number = 0 ;
	file.open( rate_table_filename.c_str(), ios::in ) ;
	while ( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#' )
		{	
			ss_buffer.str(""); 
			ss_buffer.clear();
			ss_buffer << line.substr( 0, line.find('=') ) ;
			ss_buffer >> tmp1 ;
			input_table_number++ ;
		}
	}
	file.close() ;
//	cout<<"table_number= "<< table_number << endl ;

	dummy_rate = new double [ input_table_number ] ;
	dummy_dependent_variable = new double [ input_table_number ] ;

	i = 0;
	file.open( rate_table_filename.c_str(), ios::in ) ;
	while( getline( file, line ) )
	{
		if ( line.size() > 0 && line[0] != '#' )
		{	
			ss_buffer.str("") ;
			ss_buffer.clear() ;
			ss_buffer << line ;
			ss_buffer >> dummy_dependent_variable[i] >> dummy_rate[i] ;
			i++ ;
		}
	}
	file.close() ;
	
	for( i = 0 ; i < rate_data_number ; i++ )
	{ 
		dependent_variable_table[i] = dependent_variable_interval * i ;
		
		if(  dependent_variable_table[i] <= dummy_dependent_variable[0]  )
				rate_coefficient_table[i] = dummy_rate [0] ;
		else if (  dependent_variable_table[i] > dummy_dependent_variable[ input_table_number-1 ]  )
				rate_coefficient_table[i] = dummy_rate [ input_table_number-1 ] ;
		else
		{
			for( j = 0 ; j < input_table_number-1 ; j++ )
			{
				if ( dependent_variable_table[i] >= dummy_dependent_variable[j] && dependent_variable_table[i] <dummy_dependent_variable[j+1] )
				{
					rate_coefficient_table [i] = dummy_rate [j] + ( dependent_variable_table[i] - dummy_dependent_variable[j] ) * ( dummy_rate[j+1] - dummy_rate[j] ) / ( dummy_dependent_variable[j+1] - dummy_dependent_variable[j] ) ;
				}
			}
		}
		
	}

}

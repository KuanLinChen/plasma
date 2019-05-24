#include <string>
#include <iostream>
#include <cmath>
using namespace std;

#ifndef __REACTION_H
#define __REACTION_H


class Reaction
{
	public:
		Reaction ( ) ;
		int		reactant_number ;
		int		product_number ;

		Species **reactant_species ;
		Species **product_species ;

		double	*reactant_coefficient ;
		double	*product_coefficient ;
		string  *reactant_name ;
		string  *product_name ;
		string  dependent_variable_name ;
		
		string	rate_table_filename , type ;
		double	dependent_variable_max, dependent_variable_min,  dependent_variable_interval ;
		int		rate_data_number ;
		double	*dependent_variable_table ;
		double	*rate_coefficient_table ;
		double	rate_coefficient ;

		double	threshold_energy ;

		/* pointer attached onto the memory */
		double	**reactant ;
		double	**product ;
		double 	*dependent_variable, key_dependent_variable ;
		double	reaction_rate ;
		double	energy_lost ;

		/* rate coefficient */
		double 	order , Aconst , Energy ;
		double	rate ;
		//double ( * rate_function ) ( Reaction * ) ;
		double ( * rate_function ) ( Reaction * , int ) ;


		void	init( string ) ;
		void	rate_table_init ( double  ) ;
		void	equation_init( string ) ;
		void	rate_table_init ( string,  double  ) ;
		void 	calculate_reaction_rate () 
		{
			int i ;
			rate_coefficient = rate_function( this, 0 ) ;
			reaction_rate = rate_coefficient ;
			for ( i = 0 ; i < reactant_number ; i++ )
				reaction_rate *= pow ( *(reactant[i]), reactant_coefficient[i] );
		}
		void 	calculate_reaction_rate ( int index )  
		{
			int i ;
			rate_coefficient = rate_function( this, index ) ;
			reaction_rate = rate_coefficient ;
			for ( i = 0 ; i < reactant_number ; i++ )
				reaction_rate *= pow ( *(reactant[i] + index ), reactant_coefficient[i] );
		}


	private:


};



extern double const_reaction_rate_coefficient( Reaction * , int  ) ;
extern double table_reaction_rate_coefficient( Reaction * , int  ) ;


#endif
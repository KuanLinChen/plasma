#include <string>
#include <boost/shared_array.hpp>
#include "species.h"
#include "reaction.h"

using namespace std;

#ifndef __CHEMISTRY_H
#define __CHEMISTRY_H

// Chemistry main
class Chemistry
{
	public:
		Chemistry( string ) ;
		Chemistry( ) ;
//		~Chemistry( ) ;

		void	init ( string ) ;
		int		reaction_number ;
		int		species_number ;
		int		system_number ; /* Or we can say: How many grids that need chemistry module */

		Reaction	*reaction ;
		Species		*species ;

		double		**reaction_rate_of_species ;
		double		*reaction_rate_of_equation ;
		double		*energy_lost_of_equation ;
		double		*energy_lost_of_electron ;

		boost::shared_array<int>			number_RRoS ; /* How many related reactions of species */
		double		**coefficient_RRoS ; /* The reaction coefficient of related reactions of species */
		double		**coefficient_energy_RRoS ; /* The energy lost coefficient of related reactions of species */
		double		***pointer_RRoS ; /* Pointer of related reactions of species, which will be pointed to the memory of reaction rate of chemistry reaction. */
		double		***pointer_energy_RRoS ; /* Pointer of related reactions of species, which will be pointed to the memory of reaction rate of chemistry reaction. */

		int *index_of_momentum_trans_reaction ;
		int number_of_momentum_trans_reaction ;
		double ** target_gas_of_momentum_trans_reaction ;
		double *coefficient_momentum_trans_reaction ;
		double * energy_lost_of_momentum_trans_reaction ;
		double * electron_energy_density ;

		void set_pointer_of_species ( string, double * ) ;
		void set_pointer_of_species_transport_dependent_variable ( string, double * ) ;
		void set_pointer_of_reaction_rate ( string, double * ) ;
		void set_pointer_of_reaction_rate_dependent_variable ( string, double* );
		void set_pointer_of_species_diffusivity ( string , double * ) ;
		void set_pointer_of_species_mobility ( string , double * ) ;
		void set_pointer_of_energy_lost ( double * ) ;
		void set_pointer_of_momentum_trans_energy_lost ( double * ) ;
		void set_pointer_of_electron_energy_density ( double * ) ;

		void set_background_density ( string , double  * ) ;
		int	get_electron_ID ();
		void solve_RRoS( ) ;
		void solve_transport ( int ) ;

	private:

};


#endif

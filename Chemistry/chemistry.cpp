#include <string>
#include <cstdlib>
#include <iostream>
#include "species.h"
#include "reaction.h"
#include "chemistry.h"
#include <sstream>
#include <fstream>
#include <vector>
#include "main.h"
#include <libconfig.h++>

using namespace std ;
using namespace libconfig;

/*!

@file chemistry.cpp chemistry module.

*/


/*!

\brief Constructor of the Chemistry class with argument

The argument is the input file of chemistry module.

*/

Chemistry::Chemistry( string filename )
{
	init ( filename ) ;
}

/*!

\brief Constructor of the Chemistry class without argument

Without the argument, default input file will be \it channel.txt.

*/

Chemistry::Chemistry(  )
{

	cout << "No specific reaction input file assigned. Use channel.txt as defult." << endl;
	init ( "channel.txt" ) ;
}


/*! \brief Initiallization of cemistry module

To initialize the chemistry, we have to know how many reactions and species in a set of chemistry reaction. For a given \b filename, the code will read the file line by line, and line spaces will be removed. if the first character is '#', this line will be treated as comment and ignored. The rest of non-empty lines will be counted as the number of reaction (\b reaction_number). Once we know the number of reaction, pointer of \b Reaction will be used to allocate a series of variables. The file will open again and analysis line by line using the member function Reaction::init(). The function will tell us how many reactants, products, and what species are involved.

When reaction initiallization has done, we will analysis the species information from reaction objects to find the species in this system. Since the species might appear in many reactions, the collecting of the species information from reactions will sort in order and purify by removing the same ones. Variable \ti species_number will record how many species type in this system and alocalte the a series of variable of \it species. The Species::init()

@param filename the input file of chemistry

*/

void Chemistry::init ( string filename )
{
	int					i , j , k, counter  ;
	bool				flag ;
	stringstream		ss_buffer ;
	string				line,  tmp1, tmp2 ;
	fstream				file ;
	vector<string>		species_names ;
	bool				m ;

	Config				reaction_cfg ;

	number_RRoS = boost::shared_array<int>( new int [10]) ;
	number_RRoS[3] = 1 ;
	cout << "tt " << number_RRoS[3]  << endl;

	try{
		reaction_cfg.readFile ( filename.c_str() );
		//reaction_cfg.readFile ( filename.c_str() );
	} catch (libconfig::ParseException& e)
	{
		cout << "Error reading Configuration file (" << filename << ")!" << endl;
		cout << "ParseException at Line " << e.getLine() << ": " << e.getError() << endl;
		exit(EXIT_FAILURE);
	} catch (libconfig::SettingException& e) {
		cout << "Error reading Configuration file (" << filename << ")!" << endl;
		cout << "SettingException at " << e.getPath() << endl;
		exit(EXIT_FAILURE);
	} catch (std::exception& e)
	{
		cout << "Error reading Configuration file (" << filename << ")!" << endl;
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	} catch (...) {
		cout << "Error reading Configuration file (" << filename << ")!" << endl;
		cout << "Unknown Exception" << endl;
		exit(EXIT_FAILURE);
	}

	try
	{
		string information = reaction_cfg.lookup( "information" );
		cout << "Chemistry: " << filename << " : " << information << endl << endl;
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'information' setting in configuration file." << endl;
	}

	Setting &root = reaction_cfg.getRoot() ;
	try
	{
    	const Setting &reactions = root["chemistry"]["reactions"];
    	int count = reactions.getLength();
		reaction_number = count;

		cout << reaction_number << " reaction(s) found." << "\n"  ;
		reaction = new Reaction [ reaction_number ] ;

		index_of_momentum_trans_reaction = new int [ reaction_number ] ;
		coefficient_momentum_trans_reaction = new double [ reaction_number ] ;
		target_gas_of_momentum_trans_reaction = new double * [ reaction_number ] ;
		number_of_momentum_trans_reaction = 0 ;

   		for( i = 0 ; i < reaction_number ; i++ )
   		{
   			const Setting &r = reactions[i];
			// Only output the record if all of the expected fields are present.
			string equation,rate_table;
			double threshold, rate_constant;

			threshold = 0. ;
			rate_constant = 0. ;

			// Analysis the given chemistry reaction equation.
			if ( r.lookupValue ( "equation", equation ) )
			{
				equation = analysor( equation ) ;
    				reaction[i].equation_init( equation ) ;
    		} else
    			cout << "Reaction " << i << " needs reaction equation. " ;

			// Read the threshold energy when input file gives. Otherwise, the threshold energy is zero.
			if ( ! ( r.lookupValue ( "threshold", threshold) ) )
				reaction[i].threshold_energy = 0. ;
			else
				reaction[i].threshold_energy = threshold  ;

			// Get reaction rate table or constant. User has to provide one from either rate table or rate constant.
			if ( r.lookupValue ( "rate_table", rate_table ) )
			{
				reaction[i].rate_table_init ( rate_table,  0.01 ) ;
				reaction[i].rate_function = &table_reaction_rate_coefficient  ;
				reaction[i].dependent_variable_name = "Te" ;
			} else if ( r.lookupValue ( "rate_constant", rate_constant )  )
			{
				reaction[i].rate_coefficient_table =  new double [1] ;
				reaction[i].rate_coefficient_table [0] = rate_constant ;
				reaction[i].rate_function = &const_reaction_rate_coefficient ;
			}
			else
				cout << "Reaction " << i << " needs reaction rate coefficient. " ;
			//cout << equation << endl ;

			// Find momentum transfer collision
			m = false ;
			index_of_momentum_trans_reaction [i] = -1 ;
			coefficient_momentum_trans_reaction[i] = 0 ;
			if ( r.lookupValue ( "momentum", m ) )
			{
				if ( m )
				{
					index_of_momentum_trans_reaction [ number_of_momentum_trans_reaction ] = i ;
					number_of_momentum_trans_reaction ++ ;

				}
			}
		}
		cout << endl;
	}
	catch ( const SettingNotFoundException &nfex )
	{

		//exit(EXIT_FAILURE);
	}


	/*!
		The names of reactant and product will put into the container ({\bf species_names}) from all the reactions. In order to eliminate the same species names from the list, we sort it in order and erase by standard std. By this way we can easily know how many species ({\bf species_number}) in the chemistry set.
	*/
	/* Initialization of species */
	for ( i = 0 ; i < reaction_number ; i++ )
	{
		for ( j = 0 ; j < reaction[i].reactant_number ; j ++ )
			species_names.push_back ( reaction[i].reactant_name[j] ) ;
		for ( j = 0 ; j < reaction[i].product_number; j ++ )
			species_names.push_back ( reaction[i].product_name[j] ) ;
	}

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		if ( i == mpi_id )
		{
			sort ( species_names.begin() , species_names.end()  ) ;
			species_names.erase( unique( species_names.begin(), species_names.end()), species_names.end()) ;
			species_number = species_names.size()  ;
		}
	}

	/*!
		Once we know the number of species, we have to initialize all the species information by reading the species information from file. The information file name is based on the species name. For example, a species named "e" will automatically read the information from file "e.txt".
	*/
	species = new Species [species_number ] ;
	for ( j = 0 ; j < mpi_size ; j++ )
	{
		if ( j == mpi_id )
		{
			for ( i = 0 ; i < species_number ; i++ )
			{
				ss_buffer.str("") ;
				ss_buffer.clear() ;
				ss_buffer << species_names[i] << ".txt" ;
				ss_buffer >> tmp1  ;
				species[i].init( tmp1 ) ;
			}
		}
	}
	/* End of initialization of species */

	/*!
		The variable {\bf reaction_rate_of_species} denotes to the source term \f$S\f$ in the continuity equation.
		\f[
			\frac{\partial n}{\partial t} + \nabla\cdot\Gamma = S
		\f].

		The variable {\bf reaction_rate_of_equation} denotes to the term \f$R\f$ of a reaction. For example, if we have an inoization reaction
		\f[
			e + Ar \rightarrow Ar^{+} + 2 e
		\f], the reaction rate of this equation is
		\f[
			R = k_{ionization} n_{e} n_{Ar}
		\f]

		The variable {\bf energy_lost_of_equation} denotes to the term \f$\epsilon R\f$ of a reaction. Using the same example, if we have an inoization reaction and the threshold energy of the argon ionization is 15.9 eV, the energy lost of this equation is
		\f[
			\epsilon R = q \times 15.9 \times R = 1.6 \times 10^{-19} \times 15.9 \times k_{ionization} n_{e} n_{Ar}
		\f], where the unit is in joule.

	*/
	/* Building the relation from each reaction to total rate of species */
	reaction_rate_of_species = new double * [ species_number ] ; // Storage pointers of the reaction rate of species, using for the source/sink term in continuity equation.
	reaction_rate_of_equation = new double [reaction_number ] ;
	energy_lost_of_equation = new double [reaction_number ] ;


	/*!
		For a large set of chemistry reactions, each species involves in only few reactions. To save computation time, we have to specific find the reaction information for each species for the further calculation of reaction rate. Introduce an example of chemistry set:

		\f[
			R1: e + Ar \rightarrow Ar^{+} + 2 e
		\f]
		\f[
			R2: e + Ar \rightarrow Ar_{meta} + e
		\f]
		\f[
			R3: e + Ar \rightarrow Ar_{r} + e
		\f]
		\f[
			R4: Ar_{r} \rightarrow Ar
		\f]

		The variables named with {\bf RRoS} (Reaction Rate of Species) are describing as follows:
		pointer_RRoS are the pointers point to the memory of reaction rate. For example, the reaction rate of electron is
		\f[
			R_{e} = c_{1} k_{1} n_{e} n_{Ar} + c_{2} k_{2} n_{e} n_{Ar} + c_{3} k_{3} n_{e} n_{Ar}
		\f]
		The pointer will point to a given memory from user to storage the rate of electron. {\bf pointer_energy_RRoS} is simililar one that point to a memory to storege
		\f[
			R_{energy} = \epsilon_{1} k_{1} n_{e} n_{Ar} + \epsilon_{1} k_{2} n_{e} n_{Ar} + \epsilon_{1} k_{3} n_{e} n_{Ar}
		\f].
		The difference in the above equation is \f$c_{i}\f$ ({\bfcoefficient_RRoS}). For instance, in the reaction \f$R_{2}\f$, there is no electron lost in the reaction and \f$c_{2}\f$ is zero. However, reaction \f$R_{2}\f$ causes electron energy lost, which the lost of energy is \f$\epsilon_{2}\f$

		Again, we use electron as an example. When we calculate the source and sink of electrons, we don't have to scan all the raction to see if electron involve. Only 3 reactions in the reaction set invole with electron, so we record the relation in variable {\bf number_RRoS}. In this case, {\bf number_RRoS} for electron is 3, for \f$Ar^{+}\f$ is 1, for \f$Ar_{r}\f$ is 2, and for \f$Ar_{meta}\f$ is 1.

	*/

	number_RRoS 				= boost::shared_array<int> ( new int [ species_number ] );
	coefficient_RRoS			= new double *  [ species_number ] ;
	coefficient_energy_RRoS		= new double *  [ species_number ] ;
	pointer_RRoS				= new double ** [ species_number ] ;
	pointer_energy_RRoS			= new double ** [ species_number ] ;

	for ( j = 0 ; j < species_number ; j ++ )
	{
		counter = 0 ;
		for ( i = 0 ; i < reaction_number ; i++ )
		{
			flag = true ;
			if ( flag )
			{
				for ( k = 0 ; k < reaction[i].reactant_number ;k ++ )
					if ( reaction[i].reactant_name[k] == species[j].name )
					{
						counter ++ ;
						flag = false ;
						k = reaction[i].reactant_number ;
					}
			}
			if ( flag )
			{
				for ( k = 0 ; k < reaction[i].product_number ; k ++ )
					if ( reaction[i].product_name[k] == species[j].name )
					{
						counter ++ ;
						k = reaction[i].product_number ;
					}
			}
		}
		number_RRoS[j] 				= counter ;
		//cout << "RRoS "  << species[j].name << " " << number_RRoS[j] 	 << endl ;
		coefficient_RRoS[j]			= new double [ number_RRoS[j] ];
		coefficient_energy_RRoS[j]	= new double [ number_RRoS[j] ];
		pointer_RRoS[j]				= new double *[ number_RRoS[j] ];
		pointer_energy_RRoS[j]		= new double *[ number_RRoS[j] ];
	}

	for ( j = 0 ; j < species_number ; j ++ )
		for ( i = 0 ; i < number_RRoS[j]  ; i++ )
		{
			coefficient_RRoS[j][i] = 0. ;
			coefficient_energy_RRoS[j][i] = 0. ;
		}

	for ( j = 0 ; j < species_number ; j ++ )
	{
		counter = 0 ;
		for ( i = 0 ; i < reaction_number ; i++ )
		{
			flag = true ;
			for ( k = 0 ; k < reaction[i].reactant_number ; k ++ )
				if ( reaction[i].reactant_name[k] == species[j].name )
				{
					coefficient_RRoS[j][counter] -= reaction[i].reactant_coefficient[k]	;
					if ( flag )
					{
						flag = false ;
					}
				}
			for ( k = 0 ; k < reaction[i].product_number ; k ++ )
				if ( reaction[i].product_name[k] == species[j].name )
				{
					coefficient_RRoS[j][counter] += reaction[i].product_coefficient[k] ;
					if ( flag )
					{
						flag = false ;
					}
				}
			if ( !flag )
			{
				pointer_RRoS [j][ counter ] = & (reaction_rate_of_equation[ i ] );
				counter ++ ;
			}
		}
	}

	//for ( j = 0 ; j < species_number ; j ++ )
	//{
		j = get_electron_ID () ;
		counter = 0 ;
		for ( i = 0 ; i < reaction_number ; i++ )
		{
			flag = true ;
			for ( k = 0 ; k < reaction[i].reactant_number ; k ++ )
				if ( reaction[i].reactant_name[k] == species[j].name )
				{
					coefficient_energy_RRoS[j][counter] = reaction[i].threshold_energy	;
					//cout << coefficient_energy_RRoS[j][counter] << " " << j << " "<< counter << " " << endl ;
					if ( flag )
					{
						flag = false ;
					}
				}
			if ( !flag )
			{
				pointer_energy_RRoS [j][ counter ] = & ( reaction_rate_of_equation[ i ]);
				counter ++ ;
			}
		}
	//}

	for ( i = 0 ; i < reaction_number ; i++ )
	{
		reaction[i].reactant_species = new Species * [ reaction[i].reactant_number ] ;
		reaction[i].product_species = new Species * [ reaction[i].product_number ] ;


		for ( j = 0 ; j < reaction[i].reactant_number ; j ++ )
		{
			for ( k = 0 ; k < species_number ; k ++ )
			{
				if ( reaction[i].reactant_name[j] == species[k].name )
					reaction[i].reactant_species[ j ] = & species[k] ;
			}
		}



		for ( j = 0 ; j < reaction[i].product_number ; j ++ )
		{
			for ( k = 0 ; k < species_number ; k ++ )
			{
				if ( reaction[i].product_name[j] == species[k].name )
					reaction[i].product_species [j] = & species[k] ;
			}

		}

	}

}

void Chemistry::set_pointer_of_species( string name, double *p )
{
	int		i, j ;

	for( i = 0 ; i < reaction_number ; i ++ )
	{
		for ( j = 0 ; j < reaction[i].reactant_number ; j ++ )
		{
			if ( name == reaction[i].reactant_name[j] )
			{
				reaction[i].reactant[j] = p ;
			}
		}
		for ( j = 0 ; j < reaction[i].product_number ; j ++ )
		{
			if ( name == reaction[i].product_name[j] )
			{
				reaction[i].product[j] = p ;
			}
		}
	}
}

void Chemistry::set_pointer_of_species_transport_dependent_variable ( string name, double *p )
{
	int		i, j ;

	for( i = 0 ; i < species_number ; i ++ )
	{
		if ( name == species[i].name )
		{
			species[i].diffusivity_dependent_variable = p ;
			species[i].mobility_dependent_variable = p ;
		}
	}
}

void Chemistry::set_pointer_of_reaction_rate_dependent_variable( string name, double *p )
{
	int		i, j ;

	for( i = 0 ; i < reaction_number ; i ++ )
	{
		if ( reaction[i].dependent_variable_name == name )
		{
			reaction[i].dependent_variable = p ;
			//cout << " AAAA " << reaction[i].dependent_variable[10] << endl;
		}
	}
}

void Chemistry::set_pointer_of_reaction_rate ( string name, double *p )
{
	int		i ;
	for( i = 0 ; i < species_number ; i ++ )
	{
		if ( name == species[i].name )
			reaction_rate_of_species[i] = p ;
	}

}

void Chemistry::set_pointer_of_species_diffusivity ( string name, double *p )
{
	int i ;
	for( i = 0 ; i < species_number ; i ++ )
	{
		if ( name == species[i].name )
			species[i].diffusivity = p ;
	}
}

void Chemistry::set_pointer_of_species_mobility ( string name, double  *p )
{
	int i ;
	for( i = 0 ; i < species_number ; i ++ )
	{
		if ( name == species[i].name )
			species[i].mobility = p ;
	}
}

void Chemistry::set_pointer_of_energy_lost ( double  *p )
{
	energy_lost_of_electron = p ;
}

void Chemistry::set_pointer_of_electron_energy_density ( double  *p )
{
	electron_energy_density = p ;
}

void Chemistry::set_pointer_of_momentum_trans_energy_lost ( double  *p )
{

	int i, j, k ;
	double Me, Mg ;

	for ( i = 0 ; i <  number_of_momentum_trans_reaction ; i++ )
	{
		Me = species[ get_electron_ID() ].mass_amu ;
		for ( j =  0 ; j < reaction[ index_of_momentum_trans_reaction[i]  ].reactant_number ; j++ )
			if ( reaction[ index_of_momentum_trans_reaction[i]  ].reactant_species[j]->name != "e" )
			{
				Mg = reaction[ index_of_momentum_trans_reaction[i]  ].reactant_species[j]->mass_amu ;
				target_gas_of_momentum_trans_reaction [ i ] =  reaction[ index_of_momentum_trans_reaction[i]  ].reactant[j] ;
				//cout << "ASAA  \n\n"<< target_gas_of_momentum_trans_reaction [ i ]  << " " <<  reaction[ index_of_momentum_trans_reaction[i]  ].reactant[j]  << endl;
				//for ( k = 0 ; k < system_number ; k++ )
				//{
				//	cout << reaction[ index_of_momentum_trans_reaction[i]  ].reactant[j] [k] << endl ;
				//	cout << reaction[ index_of_momentum_trans_reaction[i]  ].reactant[j] [k] << endl ;
				//}
			}

		coefficient_momentum_trans_reaction[i] = 2. * Me / Mg ;

	}

	energy_lost_of_momentum_trans_reaction = p ;
}

void Chemistry::set_background_density ( string name , double  *p )
{
	int i ;
	for( i = 0 ; i < species_number ; i ++ )
	{
		if ( name == species[i].name )
			species[i].set_background_density ( p ) ;
	}
}


int	Chemistry::get_electron_ID ()
{
	int i ;
	for ( i = 0 ; i < species_number ; i++ )
		if ( species[i].name == "e")
			return i ;
}

/*! \brief To obtain the reaction rates of species.

\f[
	S_{k} = k_{i}
\f]
*/

void Chemistry::solve_RRoS ()
{
	int	i, j, k, l;

	for ( j = 0 ; j < system_number ; j ++ )
	{
		for ( i = 0 ;  i < reaction_number ; i++ )
		{
			reaction[i].calculate_reaction_rate ( j ) ;
			reaction_rate_of_equation[i] = reaction[i].reaction_rate ;
			energy_lost_of_equation[i] = reaction[i].energy_lost  ;
			//cout << i << " R " << reaction_rate_of_equation[i] << endl ;
			/*!
				\f$ R = k \times n_1 n_2 \f$
			*/
		}

		for ( k = 0 ; k < species_number ; k ++ )
		{
			* ( reaction_rate_of_species[k] + j ) = 0. ;
			for ( l = 0 ; l < number_RRoS[k] ; l++ )
			{
				* ( reaction_rate_of_species[k] + j ) += coefficient_RRoS[k][l] * *(pointer_RRoS [k][l]) ;
			}
		}

		k = get_electron_ID() ;
		* ( energy_lost_of_electron + j ) = 0. ;
		for ( l = 0 ; l < number_RRoS[k] ; l++ )
		{
			* ( energy_lost_of_electron + j ) += coefficient_energy_RRoS[k][l] * *(pointer_RRoS [k][l]) ;

			//cout << " KK " << l << " " << k << " " <<   coefficient_energy_RRoS[k][l] << " " <<  *(pointer_RRoS [k][l]) << endl ;
		}
	}

	for ( i = 0 ; i <  number_of_momentum_trans_reaction ; i++ )
	{
		for (  j = 0 ; j < system_number ; j ++ )
		{
			reaction[i].calculate_reaction_rate ( j ) ;
			//cout << "DT " << reaction[i].dependent_variable [j] << endl ;
			energy_lost_of_momentum_trans_reaction[j] = coefficient_momentum_trans_reaction[i] * reaction[index_of_momentum_trans_reaction[i] ].rate_coefficient * target_gas_of_momentum_trans_reaction[ i ] [j] * electron_energy_density[j]  ;
		}
	}

}

void Chemistry::solve_transport ( int i )
{
	int j ;

	species[i].diffusivity_function (  & species[i] ,  system_number ) ;
	species[i].mobility_function (  & species[i] ,  system_number ) ;

}


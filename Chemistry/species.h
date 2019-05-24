#include <string>

using namespace std;

#ifndef __SPECIES_H
#define __SPECIES_H


// Species information
class Species
{
	public:
		string 	name;
		string 	official_name;

		/* Diffusivity */
		double 	*diffusivity ;
		double 	*diffusivity_dependent_variable ;
		/* Diffusivity table */
		int		diffusivity_data_number ;
		double 	*diffusivity_table ;
		double 	*diffusivity_dependent_variable_table ;
		double	diffusivity_interval ;
		int		diffusivity_init ( string, double );
		/* Mobility */
		double 	*mobility ;
		double 	*mobility_dependent_variable ;

		int		mobility_data_number ;
		double 	*mobility_table ;
		double 	*mobility_dependent_variable_table ;
		double	mobility_interval ;
		int		mobility_init ( string, double );
		string  dependent_variable_name ;

		void	init( string ) ;
		double	charge ;
		double 	mass_amu, mass_kg ;
		double 	self_diff_constant ;
		double 	self_mob_constant ;
		double 	binary_diameter ;
		double 	LJ_potential ;
		double 	polarizability ;
		double 	viscosity ;

		bool	f_diffusivity_N, f_mobility_N ;
		void ( * diffusivity_function ) ( Species *, int  ) ;
		void ( * mobility_function ) ( Species *, int ) ;
		
		void	set_background_density ( double *p ) { background_density = p ; }  
		double	*background_density ;
	private:
} ;


extern void const_diffusivity ( Species *, int ) ;
extern void const_diffusivity_N ( Species *, int ) ;
extern void const_mobility ( Species *, int ) ;
extern void const_mobility_N ( Species *, int ) ;
extern void function_mobility (  Species *, int N ) ;

#endif

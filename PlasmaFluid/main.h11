#include <string>
#include <cmath>

using namespace std;

#if !defined(__MAIN_H)
#define __MAIN_H

#if defined USEMPI
#elif defined HAVE_MPI 
	#define USEMPI
#else 
	#define NOUSEMPI
#endif


extern int gargc ;
extern char **gargv ;

extern string analysor( string ) ;

const double		kB 		= 8.617332478e-5 ; /* Boltzmann constant; unit in eV/K */
const double		q 		= 1.6e-19 ; /* Unit charge; unit in Q */
const double		vacuum_permittivity 	=  8.8541878176e-12 ; /* Unit in F/m */
const double 		PI = 4.0*atan(1.0) ;
extern	string	program_path ;
extern	string	inputfile_path ;

#endif


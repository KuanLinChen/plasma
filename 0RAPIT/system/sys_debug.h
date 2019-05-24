#include "iostream"
#include <iomanip>  

using namespace std;

#if !defined(__SYS_DEBUG_H)
#define __SYS_DEBUG_H


#ifdef DEBUG

	#define dout cout 
	#define debug_info cout <<  __FILE__ << "(" << __LINE__ << ")" << endl
	#define debug_printf DEBUG_printf 

#else

	#define dout 0 && cout
	#define debug_info 
	#define debug_printf 

#endif


void DEBUG_printf( const char *, ... ) ;

#endif
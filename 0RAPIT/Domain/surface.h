#include "domain.h"

#if !defined(__SURFACE_H)
#define __SURFACE_H

class Surface
{
	public:
		Domain *BaseDomain ;
		string ID ;
		static set<string> ID_bank ;

} ;

#endif
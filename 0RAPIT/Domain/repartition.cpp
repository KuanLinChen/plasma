
#include "domain.h"
#include "repartition.h"

using namespace std;


void repartition ( Domain *old_domain, Domain *new_domain, double *weight ) 
{
	string input_filename ;

	input_filename = old_domain->input_filename ;

	new_domain = new Domain ( input_filename ) ;
}

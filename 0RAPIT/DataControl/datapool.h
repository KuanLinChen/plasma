#include <vector>
#include <string>
#include <map>
//#include "domain.h"
#include "dataexchanger.h"
//#include "datamanager.h"

using namespace std ;

#if !defined(__DATAPOOL_H)
#define __DATAPOOL_H

class DataPool
{
	public:
		double * create_dependent( string, string * ) ;
		double * create_instantaneous ( string, int ) ;
		//double * request_memory ( string, Domain *, int ) ;
		//Scalar   request_memory ( string, Domain *, int ) ;
		void	 show() ;

		//DataManager * domain_trans_register (  Domain *,  Domain * ) ;

	private:
		static vector<double *>		dependent_variable ;
		static vector<string **>	dependent_variable_name ;
		static vector<string>		dependent_name ;
		static vector<int>		dependent_dof ; /* Degree of free */

		static vector<double *>		instantaneous_variable ;
		static vector<string>		instantaneous_name ;

		static vector<double **>	average_variable ;
		static vector<double **>	average_weight ;
		static vector<double>		average_name ;

		//static vector<Scalar *>		variables ;
		//static vector<string>		variables_name ;

		/* For data interpolation between meshes */
		//static vector<DataManager *>	vec_datamanager ;
		//static vector<string>		vec_datamanager_tags ;
} ;

extern DataPool datapool ;

#endif

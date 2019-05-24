//#include <string>
#include <iostream>
//#include <cmath>
//#include <algorithm>
//#include <set>
//#include <sstream>
//#include "main.h"
#include "datapool.h"
//#include "dataexchanger_functions.h"
using namespace std;

/*!

@file datapool.cpp Memory stroage and transfer control.

*/


vector<double *>	DataPool::dependent_variable ;
vector<string **>	DataPool::dependent_variable_name ;
vector<string>		DataPool::dependent_name ;
vector<int>			DataPool::dependent_dof ; /* Degree of free */

vector<double *>	DataPool::instantaneous_variable ;
vector<string>		DataPool::instantaneous_name ;

vector<double **>	DataPool::average_variable ;
vector<double **>	DataPool::average_weight ;
vector<double>		DataPool::average_name ;

//vector<DataManager *>	vec_datamanager ;
//vector<string>		vec_datamanager_tags ;

/*
DataManager * domain_trans_register (  Domain * source,  Domain * target )
{
	//if (  source == target ) return  ;

	vector<DataManager *>::iterator	i ;
	string tags = source->meshfile + target->meshfile ;

	//i = find ( vec_datamanager_tags.begin(), vec_datamanager_tags.end() , tags ) ;

	if ( std::find( vec_datamanager_tags.begin(), vec_datamanager_tags.end(), tags ) != vec_datamanager_tags.end()  )
	{
		cout << "Not Found" << endl ;
		vec_datamanager_tags.push_back ( tags ) ;
	} else
	{
		cout << "Found" << endl ;
	}

}
*/

double * DataPool::create_dependent( string a, string *b )
{
	double *c ;
	return c ;
}

double * DataPool::create_instantaneous( string name, int num )
{
	int		i ;
	double	*memory ;

	cout << "New version datapool is underconstrution!!!! ";
	cout << "Create " << name << " in datapool " ;
	memory = new double [ num ] ;

	instantaneous_name.push_back( name ) ;
	instantaneous_variable.push_back( memory ) ;

	cout << memory << endl ;

	for ( i = 0; i < num ; i++ )
		memory[ i ] = 0. ;

	return memory ;
}

/*
Scalar DataPool::request_memory( string name, Domain * domain_pointer , int D )
{
	int		i ;

	Scalar scalar ( domain_pointer , D ) ;

	cout << "Create " << name << " in datapool " ;

	cout << scalar.data << endl ;

	return scalar ;
}
*/

DataPool	datapool ;


#include "tools.h"
#include <algorithm>
#include  <string>

using namespace std;

string analysor( string s )
{
	s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() ) ;
	return s ;
}
double DotProduct(double *A, double *B ){
	return A[0]*B[0] + A[1]*B[1] ;
};
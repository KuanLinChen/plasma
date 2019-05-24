#include <boost/shared_ptr.hpp>

#if !defined(__CONNECTION_H)
#define __CONNECTION_H

void build_connection ( int source_number, boost::shared_array<double> position_x, boost::shared_array<double> position_y, boost::shared_array<double> position_z, boost::shared_array<int> connection_number,  boost::shared_array< vector<int>  > related_local_id, boost::shared_array< vector<int>  > related_global_id, boost::shared_array< vector<int> > related_mpi_id, boost::shared_array< vector<double> > weighing,  Domain *BaseDomain) ;

void weighing_inversed_distance ( int base_number, double *x, double *y, double *z, double *weight ) ;

#endif

#include <boost/fusion/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include "dataexchanger.h"
#include "surface.h"

#if !defined(__DATAMANAGER_H)
#define __DATAMANAGER_H

class DataManager
{
	public:
		DataManager () ;
		~DataManager () ;

		int			map_iterator ;
		vector< boost::shared_ptr<DataExchanger> >	vec_dataexchanger_ptr ;
		map<string, int>		dataexchanger_index ;

		void	remove_relation ( Domain * , Domain * , int ) ;

		vector< boost::fusion::tuple< string, string, boost::shared_ptr<DataExchanger> > > domain_domain_list ;
		vector< boost::fusion::tuple< string, string, boost::shared_ptr<DataExchanger> > > surface_domain_list ;
		vector< boost::fusion::tuple< string, string, boost::shared_ptr<DataExchanger> > > domain_surface_list ;

		void clean() ;
		void release() ; // Remember to release MPI communicator before the MPI_Finalize() 
		
		// CELL
		int if_register ( Domain *, Domain * , int ) ;
		int create_cell_exchanger ( Domain *, Domain * ) ;
		void cell_interpolation ( Domain *, Domain *, boost::shared_array<double>, boost::shared_array<double> ) ;

		// FACE
		int if_register ( Domain *, Domain *, string , int ) ;
		int create_face_exchanger ( Domain *, Domain *, string, const vector<int> & ) ;
		void face_interpolation ( Domain *, Domain *, string, boost::shared_array<double>, boost::shared_array<double>, const vector<int> & ) ;
		void face_interpolation ( Domain *, Domain *, boost::shared_array<double>, boost::shared_array<double> ) ;

};


#endif
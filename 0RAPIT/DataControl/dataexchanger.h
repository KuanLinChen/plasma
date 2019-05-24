#include <map>
#include <iostream>
#include <cmath>
#include <sstream>
#include <memory>
#include <boost/shared_array.hpp>
#include "main.h"
#include "domain.h"
//#include <mpi.h>

using namespace std ;

#if !defined(__DATAEXCHANGER_H)
#define __DATAEXCHANGER_H

#define	CELL -100
#define	NODE -101
#define	FACE -102
#define	REFCELL -103
#define	REFNODE -104
#define	REFFACE -105

#define INNERDOMAIN	-1001
#define INTERDOMAIN	-1002


extern void empty_function ( Domain * , Domain * ,boost::shared_array<double>  , boost::shared_array<double>  , void * ) ;
extern void ghost_interpolation_single ( Domain * , Domain * ,boost::shared_array<double>  , boost::shared_array<double>  , void * ) ;
extern void ghost_interpolation ( Domain * , Domain * ,boost::shared_array<double>  , boost::shared_array<double> , void * ) ;
extern void H_interpolation_single_processor ( Domain * , Domain * , boost::shared_array<double>  , boost::shared_array<double>  , void *  ) ;
extern void H_interpolation ( Domain * , Domain * , boost::shared_array<double> , boost::shared_array<double>  , void * ) ;
extern void NH_interpolation_single_processor ( Domain * , Domain * , boost::shared_array<double>  , boost::shared_array<double> , void *  ) ;
extern void NH_interpolation ( Domain * , Domain * , boost::shared_array<double>  , boost::shared_array<double> , void * ) ;
extern void identical_face_copy ( Domain * , Domain * , boost::shared_array<double>  , boost::shared_array<double> , void * ) ;

typedef void (*InterpolationFuncPtr)( Domain * , Domain * , boost::shared_array<double> , boost::shared_array<double> , void *  ); 

class DataExchanger
{
	public:
		DataExchanger () ;
		~DataExchanger () ;

		DataExchanger( Domain * ) ;					// Ghost 
		DataExchanger( Domain *, Domain * ) ; 		//NH


		DataExchanger( Domain *, const vector<int> &, int, string ) ;					// Ghost with spacial index 
		DataExchanger( Domain *, Domain *, const vector<int> &, int, string ) ; 		//NH with spacial index

		MPI_Comm comm ; 

		string DataExchanger_tag ;

		// For ghost/homogeneous updating
		int		mpiSize, mpiId ; 
		int		total_recv_number, total_send_number ;
		boost::shared_array<int> MPISEND_number, MPIRECV_number, MPISEND_index, MPIRECV_index;
		int LOCAL_TRANS_number;
		boost::shared_array<int> LOCAL_TRANS_source_index, LOCAL_TRANS_target_index ;
		boost::shared_array< int >  SENDCargo_index, RECVCargo_index ;
		boost::shared_array< double > SENDCargo, RECVCargo ;

		boost::shared_array< vector<double> > WeightLocal, WeightNonLocal ;
		boost::shared_array< vector<int> > WeightLocalID, WeightNonLocalID ;

		// DataExchangerGhost
		//int	*index_global_to_local, *index_local_to_global ;
		void 	buildCellRelation_Ghost ( Domain *domain  ) ;
		//void	buildREFFaceRelation_Ghost ( Domain *, int, int * ) ; //  domain, data_number, index array

		// DataExchangerNH (NonHomogeneous)
		//int		*source_index_global_to_local, *source_index_local_to_global ;
		//int		*target_index_global_to_local, *target_index_local_to_global ;
		void 	buildCellRelation_H ( Domain * , Domain * ) ;
		void 	buildCellRelation_NH ( Domain * , Domain * ) ;
		

		void	buildFaceRelation_Identical ( Domain *, const vector<int> & , string  ) ;
		void	buildFaceRelation_H ( Domain *, Domain *, const vector<int> & , string ) ;

		
		InterpolationFuncPtr	interpolationFunc ;
		
		//void	remove_relation ( Domain *, Domain * ) ;

		static	int 			map_iterator ;
		static 	map<string, int> exchanger_index ;

	private:
} ;

class DataExchangerL2R: public DataExchanger 
{
	public:
		void buildRelation ( Domain * domain_target, Domain * domain_source ) ;
	//void interpolation ( Domain * domain_target, Domain * domain_source, double * target, double * source ) ;
} ;

extern DataExchanger dataexchanger ;


#endif

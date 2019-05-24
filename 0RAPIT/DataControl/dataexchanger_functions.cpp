#include <mpi.h>
#include "dataexchanger.h"
#include "sys_log.h"


void ghost_interpolation_single ( Domain * domain_target, Domain * domain_source, boost::shared_array<double>  target, boost::shared_array<double>  source ,  void *ptr )  
{
	int i ;

	if ( source != target )
	{
		for ( i = 0 ; i < domain_target->local_cell_number; i++ )
		{
			target [i] = source [i] ;
		}
	}
}


/*! \brief Update the ghost cell data for given array and related domain.

The idea is each processor collects the data that will send to processor i, and, then, use MPI_Gatherv to do so. The processor i will get the organized data in array receivebuffer, which is exactly the same order to the ghost_cell_data in processor i. In this case, the easy way is to attached the pointer receivebuffer to the first element of the array of ghost data of processor i.

*/

void ghost_interpolation ( Domain * domain_target, Domain * domain_source, boost::shared_array<double> target, boost::shared_array<double> source,  void *ptr )
{
	int i, j, *k ;
	DataExchanger *DE  ;
	double *sendbuffer, *recvbuffer ;

	DE = ( DataExchanger * ) ptr ;
	if ( DE->comm == MPI_COMM_NULL ) return;

	recvbuffer = DE->RECVCargo.get() ;
	sendbuffer = DE->SENDCargo.get() ;

	for ( i = 0; i < DE->total_send_number ; i ++ )
	{
		sendbuffer[i] = source[ DE->SENDCargo_index[i] ] ;
	}

	// If we check memory by valgrind, there are uninitialized issue here, which seems the bugs of mpich-3.2 that happens to uninitial variable in MPI. 
	MPI_Alltoallv ( sendbuffer, DE->MPISEND_number.get(), DE->MPISEND_index.get(), MPI_DOUBLE, recvbuffer, DE->MPIRECV_number.get() , DE->MPIRECV_index.get(), MPI_DOUBLE , DE->comm ) ;

	for ( i = 0 ; i < domain_target->ghost_cell_number ; i++ )
	{
		target[ DE->RECVCargo_index[i]  ] = recvbuffer[ i ] ;
	}

	if ( source != target )
	{
		for ( i = 0 ; i < domain_target->local_cell_number; i++ )
			target[i] = source[i] ;
	}

}

/*! \brief Homogeneous mesh interpolation for single processor.

No needs to transfer data from cpu to cpu. Simply copy the data from one to one. 

*/
void H_interpolation_single_processor ( Domain * domain_target, Domain * domain_source, boost::shared_array<double>  target, boost::shared_array<double> source ,  void *ptr )
{
	int i ;

	if ( source != target )
	{
		for ( i = 0 ; i < domain_target->local_cell_number; i++ )
		{
			target [i] = source [i] ;
		}
	}
}

/*! \brief Homogeneous mesh interpolation for multiple processors.



*/
void H_interpolation ( Domain * domain_target, Domain * domain_source, boost::shared_array<double> target, boost::shared_array<double>  source ,  void *ptr )
{
	int i, j, *k ;
	DataExchanger *DE  ;
	double *sendbuffer, *recvbuffer ;


	DE = ( DataExchanger * ) ptr ;

	if ( DE->comm == MPI_COMM_NULL ) return;
	recvbuffer = DE->RECVCargo.get() ;
	sendbuffer = DE->SENDCargo.get() ;

 
	for ( i = 0; i < DE->total_send_number ; i ++ )
	{
		sendbuffer[i] = source[ DE->SENDCargo_index[i] ] ;
	}

	// If we check memory by valgrind, there are uninitialized issue here, which seems the bugs of mpich-3.2 that happens to uninitial variable in MPI. 
	MPI_Alltoallv ( sendbuffer, DE->MPISEND_number.get(), DE->MPISEND_index.get(), MPI_DOUBLE, recvbuffer, DE->MPIRECV_number.get() , DE->MPIRECV_index.get(), MPI_DOUBLE , DE->comm ) ;

	for ( i = 0; i < DE->total_recv_number ; i ++ )
	{
		target[ DE->RECVCargo_index[i] ] = recvbuffer[i] ;
	}

	// Local part
	for ( i = 0; i < DE->LOCAL_TRANS_number ; i ++ )
	{
		target[ DE->LOCAL_TRANS_target_index[i] ] = source[ DE->LOCAL_TRANS_source_index[i] ] ; 
	}

}


void NH_interpolation_single_processor ( Domain * domain_target, Domain * domain_source,boost::shared_array<double>  target, boost::shared_array<double>  source ,  void *ptr )
{
	int i, j;
	DataExchanger	*DE ;

	DE = ( DataExchanger * ) ptr ;

	if ( DE->comm == MPI_COMM_NULL ) return;
	
	// Do data interpolation
	//cout << "data interpolation single processor..." <<  domain_target->meshfile + domain_source->meshfile << endl ;

	for( i = 0 ; i < domain_target->local_cell_number ; i ++  )
	{
		target[i] = 0. ;
		for ( j = 0 ; j < DE->WeightLocalID[i].size() ; j++ )
		{
			target[i] += DE->WeightLocal[i][j] * source[ DE->WeightLocalID[i][j] ] ;
		}
	}

}

/*! \brief Interpolation the data. Send and then interpolate.


*/

void NH_interpolation ( Domain * domain_target, Domain * domain_source, boost::shared_array<double> target, boost::shared_array<double>  source  , void * ptr )
{
	int i, j, k ;
	DataExchanger *DE  ;
	double *sendbuffer, *recvbuffer ;
	int *number_send_from_proc;

	DE = ( DataExchanger * ) ptr ;

	if ( DE->comm == MPI_COMM_NULL ) return;

	recvbuffer = new double [  DE->total_recv_number ] ;

	for ( i = 0; i < DE->total_send_number ; i ++ )
	{
		DE->SENDCargo[i] = source[ DE->SENDCargo_index[i] ] ;
	}
	sendbuffer = DE->SENDCargo.get() ;

	// If we check memory by valgrind, there are uninitialized issue here, which seems the bugs of mpich-3.2 that happens to uninitial variable in MPI. 
	MPI_Alltoallv ( sendbuffer, DE->MPISEND_number.get(), DE->MPISEND_index.get(), MPI_DOUBLE, recvbuffer, DE->MPIRECV_number.get() , DE->MPIRECV_index.get(), MPI_DOUBLE , DE->comm ) ;

	for ( i = 0 ; i < domain_target->local_cell_number ; i ++  )
	{
		target[ i ] = 0. ;
		for ( j = 0 ; j < DE->WeightLocalID[i].size()  ; j++ )
		{
			target[i] +=  DE->WeightLocal[i][j] * source[  DE->WeightLocalID[i][j] ] ;
		}

		for ( j = 0 ; j < DE->WeightNonLocalID[i].size()  ; j++ )
		{
			target[i] += DE->WeightNonLocal[i][j] * recvbuffer [  DE->WeightNonLocalID[i][j] ] ;
		}

	}

	delete [] recvbuffer ;
}

/*! \brief Copy the data from two Scalar face, with the identical domain. 

*/
void identical_face_copy ( Domain * domain_target, Domain * domain_source, boost::shared_array<double> target, boost::shared_array<double>  source ,  void *ptr )
{
	int i ;
	DataExchanger *DE  ;

	DE = ( DataExchanger * ) ptr ;

	if ( DE->comm == MPI_COMM_NULL ) return;
	
	// Local part
	for ( i = 0; i < DE->LOCAL_TRANS_number ; i ++ )
	{
		target[ DE->LOCAL_TRANS_target_index[i] ] = source[ DE->LOCAL_TRANS_source_index[i] ] ; 

		//cout << i << " " << DE->LOCAL_TRANS_target_index[i]  << " "<< DE->LOCAL_TRANS_source_index[i] << endl;
	}

}


void empty_function ( Domain * domain_target, Domain * domain_source, boost::shared_array<double> target, boost::shared_array<double>  source ,  void *ptr )
{

}
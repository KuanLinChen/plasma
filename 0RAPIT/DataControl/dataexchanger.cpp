#include <mpi.h>
#include "dataexchanger.h"
#include "cellbycell_tracking.h"
#include "sys_log.h"
#include <boost/shared_array.hpp>

//DataExchanger		dataexchanger ;
//map<string, int> 	DataExchanger::exchanger_index ;
//vector < unique_ptr< DataExchanger > >	DataExchanger::vec_dataexchanger_ptr ;
//int 				DataExchanger::map_iterator = 0 ;


DataExchanger::DataExchanger()
{
}

DataExchanger::~DataExchanger()
{
}

DataExchanger::DataExchanger( Domain *domain )
{
	//if ( mpi_id == 0 ) cout << "Ghost" << endl ;
	buildCellRelation_Ghost ( domain ) ;
}

/*
DataExchangerNH::DataExchangerNH( Domain * domain_target, Domain * domain_source )
{
	buildCellRelation ( domain_target, domain_source ) ;
}
*/
DataExchanger::DataExchanger ( Domain * domain_target, Domain * domain_source )
{
	if ( domain_target->meshfile == domain_source->meshfile )
	{
		buildCellRelation_H (  domain_target, domain_source ) ;
	}
	else
	{
		buildCellRelation_NH ( domain_target, domain_source ) ;
	}

}

DataExchanger::DataExchanger( Domain *domain, const vector<int> &collectrd_id, int status , string additional_key)
{
	//if ( mpi_id == 0 ) cout << "Ghost" << endl ;
	
	if ( status == CELL )
	{
		buildCellRelation_Ghost ( domain ) ;
	} else if ( status == FACE )
	{
		buildFaceRelation_Identical ( domain, collectrd_id, additional_key ) ;
	}
}

DataExchanger::DataExchanger( Domain *domain_target, Domain *domain_source, const vector<int> &collectrd_id, int status  , string additional_key )
{
	//if ( mpi_id == 0 ) cout << "Ghost" << endl ;
	
	if ( status == CELL )
	{
		if ( domain_target->meshfile == domain_source->meshfile )
		{
			buildCellRelation_H (  domain_target, domain_source ) ;
		} else
		{
			buildCellRelation_NH ( domain_target, domain_source ) ;
		}
	} else if ( status == FACE )
	{
		if ( domain_target->meshfile == domain_source->meshfile )
		{
			buildFaceRelation_H (  domain_target, domain_source , collectrd_id, additional_key ) ;
		}
		else
		{
			//buildCellRelation_NH ( domain_target, domain_source ) ;
		}
			
	}
}

void DataExchanger::buildCellRelation_Ghost ( Domain *domain )
{
	int i, j, k, tag;
	int *buf ;
	MPI_Status status ;

	DataExchanger_tag = "CELL_" + domain->ID + '_' + domain->ID ;

	boost::shared_array< vector<int> > MPIRECV_ID  ;
	boost::shared_array< vector<int> > MPIRECV_LOCAL_ID  ;
	boost::shared_array< vector<int> > MPISEND_ID  ;

	comm = domain->comm ;
	mpiId = -1 ;
	mpiSize = 0 ;

	if (  comm == MPI_COMM_NULL )
	{
		// Out here if the continuation is not my business in this communicator
		interpolationFunc = empty_function ;
		return ;
	}
	else
	{
		MPI_Comm_size ( comm, &mpiSize) ;
		MPI_Comm_rank ( comm, &mpiId ) ;
	}

	MPISEND_number = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_number = boost::shared_array<int> ( new int [mpiSize] );
	MPISEND_index = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_index = boost::shared_array<int> ( new int [mpiSize] );

	for ( i = 0 ; i < mpiSize ; i++ )
	{
		MPISEND_number[i] = 0;
		MPIRECV_number[i] = 0;
		MPISEND_index[i] = 0;
		MPIRECV_index[i] = 0;
	}

	if ( mpiSize > 1)
	{
		/// There are domain->ghost_cell_number belong to this local process ( = mpiID ), and each ghost cell requires the data from other processors ( != mpiID ). number_send_from_proc[i] record how many data cell from processor No. i.

		MPIRECV_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;
		MPIRECV_LOCAL_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;
		MPISEND_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;

		for ( j = domain->local_cell_number ; j < domain->local_cell_number + domain->ghost_cell_number ; j++ )
		{
			k = domain->cell[j].mpi_id ;
			MPIRECV_number[k] ++ ;
			MPIRECV_ID[k].push_back( domain->cell[j].id ) ;
			MPIRECV_LOCAL_ID[k].push_back( domain->cell[j].local_id ) ;
			//if ( domain->comm_rank == 0 && k == 3 ) cout << j << " " << domain->cell[j].id  << endl;
		}

		// To tell other processors there will be number_send_from_proc[i] will be obtained from processor i. Also, to know how many data of cells should send to processor j (number_send_to_proc[j]).

		for ( i = 0 ; i < mpiSize ; i++ )
		{
			tag = 100 ;
			if ( i == mpiId )
			{
				for ( j = i + 1 ; j < mpiSize ; j++ )
				{
					MPI_Recv( &k, 1, MPI_INT, j, tag, comm, &status ) ;
					MPISEND_number[j] = k ;

					k = MPIRECV_number[j] ;
					MPI_Send( &k , 1, MPI_INT, j , tag, comm ) ;
				}
			}
			else if ( mpiId > i )
			{
				k = MPIRECV_number[i] ;
				MPI_Send( &k , 1, MPI_INT, i , tag, comm ) ;

				MPI_Recv( &k, 1, MPI_INT, i, tag, comm, &status ) ;
				MPISEND_number[i] = k ;
			}
		}

		for ( i = 0 ; i < mpiSize ; i++ )
		{
			tag = 100 ;
			if ( i == mpiId )
			{
				for ( j = i + 1 ; j < mpiSize ; j++ )
				{
					buf = new int [ MPISEND_number[j] ] ;
					MPI_Recv( buf, MPISEND_number[j], MPI_INT, j, tag, comm, &status ) ;
					MPISEND_ID[j].reserve( MPISEND_number[j] ) ;
					for ( k = 0 ; k < MPISEND_number[j] ; k++ )
					{
						MPISEND_ID[j].push_back( buf[k] ) ;
					}
					delete [] buf ;

					buf = new int [ MPIRECV_number[j] ] ;
					for ( k = 0 ; k < MPIRECV_number[j] ; k++ )
					{
						buf[k] = MPIRECV_ID[j][k] ;
					}
					MPI_Send( buf , MPIRECV_number[j], MPI_INT, j , tag, comm ) ;
					delete [] buf ;
				}
			}
			else if ( mpiId > i )
			{
				buf = new int [  MPIRECV_number[i] ] ;
				for ( k = 0 ; k < MPIRECV_number[i] ; k++ )
				{
					buf[k] = MPIRECV_ID[i][k] ;
				}
				MPI_Send( buf , MPIRECV_number[i], MPI_INT, i , tag, comm ) ;
				delete [] buf ;

				buf = new int [MPISEND_number[i]] ;
				MPI_Recv( buf, MPISEND_number[i], MPI_INT, i, tag, comm, &status ) ;
				MPISEND_ID[i].reserve( MPISEND_number[i] ) ;
				for ( k = 0 ; k < MPISEND_number[i] ; k++ )
				{
					MPISEND_ID[i].push_back( buf[k] ) ;
				}
				delete [] buf;
			}
		}

		// Since the target and source domain are exactly the same, we take the advantage by using the exist memory without allocating new one. MPIRECV_index will start from the first memory position, exactly the position after the local data.
		// The data that needs to send to other cpu is not continued and we need to collect them before sending. MPISEND_index will start from zero.

		MPIRECV_index [0] = 0;
		MPISEND_index [0] = 0;
		for ( i = 1 ; i < mpiSize ; i++ )
		{
			MPIRECV_index [i] = MPIRECV_number [i - 1] + MPIRECV_index [i - 1] ;
			MPISEND_index [i] = MPISEND_number [i - 1] + MPISEND_index [i - 1] ;
		}

		total_recv_number = domain->ghost_cell_number ;
		total_send_number = 0  ;
		for ( i = 0 ; i < mpiSize ; i++ )
		{
			total_send_number += MPISEND_number[i] ;
		}

		SENDCargo_index = boost::shared_array< int > ( new int [total_send_number ] ) ;
		RECVCargo_index = boost::shared_array< int > ( new int [ total_recv_number ] ) ;


		k = 0 ;
		for ( i = 0 ; i < mpiSize ; i++ )
		{
			for ( j = 0 ; j < MPISEND_ID[i].size() ; j++ )
			{
				//SENDCargo_index[k] = domain->pCell_Cell[ MPISEND_ID[i][j] ] ;
				SENDCargo_index[ k ] = domain->GlobalCell_LocalCellNo[ MPISEND_ID[ i ][ j ] ] ;
				k++ ;
			}
		}

		k = 0 ;
		for ( i = 0 ; i < mpiSize ; i++ )
		{
			for ( j = 0 ; j < MPIRECV_LOCAL_ID[i].size() ; j++ )
			{
				RECVCargo_index[k] = MPIRECV_LOCAL_ID[i][j] ;
				k++ ;
			}
		}

		SENDCargo = boost::shared_array< double > ( new double [total_send_number] ) ;
		RECVCargo = boost::shared_array< double > ( new double [total_recv_number] ) ;

		for ( i = 0 ; i < total_send_number ; i++ )
		{
			SENDCargo [i] = 0. ;
		}

		for ( i = 0 ; i < total_recv_number ; i++ )
		{
			RECVCargo [i] = 0. ;
		}

	}

	if ( mpiSize == 1 )
	{
		interpolationFunc = ghost_interpolation_single ;
	}
	else if (  mpiSize > 1 )
	{
		interpolationFunc = ghost_interpolation ;
	}

}

void DataExchanger::buildCellRelation_H ( Domain *domain_target , Domain *domain_source )
{
	int i, j, k, id, tag, *buf;
	MPI_Status status ;
	MPI_Group	world_group, mpi_group ;

	boost::shared_array <int>  group_rank ;

	int target_size, target_rank, target_root ;
	int source_size, source_rank, source_root ;
	boost::shared_array <int>  source_rank_2_comm_rank, target_rank_2_comm_rank ;
	boost::shared_array <int>  buffer ;
	boost::shared_array <int>  target_CellID2MeshOrdering, target_MeshOrdering2CellID, source_CellID2MeshOrdering, source_MeshOrdering2CellID ;
	boost::shared_array <int>  source_cell_mpi_id, target_cell_mpi_id;

	DataExchanger_tag = "CELL_" + domain_target->ID + '_' + domain_source->ID ;

	// Begin of the communicator construction. The class variable 'comm' will collect cpus from both domain_target->comm and domain_source->comm.

	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		//cout << "Different parent communicator of target and source. Might be dangerous." << endl;
		MPI_Comm_size ( MPI_COMM_WORLD, &mpiSize ) ;
		MPI_Comm_rank ( MPI_COMM_WORLD, &mpiId ) ;
	}
	else
	{
		MPI_Comm_size ( domain_target->parent_comm, &mpiSize ) ;
		MPI_Comm_rank ( domain_target->parent_comm, &mpiId ) ;
	}

	// if the exchanger uses only single CPU, assign the function and leave
	if ( mpiSize == 1 )
	{
		interpolationFunc = H_interpolation_single_processor ;
		return;
	}

	target_size = -1 ;
	target_rank = -1 ;
	source_size = -1 ;
	source_rank = -1 ;

	if (  domain_target->comm != MPI_COMM_NULL ) { MPI_Comm_size ( domain_target->comm, &target_size ) ; }
	if (  domain_target->comm != MPI_COMM_NULL ) { MPI_Comm_rank ( domain_target->comm, &target_rank ) ; }
	if (  domain_source->comm != MPI_COMM_NULL ) { MPI_Comm_size ( domain_source->comm, &source_size ) ; }
	if (  domain_source->comm != MPI_COMM_NULL ) { MPI_Comm_rank ( domain_source->comm, &source_rank ) ; }

	i = -1 ;
	if ( target_rank != -1  )
	{
		i = mpiId ;
	}
	else if (  source_rank != -1 )
	{
		i = mpiId ;
	}
	else
	{
		i = -1 ;
	}

	group_rank = boost::shared_array <int> ( new int [mpiSize] );

	// Root process will collect the array from others. If the process is involving in this exchanger, the mpiId will return. Otherwise, return -1.
	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		MPI_Gather ( &i, 1, MPI_INT, group_rank.get(), 1, MPI_INT, 0, MPI_COMM_WORLD )  ;
		MPI_Bcast ( group_rank.get(), mpiSize, MPI_INT, 0, MPI_COMM_WORLD ) ;
	}
	else
	{
		MPI_Gather ( &i, 1, MPI_INT, group_rank.get(), 1, MPI_INT, 0, domain_target->parent_comm )  ;
		MPI_Bcast ( group_rank.get(), mpiSize, MPI_INT, 0, domain_target->parent_comm  ) ;
	}

	// Bobble sorting to make the involving CPU IDs in order
	i = mpiSize - 1 ;
	while ( i >= 0 )
	{
		for ( j = 0 ; j < i - 1 ; j++ )
		{
			if ( group_rank[j + 1] == -1 )
			{
				// do nothing
			}
			else if ( group_rank[j] == -1 )
			{
				group_rank[j] = group_rank[j + 1] ;
				group_rank[j + 1] = -1 ;
			}
			else if ( group_rank[j] > group_rank[j + 1] )
			{
				k = group_rank[j + 1] ;
				group_rank[j] = group_rank[j + 1] ;
				group_rank[j + 1] = k;
			}
		}
		i-- ;
	}

	j = mpiSize;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		if ( group_rank[i] == -1 )
		{
			j = i;
			i = mpiSize;   // to stop the loop;
		}
	}
	mpiSize = j ;
	
	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		MPI_Comm_group( MPI_COMM_WORLD, &world_group ) ;
		MPI_Group_incl( world_group, mpiSize, group_rank.get(), &mpi_group );
		MPI_Comm_create_group( MPI_COMM_WORLD, mpi_group, 0, &comm );
		MPI_Group_free ( &world_group ) ;
		MPI_Group_free ( &mpi_group ) ;
	}
	else
	{
		MPI_Comm_group( domain_target->parent_comm , &world_group ) ;
		MPI_Group_incl( world_group, mpiSize, group_rank.get(), &mpi_group );
		MPI_Comm_create_group( domain_target->parent_comm, mpi_group, 0, &comm );
		MPI_Group_free ( &world_group ) ;
		MPI_Group_free ( &mpi_group ) ;
	}
	//comm = domain_target->parent_comm;

	// Leave if this process is not involved in the new communicator
	if (  comm == MPI_COMM_NULL )
	{
		// Out here if the continuation is not my business in this communicator
		interpolationFunc = empty_function ;
		return ;
	}
	else
	{
		MPI_Comm_size ( comm, &mpiSize) ;
		MPI_Comm_rank ( comm, &mpiId ) ;
	}

	buffer = boost::shared_array<int> ( new int [mpiSize] ) ;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		buffer[i] =  -1 ;
	}
	MPI_Gather ( &source_size, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	if ( mpiId == 0  )
	{
		for ( i = 0 ; i < mpiSize ; i++)
		{
			if ( buffer[i] != -1 )
			{
				source_size = buffer[i] ;
			}
		}
	}
	MPI_Bcast ( &source_size, 1, MPI_INT, 0, comm );

	source_rank_2_comm_rank = boost::shared_array<int> ( new int [source_size] ) ;
	for ( i = 0 ; i < source_size ; i++ )
	{
		source_rank_2_comm_rank[ i ] = -1 ;
	}

	MPI_Gather ( &source_rank, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	MPI_Bcast ( buffer.get(), mpiSize, MPI_INT, 0, comm );
	for ( i = 0 ; i < mpiSize; i++ )
	{
		if ( buffer[i] != -1 )
		{
			source_rank_2_comm_rank[ buffer[i] ] = i ;
		}
	}

	MPI_Gather ( &target_size, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	if ( mpiId == 0 )
	{
		for ( i = 0 ; i < mpiSize ; i++)
			if ( buffer[i] != -1 )
			{
				target_size = buffer[i] ;
			}
	}
	MPI_Bcast ( &target_size, 1, MPI_INT, 0, comm );

	target_rank_2_comm_rank = boost::shared_array<int> ( new int [target_size] ) ;
	for ( i = 0 ; i < target_size ; i++ )
	{
		target_rank_2_comm_rank[i] = -1;
	}

	MPI_Gather ( &target_rank, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	MPI_Bcast ( buffer.get(), mpiSize, MPI_INT, 0, comm );
	for ( i = 0 ; i < mpiSize; i++ )
	{
		if ( buffer[i] != -1 )
		{
			target_rank_2_comm_rank[ buffer[i] ] = i ;
		}
	}

	MPI_Bcast ( &(domain_target->global_cell_number), 1, MPI_INT, target_rank_2_comm_rank[0], comm );
	MPI_Bcast ( &(domain_source->global_cell_number), 1, MPI_INT, source_rank_2_comm_rank[0], comm );

	target_CellID2MeshOrdering = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
	target_MeshOrdering2CellID = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
	source_CellID2MeshOrdering = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
	source_MeshOrdering2CellID = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
	if ( target_rank == 0 )
	{
		for ( i = 0 ; i < domain_target->global_cell_number ; i++ )
		{
			//target_CellID2MeshOrdering[i] = domain_target->CellID2MeshOrdering[i] ;
			//target_MeshOrdering2CellID[i] = domain_target->MeshOrdering2CellID[i] ;
			target_CellID2MeshOrdering[ i ]	=	domain_target->GlobalCell_MeshCellNo[ i ] ;
			target_MeshOrdering2CellID[ i ]	=	domain_target->MeshCell_GlobalCellNo[ i ] ;
		}
	}
	MPI_Bcast ( target_CellID2MeshOrdering.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm );
	MPI_Bcast ( target_MeshOrdering2CellID.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm );

	if ( source_rank == 0 )
	{
		for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
		{
			//source_CellID2MeshOrdering[i] = domain_source->CellID2MeshOrdering[i] ;
			//source_MeshOrdering2CellID[i] = domain_source->MeshOrdering2CellID[i] ;
			source_CellID2MeshOrdering[ i ] =	domain_source->GlobalCell_MeshCellNo[ i ] ;
			source_MeshOrdering2CellID[ i ] =	domain_source->MeshCell_GlobalCellNo[ i ] ;
		}
	}

	MPI_Bcast ( source_CellID2MeshOrdering.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );
	MPI_Bcast ( source_MeshOrdering2CellID.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );

	source_cell_mpi_id = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
	if ( source_rank == 0 )
	{
		for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
		{
			//source_cell_mpi_id[i] = domain_source->global_cell[i].mpi_id ;
			source_cell_mpi_id[ i ]	=	domain_source->Mesh.Cell_ProcessorNo[ domain_source->GlobalCell_MeshCellNo[ i ] ] ;
		}
	}
	MPI_Bcast ( source_cell_mpi_id.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );

	target_cell_mpi_id = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
	if ( target_rank == 0 )
	{
		for ( i = 0 ; i < domain_target->global_cell_number ; i++ )
		{
			//target_cell_mpi_id[i] = domain_target->global_cell[i].mpi_id  ;
			target_cell_mpi_id[ i ] = domain_target->Mesh.Cell_ProcessorNo[ domain_target->GlobalCell_MeshCellNo[ i ] ]  ;
		}
	}
	MPI_Bcast ( target_cell_mpi_id.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm ) ;
	// End of communicator construction. Above code give the mpiId, mpiSize and comm. The cpu who is not involving the communicator will not participate the following code.

	boost::shared_array< vector<int> > MPIRECV_ID  ;
	boost::shared_array< vector<int> > MPISEND_ID  ;
	vector<int> SOURCE_LOCALTRANS_ID, TARGET_LOCALTRANS_ID   ;

	MPISEND_number = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_number = boost::shared_array<int> ( new int [mpiSize] );
	MPISEND_index = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_index = boost::shared_array<int> ( new int [mpiSize] );
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		MPISEND_number[i] = 0;
		MPIRECV_number[i] = 0;
		MPISEND_index[i] = 0;
		MPIRECV_index[i] = 0;
	}

	/// There are domain->ghost_cell_number belong to this local process ( = mpiID ), and each ghost cell requires the data from other processors ( != mpiID ). number_send_from_proc[i] record how many data cell from processor No. i.

	MPIRECV_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;
	MPISEND_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;

	// Since target and source are using the same meshfile, we use the CellID2MeshOrdering and MeshOrdering2CellID to find the relation between these 2 domain.
	if ( domain_target->comm != MPI_COMM_NULL )
	{
		for ( j = 0 ; j < domain_target->local_cell_number ; j++ )
		{
			k = target_CellID2MeshOrdering [ domain_target->cell[j].id ] ;
			i = source_MeshOrdering2CellID[k] ;

			id = source_rank_2_comm_rank[ source_cell_mpi_id[i] ] ;

			if ( id != target_rank_2_comm_rank [ domain_target->cell[j].mpi_id] )   // From other CPU
			{
				MPIRECV_number[id] ++ ;
				MPIRECV_ID[id].push_back( target_CellID2MeshOrdering [ domain_target->cell[j].id ] ) ;
			} else     // From local CPU
			{
				//SOURCE_LOCALTRANS_ID.push_back( domain_source->pCell_Cell[ i ] ) ;
				SOURCE_LOCALTRANS_ID.push_back( domain_source->GlobalCell_LocalCellNo[ i ] ) ;
				TARGET_LOCALTRANS_ID.push_back( j ) ;
			}
		}
	}

	LOCAL_TRANS_number = SOURCE_LOCALTRANS_ID.size() ;
	LOCAL_TRANS_source_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
	LOCAL_TRANS_target_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
	for ( i = 0 ; i < SOURCE_LOCALTRANS_ID.size() ; i++ )
	{
		LOCAL_TRANS_source_index[i] = SOURCE_LOCALTRANS_ID[i] ;
		LOCAL_TRANS_target_index[i] = TARGET_LOCALTRANS_ID[i] ;
	}

	// To tell other processors there will be number_send_from_proc[i] will be obtained from processor i. Also, to know how many data of cells should send to processor j (number_send_to_proc[j]).
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		tag = 100 ;
		if ( i == mpiId )
		{
			for ( j = i + 1 ; j < mpiSize ; j++ )
			{
				MPI_Recv( &k, 1, MPI_INT, j, tag, comm, &status ) ;
				MPISEND_number[j] = k ;

				k = MPIRECV_number[j] ;
				MPI_Send( &k , 1, MPI_INT, j , tag, comm ) ;
			}
		}
		else if ( mpiId > i )
		{
			k = MPIRECV_number[i] ;
			MPI_Send( &k , 1, MPI_INT, i , tag, comm ) ;

			MPI_Recv( &k, 1, MPI_INT, i, tag, comm, &status ) ;
			MPISEND_number[i] = k ;
		}
	}

	for ( i = 0 ; i < mpiSize ; i++ )
	{
		tag = 101 ;
		if ( i == mpiId )
		{
			for ( j = i + 1 ; j < mpiSize ; j++ )
			{
				buf = new int [ MPISEND_number[j] ] ;
				MPI_Recv( buf, MPISEND_number[j], MPI_INT, j, tag, comm, &status ) ;
				MPISEND_ID[j].reserve( MPISEND_number[j] ) ;
				for ( k = 0 ; k < MPISEND_number[j] ; k++ )
				{
					MPISEND_ID[j].push_back( buf[k]) ;
				}
				delete [] buf ;

				buf = new int [ MPIRECV_number[j] ] ;
				for ( k = 0 ; k < MPIRECV_number[j] ; k++ )
				{
					buf[k] = MPIRECV_ID[j][k] ;
				}
				MPI_Send( buf , MPIRECV_number[j], MPI_INT, j , tag, comm ) ;
				delete [] buf ;
			}
		}
		else if ( mpiId > i )
		{
			buf = new int [  MPIRECV_number[i] ] ;
			for ( k = 0 ; k < MPIRECV_number[i] ; k++ )
			{
				buf[k] = MPIRECV_ID[i][k] ;
			}
			MPI_Send( buf , MPIRECV_number[i], MPI_INT, i , tag, comm ) ;
			delete [] buf ;

			buf = new int [MPISEND_number[i]] ;
			MPI_Recv( buf, MPISEND_number[i], MPI_INT, i, tag, comm, &status ) ;
			MPISEND_ID[i].reserve( MPISEND_number[i] ) ;
			for ( k = 0 ; k < MPISEND_number[i] ; k++ )
			{
				MPISEND_ID[i].push_back(  buf[k] ) ;
			}
			delete [] buf;
		}
	}

	// Since the target and source domain are exactly the same, we take the advantage by using the exist memory without allocating new one. MPIRECV_index will start from the first memory position, exactly the position after the local data.
	// The data that needs to send to other cpu is not continued and we need to collect them before sending. MPISEND_index will start from zero.

	MPIRECV_index [0] = 0;
	MPISEND_index [0] = 0;
	for ( i = 1 ; i < mpiSize ; i++ )
	{
		MPIRECV_index [i] = MPIRECV_number [i - 1] + MPIRECV_index [i - 1] ;
		MPISEND_index [i] = MPISEND_number [i - 1] + MPISEND_index [i - 1] ;
	}

	total_recv_number = 0 ;
	total_send_number = 0 ;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		total_send_number += MPISEND_number[i] ;
		total_recv_number += MPIRECV_number[i] ;
	}

	SENDCargo_index = boost::shared_array< int > ( new int [total_send_number ] ) ;
	RECVCargo_index = boost::shared_array< int > ( new int [ total_recv_number ] ) ;

	SENDCargo = boost::shared_array< double > ( new double [total_send_number] ) ;
	RECVCargo = boost::shared_array< double > ( new double [total_recv_number] ) ;

	k = 0 ;
	for ( i = 0 ; i < mpiSize ; i++ )
		for ( j = 0 ; j < MPISEND_ID[i].size() ; j++ )
		{
			//SENDCargo_index[k] = domain_source->pCell_Cell[ source_MeshOrdering2CellID[ MPISEND_ID[i][j] ] ] ;
			SENDCargo_index[k] = domain_source->GlobalCell_LocalCellNo[ source_MeshOrdering2CellID[ MPISEND_ID[ i ][ j ] ] ] ;
			k++ ;
		}

	k = 0 ;
	for ( i = 0 ; i < mpiSize ; i++ )
		for ( j = 0 ; j < MPIRECV_ID[i].size() ; j++ )
		{
			//RECVCargo_index[k] = domain_target->pCell_Cell[ target_MeshOrdering2CellID[ MPIRECV_ID[i][j] ] ];
			RECVCargo_index[ k ] = domain_target->GlobalCell_LocalCellNo[ target_MeshOrdering2CellID[ MPIRECV_ID[ i ][ j ] ] ] ;
			k++ ;
		}

	interpolationFunc = H_interpolation ;


}


void DataExchanger::buildCellRelation_NH ( Domain * domain_target, Domain * domain_source )
{
	int i, j, k, l, m, tag, *buf;
	int target_size, target_rank, source_size, source_rank ;
	MPI_Group world_group, mpi_group;
	MPI_Status 	status ;

	boost::shared_array <int> buffer;
	boost::shared_array <int> group_rank;
	boost::shared_array <int> source_rank_2_comm_rank, target_rank_2_comm_rank;
	boost::shared_array <int> target_CellID2MeshOrdering, target_MeshOrdering2CellID, source_CellID2MeshOrdering, source_MeshOrdering2CellID;
	boost::shared_array <int> source_cell_mpi_id, target_cell_mpi_id;

	CellByCellRayTracking CCRT( domain_source ) ;
	boost::shared_array <double> x, y, z;
	boost::shared_array <int> CCRT_id ;

	int closest_node_id;
	double dummy, lx2, ly2, lz2, distance, total_weight;

	boost::shared_array< set<int> >	set_from_proc ;
	boost::shared_array< vector<double> > interpolation_weight ;
	boost::shared_array< vector<int> > interpolation_cell_id ;
	boost::shared_array< vector<int> > interpolation_proc_id ;
	//set<int>	set_local_index, set_nonlocal_index ;

	DataExchanger_tag = "CELL_" + domain_target->ID + '_' + domain_source->ID ;

	// Begin of the communicator construction. The class variable 'comm' will collect cpus from both domain_target->comm and domain_source->comm.
	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		Log().TagDump( logLEVEL4 ) << "Different parent communicator: [t] " << domain_target->parent_comm << " [s] " <<  domain_source->parent_comm ;
		MPI_Comm_size ( MPI_COMM_WORLD, &mpiSize ) ;
		MPI_Comm_rank ( MPI_COMM_WORLD, &mpiId ) ;
	}
	else
	{
		MPI_Comm_size ( domain_target->parent_comm, &mpiSize ) ;
		MPI_Comm_rank ( domain_target->parent_comm, &mpiId ) ;
	}

	// if the exchanger uses only single CPU, assign the function and leave
	if ( mpiSize == 1 )
	{
		comm = MPI_COMM_WORLD ;
		source_rank_2_comm_rank = boost::shared_array<int> ( new int [1] ) ;
		source_rank_2_comm_rank[0] = 0 ;

		source_cell_mpi_id = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
		for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
		{
			source_cell_mpi_id[ i ] = 0 ;
		}
	}
	else
	{
		target_size = -1 ;
		target_rank = -1 ;
		source_size = -1 ;
		source_rank = -1 ;

		if (  domain_target->comm != MPI_COMM_NULL ) { MPI_Comm_size ( domain_target->comm, &target_size ) ; }
		if (  domain_target->comm != MPI_COMM_NULL ) { MPI_Comm_rank ( domain_target->comm, &target_rank ) ; }
		if (  domain_source->comm != MPI_COMM_NULL ) { MPI_Comm_size ( domain_source->comm, &source_size ) ; }
		if (  domain_source->comm != MPI_COMM_NULL ) { MPI_Comm_rank ( domain_source->comm, &source_rank ) ; }

		i = -1 ;
		if ( target_rank != -1  )
		{
			i = mpiId ;
		}
		else if (  source_rank != -1 )
		{
			i = mpiId ;
		}
		else
		{
			i = -1 ;
		}

		group_rank = boost::shared_array <int> ( new int [mpiSize] );

		// Root process will collect the array from others. If the process is involving in this exchanger, the mpiId will return. Otherwise, return -1.
		if ( domain_target->parent_comm != domain_source->parent_comm )
		{
			MPI_Gather ( &i, 1, MPI_INT, group_rank.get(), 1, MPI_INT, 0, MPI_COMM_WORLD )  ;
			MPI_Bcast ( group_rank.get(), mpiSize, MPI_INT, 0, MPI_COMM_WORLD ) ;
		}
		else
		{
			MPI_Gather ( &i, 1, MPI_INT, group_rank.get(), 1, MPI_INT, 0, domain_target->parent_comm )  ;
			MPI_Bcast ( group_rank.get(), mpiSize, MPI_INT, 0, domain_target->parent_comm  ) ;
		}

		// Bobble sorting to make the involving CPU IDs in order
		i = mpiSize - 1 ;
		while ( i >= 0 )
		{
			for ( j = 0 ; j < i - 1 ; j++ )
			{
				if ( group_rank[j + 1] == -1 )
				{
					// do nothing
				}
				else if ( group_rank[j] == -1 )
				{
					group_rank[j] = group_rank[j + 1] ;
					group_rank[j + 1] = -1 ;
				}
				else if ( group_rank[j] > group_rank[j + 1] )
				{
					k = group_rank[j + 1] ;
					group_rank[j] = group_rank[j + 1] ;
					group_rank[j + 1] = k;
				}
			}
			i-- ;
		}

		j = mpiSize;
		for ( i = 0 ; i < mpiSize ; i++ )
		{
			if ( group_rank[i] == -1 )
			{
				j = i;
				i = mpiSize;   // to stop the loop;
			}
		}
		mpiSize = j ;

		/*
		if ( domain_target->parent_comm != domain_source->parent_comm )
		{
			MPI_Comm_group( MPI_COMM_WORLD, &world_group ) ;
			MPI_Group_incl( world_group, mpiSize, group_rank.get(), &mpi_group );
			MPI_Comm_create_group( MPI_COMM_WORLD, mpi_group, 0, &comm );
			MPI_Group_free ( &world_group ) ;
			MPI_Group_free ( &mpi_group ) ;
		}
		else
		{
			MPI_Comm_group( domain_target->parent_comm , &world_group ) ;
			MPI_Group_incl( world_group, mpiSize, group_rank.get(), &mpi_group );
			MPI_Comm_create_group( domain_target->parent_comm, mpi_group, 0, &comm );
			MPI_Group_free ( &world_group ) ;
			MPI_Group_free ( &mpi_group ) ;
		}
		*/
		comm = domain_target->parent_comm;

		// Leave if this process is not involved in the new communicator
		if (  comm == MPI_COMM_NULL )
		{
			// Out here if the continuation is not my business in this communicator
			interpolationFunc = empty_function ;
			return ;
		}
		else
		{
			MPI_Comm_size ( comm, &mpiSize) ;
			MPI_Comm_rank ( comm, &mpiId ) ;
		}

		buffer = boost::shared_array<int> ( new int [mpiSize] ) ;
		for ( i = 0 ; i < mpiSize ; i++ )
		{
			buffer[i] =  -1 ;
		}
		MPI_Gather ( &source_size, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
		if ( mpiId == 0  )
		{
			for ( i = 0 ; i < mpiSize ; i++)
			{
				if ( buffer[i] != -1 )
				{
					source_size = buffer[i] ;
				}
			}
		}
		MPI_Bcast ( &source_size, 1, MPI_INT, 0, comm );

		source_rank_2_comm_rank = boost::shared_array<int> ( new int [source_size] ) ;
		for ( i = 0 ; i < source_size ; i++ )
		{
			source_rank_2_comm_rank[ i ] = -1 ;
		}

		MPI_Gather ( &source_rank, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
		MPI_Bcast ( buffer.get(), mpiSize, MPI_INT, 0, comm );
		for ( i = 0 ; i < mpiSize; i++ )
		{
			if ( buffer[i] != -1 )
			{
				source_rank_2_comm_rank[ buffer[i] ] = i ;
			}
		}

		MPI_Gather ( &target_size, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
		if ( mpiId == 0 )
		{
			for ( i = 0 ; i < mpiSize ; i++)
				if ( buffer[i] != -1 )
				{
					target_size = buffer[i] ;
				}
		}
		MPI_Bcast ( &target_size, 1, MPI_INT, 0, comm );

		target_rank_2_comm_rank = boost::shared_array<int> ( new int [target_size] ) ;
		for ( i = 0 ; i < target_size ; i++ )
		{
			target_rank_2_comm_rank[i] = -1;
		}

		MPI_Gather ( &target_rank, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
		MPI_Bcast ( buffer.get(), mpiSize, MPI_INT, 0, comm );
		for ( i = 0 ; i < mpiSize; i++ )
		{
			if ( buffer[i] != -1 )
			{
				target_rank_2_comm_rank[ buffer[i] ] = i ;
			}
		}

		MPI_Bcast ( &(domain_target->global_cell_number), 1, MPI_INT, target_rank_2_comm_rank[0], comm );
		MPI_Bcast ( &(domain_source->global_cell_number), 1, MPI_INT, source_rank_2_comm_rank[0], comm );

		target_CellID2MeshOrdering = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
		target_MeshOrdering2CellID = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
		source_CellID2MeshOrdering = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
		source_MeshOrdering2CellID = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;

		if ( target_rank == 0 )
		{
			for ( i = 0 ; i < domain_target->global_cell_number ; i++ )
			{
				//target_CellID2MeshOrdering[i] = domain_target->CellID2MeshOrdering[i] ;
				//target_MeshOrdering2CellID[i] = domain_target->MeshOrdering2CellID[i] ;
				target_CellID2MeshOrdering[ i ]	=	domain_target->GlobalCell_MeshCellNo[ i ] ;
				target_MeshOrdering2CellID[ i ]	=	domain_target->MeshCell_GlobalCellNo[ i ] ;
			}
		}
		MPI_Bcast ( target_CellID2MeshOrdering.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm );
		MPI_Bcast ( target_MeshOrdering2CellID.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm );

		if ( source_rank == 0 )
		{
			for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
			{
				//source_CellID2MeshOrdering[i] = domain_source->CellID2MeshOrdering[i] ;
				//source_MeshOrdering2CellID[i] = domain_source->MeshOrdering2CellID[i] ;
				source_CellID2MeshOrdering[ i ] =	domain_source->GlobalCell_MeshCellNo[ i ] ;
				source_MeshOrdering2CellID[ i ] =	domain_source->MeshCell_GlobalCellNo[ i ] ;	
			}
		}
		MPI_Bcast ( source_CellID2MeshOrdering.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );
		MPI_Bcast ( source_MeshOrdering2CellID.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );

		source_cell_mpi_id = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
		if ( source_rank == 0 )
		{
			for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
			{
				//source_cell_mpi_id[i] = domain_source->global_cell[i].mpi_id  ;
				source_cell_mpi_id[ i ] =	domain_source->Mesh.Cell_ProcessorNo[ domain_source->GlobalCell_MeshCellNo[ i ] ]  ;
			}
		}
		MPI_Bcast ( source_cell_mpi_id.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );

		target_cell_mpi_id = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
		if ( target_rank == 0 )
		{
			for ( i = 0 ; i < domain_target->global_cell_number ; i++ )
			{
				//target_cell_mpi_id[i] = domain_target->global_cell[i].mpi_id  ;
				target_cell_mpi_id[ i ]	=	domain_target->Mesh.Cell_ProcessorNo[ domain_target->GlobalCell_MeshCellNo[ i ] ]  ;
			}
		}
		MPI_Bcast ( target_cell_mpi_id.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm ) ;
	}
	// End of communicator construction. Above code give the mpiId, mpiSize and comm. The cpu who is not involving the communicator will not participate the following code.

	x = boost::shared_array <double> ( new double [ domain_target->local_cell_number ] ) ;
	y = boost::shared_array <double> ( new double [ domain_target->local_cell_number ] ) ;
	z = boost::shared_array <double> ( new double [ domain_target->local_cell_number ] ) ;
	CCRT_id = boost::shared_array <int> ( new int [ domain_target->local_cell_number ] ) ;

	for ( i = 0 ; i < domain_target->local_cell_number ; i++ )
	{
		x[i] = domain_target->cell[i].x ;
		y[i] = domain_target->cell[i].y ;
		z[i] = domain_target->cell[i].z ;
		CCRT_id[i] = -1 ;
	}

	//CCRT << 10 ;
	CCRT.PointsMapping( domain_target->local_cell_number, x.get(), y.get(), z.get(), CCRT_id.get() ) ;

	interpolation_weight = boost::shared_array< vector<double> > ( new vector<double> [domain_target->local_cell_number] ) ;
	interpolation_cell_id = boost::shared_array< vector<int> > ( new vector<int> [domain_target->local_cell_number] ) ;
	interpolation_proc_id = boost::shared_array< vector<int> > ( new vector<int> [domain_target->local_cell_number] ) ;

	for ( i = 0 ; i < domain_target->local_cell_number ; i++ )
	{
		if ( CCRT_id[i] < 0  )
		{
			if ( domain_target->cell[ i ].type == -1 )
			{
				Log().TagDump( logLEVEL4 ) << "Cell ID: " << i << " is ignored since the cell type is SOLID." ;
			}
			else
			{
				Log().TagDump( logERROR ) << "Something wrong at CCRT. ID = " << i << ". " ;
			}
		} else
		{

			dummy	=	1.e20 ;
			closest_node_id	=	-1 ;

			for ( j = 0 ; j < domain_source->Mesh.Cell_Node[ domain_source->GlobalCell_MeshCellNo[ CCRT_id[ i ] ] ].size() ; j++ )
			{
				k 	=	domain_source->Mesh.Cell_Node[ domain_source->GlobalCell_MeshCellNo[ CCRT_id[ i ] ] ][ j ] ;
				lx2 =	( x[ i ] - domain_source->Mesh.Node_Position[ 0 ][ k ] ) * ( x[ i ] - domain_source->Mesh.Node_Position[ 0 ][ k ] ) ;
				ly2 =	( y[ i ] - domain_source->Mesh.Node_Position[ 1 ][ k ] ) * ( y[ i ] - domain_source->Mesh.Node_Position[ 1 ][ k ] ) ;
				lz2 =	( z[ i ] - domain_source->Mesh.Node_Position[ 2 ][ k ] ) * ( z[ i ] - domain_source->Mesh.Node_Position[ 2 ][ k ] ) ;

				distance =	sqrt ( lx2 + ly2 + lz2 ) ;
				if ( dummy > distance )
				{
					dummy 	=	distance ;
					closest_node_id =	k ;
				}
			}			
		}

		//interpolation_weight[i].reserve ( domain_source->pre_Node[ closest_node_id ].CellRelation.size() ) ;
		//interpolation_cell_id[i].reserve ( domain_source->pre_Node[ closest_node_id ].CellRelation.size() ) ;
		//interpolation_proc_id[i].reserve ( domain_source->pre_Node[ closest_node_id ].CellRelation.size() ) ;		
		interpolation_weight[ i ].reserve ( domain_source->Mesh.Node_Cell[ closest_node_id ].size() ) ;
		interpolation_cell_id[ i ].reserve ( domain_source->Mesh.Node_Cell[ closest_node_id ].size() ) ;
		interpolation_proc_id[ i ].reserve ( domain_source->Mesh.Node_Cell[ closest_node_id ].size() ) ;

		for ( j = 0 ; j < domain_source->Mesh.Node_Cell[ closest_node_id ].size() ; j++ )
		{
			//interpolation_cell_id[ i ].push_back( domain_source->pre_Cell[ domain_source->pre_Node[ closest_node_id ].CellRelation[j] ].id ) ;
			//interpolation_proc_id[ i ].push_back( source_rank_2_comm_rank [ source_cell_mpi_id [ interpolation_cell_id[i][j] ] ] ) ;
			interpolation_cell_id[ i ].push_back( domain_source->MeshCell_GlobalCellNo[ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] ) ;
			interpolation_proc_id[ i ].push_back( source_rank_2_comm_rank [ source_cell_mpi_id [ interpolation_cell_id[ i ][ j ] ] ] ) ;			
		}

		
		total_weight =	0. ;
		for ( j = 0 ; j < domain_source->Mesh.Node_Cell[ closest_node_id ].size(); j++ )
		{
			lx2 	=	( domain_source->Mesh.Cell_Position[ 0 ][ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] - domain_target->Mesh.Cell_Position[ 0 ][ domain_target->LocalCell_MeshCellNo[ i ] ] ) * ( domain_source->Mesh.Cell_Position[ 0 ][ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] - domain_target->Mesh.Cell_Position[ 0 ][ domain_target->LocalCell_MeshCellNo[ i ] ] ) ;
			ly2 	=	( domain_source->Mesh.Cell_Position[ 1 ][ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] - domain_target->Mesh.Cell_Position[ 1 ][ domain_target->LocalCell_MeshCellNo[ i ] ] ) * ( domain_source->Mesh.Cell_Position[ 1 ][ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] - domain_target->Mesh.Cell_Position[ 1 ][ domain_target->LocalCell_MeshCellNo[ i ] ] ) ;
			lz2 	=	( domain_source->Mesh.Cell_Position[ 2 ][ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] - domain_target->Mesh.Cell_Position[ 2 ][ domain_target->LocalCell_MeshCellNo[ i ] ] ) * ( domain_source->Mesh.Cell_Position[ 2 ][ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] - domain_target->Mesh.Cell_Position[ 2 ][ domain_target->LocalCell_MeshCellNo[ i ] ] ) ;
			dummy 	=	sqrt ( lx2 + ly2 + lz2 ) ;

			if ( dummy < 1.e-20 )
			{
				// This is the case while target cell is exactly at the source cell.
				interpolation_cell_id[ i ].clear() ;
				interpolation_proc_id[ i ].clear() ;

				interpolation_cell_id[ i ].push_back( domain_source->MeshCell_GlobalCellNo[ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] ) ;
				interpolation_proc_id[ i ].push_back( domain_source->Mesh.Cell_ProcessorNo[ domain_source->Mesh.Node_Cell[ closest_node_id ][ j ] ] ) ;
				j =  domain_source->Mesh.Node_Cell[ closest_node_id ].size() ; // Stop the looping
			} else
			{
				interpolation_weight[ i ].push_back ( 1. / dummy ) ;
				total_weight += 1.0 / dummy ;
			}
		}		

		if ( interpolation_cell_id[i].size() == 1 )
		{
			interpolation_weight[i].clear() ;
			interpolation_weight[i].push_back( 1.0 ) ;
		}
		else
		{
			for ( j = 0 ; j <  domain_source->Mesh.Node_Cell[ closest_node_id ].size(); j++ )
			{
				interpolation_weight[ i ][ j ]	= 	interpolation_weight[ i ][ j ] / total_weight ;
			}
		}
	}

	if ( mpiSize == 1 )
	{
		WeightLocalID	=	boost::shared_array< vector<int> > ( new vector<int> [ domain_target->local_cell_number ] ) ;
		WeightLocal		=	boost::shared_array< vector<double> > ( new vector<double> [ domain_target->local_cell_number ] ) ;
		for ( i = 0 ; i < domain_target->local_cell_number ; i++ )
		{
			for ( j = 0 ; j < interpolation_weight[ i ].size() ; j++ )
			{
				if ( domain_source->Mesh.Cell_ProcessorNo[ domain_source->GlobalCell_MeshCellNo[ interpolation_cell_id[ i ][ j ] ] ] == mpiId )
				{
					WeightLocalID[ i ].push_back ( domain_source->GlobalCell_LocalCellNo[ interpolation_cell_id[ i ][ j ] ] ) ;
					WeightLocal[ i ].push_back( interpolation_weight[ i ][ j ] ) ;
				}
			}

		}

		interpolationFunc = NH_interpolation_single_processor ;
		return ;
	}

	boost::shared_array< vector<int> > MPIRECV_ID  ;
	boost::shared_array< vector<int> > MPISEND_ID  ;
	vector<int> SOURCE_LOCALTRANS_ID, TARGET_LOCALTRANS_ID   ;

	MPISEND_number = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_number = boost::shared_array<int> ( new int [mpiSize] );
	MPISEND_index = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_index = boost::shared_array<int> ( new int [mpiSize] );
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		MPISEND_number[i] = 0;
		MPIRECV_number[i] = 0;
	}

	/// There are domain->ghost_cell_number belong to this local process ( = mpiID ), and each ghost cell requires the data from other processors ( != mpiID ). number_send_from_proc[i] record how many data cell from processor No. i.
	MPIRECV_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;
	MPISEND_ID = boost::shared_array< vector<int> > ( new vector<int> [mpiSize] ) ;

	set_from_proc = boost::shared_array< set<int> > ( new set<int> [ mpiSize ] ) ;
	for ( i = 0 ; i < domain_target->local_cell_number ; i++ )
	{
		for ( j = 0 ; j < interpolation_weight[i].size() ; j++ )
		{
			if ( ( mpiId != interpolation_proc_id [i][j]) )
			{
				set_from_proc[ interpolation_proc_id [i][j] ].insert( interpolation_cell_id [i][j] ) ;
			}
		}
	}

	for ( i = 0 ; i < mpiSize ; i++ )
	{
		MPIRECV_number[i] = set_from_proc[i].size() ;
		for ( set<int>::iterator it = set_from_proc[i].begin() ; it != set_from_proc[i].end()  ; it++  )
		{
			MPIRECV_ID[i].push_back ( *it ) ;
		}
	}

	// To tell other processors there will be number_send_from_proc[i] will be obtained from processor i. Also, to know how many data of cells should send to processor j (number_send_to_proc[j]).
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		tag = 100 ;
		if ( i == mpiId )
		{
			for ( j = i + 1 ; j < mpiSize ; j++ )
			{
				MPI_Recv( &k, 1, MPI_INT, j, tag, comm, &status ) ;
				MPISEND_number[j] = k ;

				k = MPIRECV_number[j] ;
				MPI_Send( &k , 1, MPI_INT, j , tag, comm ) ;
			}
		}
		else if ( mpiId > i )
		{
			k = MPIRECV_number[i] ;
			MPI_Send( &k , 1, MPI_INT, i , tag, comm ) ;
			MPI_Recv( &k, 1, MPI_INT, i, tag, comm, &status ) ;
			MPISEND_number[i] = k ;
		}
	}

	for ( i = 0 ; i < mpiSize ; i++ )
	{
		tag = 101 ;
		if ( i == mpiId )
		{
			for ( j = i + 1 ; j < mpiSize ; j++ )
			{
				buf = new int [ MPISEND_number[j] ] ;
				MPI_Recv( buf, MPISEND_number[j], MPI_INT, j, tag, comm, &status ) ;
				MPISEND_ID[j].reserve( MPISEND_number[j] ) ;
				for ( k = 0 ; k < MPISEND_number[j] ; k++ )
				{
					MPISEND_ID[j].push_back( buf[k] ) ;
				}
				delete [] buf ;

				buf = new int [ MPIRECV_number[j] ] ;
				for ( k = 0 ; k < MPIRECV_number[j] ; k++ )
				{
					buf[k] = MPIRECV_ID[j][k] ;
				}
				MPI_Send( buf , MPIRECV_number[j], MPI_INT, j , tag, comm ) ;
				delete [] buf ;
			}
		}
		else if ( mpiId > i )
		{
			buf = new int [  MPIRECV_number[i] ] ;
			for ( k = 0 ; k < MPIRECV_number[i] ; k++ )
			{
				buf[k] = MPIRECV_ID[i][k] ;
			}
			MPI_Send( buf , MPIRECV_number[i], MPI_INT, i , tag, comm ) ;
			delete [] buf ;

			buf = new int [MPISEND_number[i]] ;
			MPI_Recv( buf, MPISEND_number[i], MPI_INT, i, tag, comm, &status ) ;
			MPISEND_ID[i].reserve( MPISEND_number[i] ) ;
			for ( k = 0 ; k < MPISEND_number[i] ; k++ )
			{
				MPISEND_ID[i].push_back(  buf[k] ) ;
			}
			delete [] buf;
		}
	}

	MPIRECV_index [0] = 0;
	MPISEND_index [0] = 0;
	for ( i = 1 ; i < mpiSize ; i++ )
	{
		MPIRECV_index [i] = MPIRECV_number [i - 1] + MPIRECV_index [i - 1] ;
		MPISEND_index [i] = MPISEND_number [i - 1] + MPISEND_index [i - 1] ;
	}

	total_recv_number = 0 ;
	total_send_number = 0 ;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		total_send_number += MPISEND_number[i] ;
		total_recv_number += MPIRECV_number[i] ;
	}

	SENDCargo_index = boost::shared_array< int > ( new int [total_send_number ] ) ;
	RECVCargo_index = boost::shared_array< int > ( new int [ total_recv_number ] ) ;

	k = 0 ;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		for ( j = 0 ; j < MPISEND_ID[i].size() ; j++ )
		{
			//SENDCargo_index[k] = domain_source->pCell_Cell[  MPISEND_ID[i][j] ];
			SENDCargo_index[ k ]	=	domain_source->GlobalCell_LocalCellNo[  MPISEND_ID[ i ][ j ] ] ;
			k++ ;
		}
	}

	k = 0 ;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		for ( j = 0 ; j < MPIRECV_ID[i].size() ; j++ )
		{
			RECVCargo_index[k] = MPIRECV_ID[i][j] ;
			k++ ;
		}
	}

	SENDCargo = boost::shared_array< double > ( new double [total_send_number] ) ;
	RECVCargo = boost::shared_array< double > ( new double [total_recv_number] ) ;


	WeightLocalID	= boost::shared_array< vector<int> > ( new vector<int> [domain_target->local_cell_number] ) ;
	WeightNonLocalID	= boost::shared_array< vector<int> > ( new vector<int> [domain_target->local_cell_number] ) ;
	WeightLocal	= boost::shared_array< vector<double> > ( new vector<double> [domain_target->local_cell_number] ) ;
	WeightNonLocal	= boost::shared_array< vector<double> > ( new vector<double> [domain_target->local_cell_number] ) ;

	for ( i = 0 ; i < domain_target->local_cell_number ; i++ )
	{
		for ( j = 0 ; j < interpolation_weight[i].size() ; j++ )
		{
			//if ( domain_source->pre_Cell[ interpolation_cell_id[ i ][ j ] ].ProcessorNo == mpiId )
			if ( domain_source->Mesh.Cell_ProcessorNo[ domain_source->GlobalCell_LocalCellNo[ interpolation_cell_id[ i ][ j ] ] ] == mpiId )
			{
				// source from the same cpu
				//WeightLocalID[ i ].push_back ( domain_source->pCell_Cell[ interpolation_cell_id[i][j] ] ) ;
				WeightLocalID[ i ].push_back ( domain_source->GlobalCell_LocalCellNo[ interpolation_cell_id[ i ][ j ] ] ) ;
				WeightLocal[ i ].push_back( interpolation_weight[ i ][ j ] ) ;
				//Log( ).TagDump( logLEVEL4 ) << i << " " << domain_source->pCell_Cell[ interpolation_cell_id[i][j] ] << " " <<interpolation_weight[i][j] ;
			} else
			{
				m = 0 ;
				// source from the different cpu
				for ( k = 0 ; k < mpiSize ; k++ )
				{
					tag =  MPIRECV_ID[k].size() ;
					for ( l = 0 ; l < tag ; l++ )
					{
						if ( MPIRECV_ID[ k ][ l ] == interpolation_cell_id[ i ][ j ] )
						{
							WeightNonLocalID[ i ].push_back( m ) ;
							l = tag ;
							k = mpiSize ;
						}
						m++ ;
					}
				}
				WeightNonLocal[i].push_back( interpolation_weight[i][j] ) ;
				//Log( ).TagDump( logLEVEL4 ) << i << " " << interpolation_cell_id[i][j] << " " << m << " " << interpolation_weight[i][j] ;
			}
		}

	}

	// NonHomogeneous interpolation. Send and then Interpolate.
	interpolationFunc = NH_interpolation ;

}

/*! \brief Face data copy.

This is for the case that both target and source are belong to identical domain.

*/
void DataExchanger::buildFaceRelation_Identical ( Domain *domain, const vector<int>  &collectrd_id, string additional_key )
{
	int i, j, k, tag;
	int *buf ;
	MPI_Status status ;

	DataExchanger_tag = "FACE_" + domain->ID + '_' + domain->ID + '_' + additional_key ;


	boost::shared_array< vector<int> > MPIRECV_ID  ;
	boost::shared_array< vector<int> > MPIRECV_LOCAL_ID  ;
	boost::shared_array< vector<int> > MPISEND_ID  ;

	comm = domain->comm ;
	mpiId = -1 ;
	mpiSize = 0 ;

	if (  comm == MPI_COMM_NULL )
	{
		// Out here if the continuation is not my business in this communicator
		interpolationFunc = empty_function ;
		return ;
	} else
	{
		MPI_Comm_size ( comm, &mpiSize) ;
		MPI_Comm_rank ( comm, &mpiId ) ;
	}

	MPISEND_number = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_number = boost::shared_array<int> ( new int [mpiSize] );
	MPISEND_index = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_index = boost::shared_array<int> ( new int [mpiSize] );

	for ( i = 0 ; i < mpiSize ; i++ )
	{
		MPISEND_number[i] 	= 0;
		MPIRECV_number[i] 	= 0;
		MPISEND_index[i]	= 0;
		MPIRECV_index[i] 	= 0;
	}

	// Ghost face is playing passive role in the program. In the case of Ghost updating, we actually don't need ghost face exchanging. Local face copy is enough. 
	LOCAL_TRANS_number = collectrd_id.size() ;
	LOCAL_TRANS_source_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
	LOCAL_TRANS_target_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
	for ( i = 0 ; i < LOCAL_TRANS_number ; i++ )
	{
		LOCAL_TRANS_source_index[i] = collectrd_id[i] ;
		LOCAL_TRANS_target_index[i] = collectrd_id[i] ;
	}
	
	// Ghost interpolation. Send and then Interpolate.
	interpolationFunc = identical_face_copy ;
}


/*! \brief Face Update the ghost cell data for given array and related domain.



*/
void DataExchanger::buildFaceRelation_H ( Domain *domain_target,Domain *domain_source, const vector<int>  &collectrd_id, string additional_key )
{
	int i, j, k, tag;
	int *buf ;
	MPI_Status status ;
	MPI_Group	world_group, mpi_group ;

	boost::shared_array <int>  group_rank ;

	int target_size, target_rank, target_root ;
	int source_size, source_rank, source_root ;
	boost::shared_array <int>  source_rank_2_comm_rank, target_rank_2_comm_rank ;
	boost::shared_array <int>  buffer ;
	boost::shared_array <int>  target_CellID2MeshOrdering, target_MeshOrdering2CellID, source_CellID2MeshOrdering, source_MeshOrdering2CellID ;
	boost::shared_array <int>  source_cell_mpi_id, target_cell_mpi_id;

	DataExchanger_tag = "FACE_" + domain_target->ID + '_' + domain_source->ID + '_' + additional_key ;

	boost::shared_array< vector<int> > MPIRECV_ID  ;
	boost::shared_array< vector<int> > MPIRECV_LOCAL_ID  ;
	boost::shared_array< vector<int> > MPISEND_ID  ;

	// Begin of the communicator construction. The class variable 'comm' will collect cpus from both domain_target->comm and domain_source->comm.

	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		//cout << "Different parent communicator of target and source. Might be dangerous." << endl;
		MPI_Comm_size ( MPI_COMM_WORLD, &mpiSize ) ;
		MPI_Comm_rank ( MPI_COMM_WORLD, &mpiId ) ;
	} else
	{
		MPI_Comm_size ( domain_target->parent_comm, &mpiSize ) ;
		MPI_Comm_rank ( domain_target->parent_comm, &mpiId ) ;
	}

	// if the exchanger uses only single CPU, assign the function and leave
	if ( mpiSize == 1 )
	{
		LOCAL_TRANS_number = collectrd_id.size() ;
		LOCAL_TRANS_source_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
		LOCAL_TRANS_target_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
		for ( i = 0 ; i < LOCAL_TRANS_number ; i++ )
		{
			LOCAL_TRANS_source_index[i] = collectrd_id[i] ;
			LOCAL_TRANS_target_index[i] = domain_target->MeshFace_LocalFaceNo [ domain_source->LocalFace_MeshFaceNo [ collectrd_id[i] ] ] ;
		}
		interpolationFunc = identical_face_copy ;
		return;
	}

	target_size = -1 ;
	target_rank = -1 ;
	source_size = -1 ;
	source_rank = -1 ;

	if (  domain_target->comm != MPI_COMM_NULL ) { MPI_Comm_size ( domain_target->comm, &target_size ) ; }
	if (  domain_target->comm != MPI_COMM_NULL ) { MPI_Comm_rank ( domain_target->comm, &target_rank ) ; }
	if (  domain_source->comm != MPI_COMM_NULL ) { MPI_Comm_size ( domain_source->comm, &source_size ) ; }
	if (  domain_source->comm != MPI_COMM_NULL ) { MPI_Comm_rank ( domain_source->comm, &source_rank ) ; }

	i = -1 ;
	if ( target_rank != -1  )
	{
		i = mpiId ;
	} else if (  source_rank != -1 )
	{
		i = mpiId ;
	} 

	group_rank = boost::shared_array <int> ( new int [mpiSize] );

	// Root process will collect the array from others. If the process is involving in this exchanger, the mpiId will return. Otherwise, return -1.
	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		MPI_Gather ( &i, 1, MPI_INT, group_rank.get(), 1, MPI_INT, 0, MPI_COMM_WORLD )  ;
		MPI_Bcast ( group_rank.get(), mpiSize, MPI_INT, 0, MPI_COMM_WORLD ) ;
	}
	else
	{
		MPI_Gather ( &i, 1, MPI_INT, group_rank.get(), 1, MPI_INT, 0, domain_target->parent_comm )  ;
		MPI_Bcast ( group_rank.get(), mpiSize, MPI_INT, 0, domain_target->parent_comm  ) ;
	}

	// Bobble sorting to make the involving CPU IDs in order
	i = mpiSize - 1 ;
	while ( i >= 0 )
	{
		for ( j = 0 ; j < i - 1 ; j++ )
		{
			if ( group_rank[j + 1] == -1 )
			{
				// do nothing
			}
			else if ( group_rank[j] == -1 )
			{
				group_rank[j] = group_rank[j + 1] ;
				group_rank[j + 1] = -1 ;
			}
			else if ( group_rank[j] > group_rank[j + 1] )
			{
				k = group_rank[j + 1] ;
				group_rank[j] = group_rank[j + 1] ;
				group_rank[j + 1] = k;
			}
		}
		i-- ;
	}

	// Ignore the cpu who is not involving in both source and target communicator
	j = mpiSize;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		if ( group_rank[i] == -1 )
		{
			j = i;
			i = mpiSize;   // to stop the loop;
		}
	}
	mpiSize = j ;
	
	if ( domain_target->parent_comm != domain_source->parent_comm )
	{
		MPI_Comm_group( MPI_COMM_WORLD, &world_group ) ;
		MPI_Group_incl( world_group, mpiSize, group_rank.get(), &mpi_group );
		MPI_Comm_create_group( MPI_COMM_WORLD, mpi_group, 0, &comm );
		MPI_Group_free ( &world_group ) ;
		MPI_Group_free ( &mpi_group ) ;
	}
	else
	{
		MPI_Comm_group( domain_target->parent_comm , &world_group ) ;
		MPI_Group_incl( world_group, mpiSize, group_rank.get(), &mpi_group );
		MPI_Comm_create_group( domain_target->parent_comm, mpi_group, 0, &comm );
		MPI_Group_free ( &world_group ) ;
		MPI_Group_free ( &mpi_group ) ;
	}

	// Leave if this process is not involved in the new communicator
	if (  comm == MPI_COMM_NULL )
	{
		// Out here if the continuation is not my business in this communicator
		interpolationFunc = empty_function ;
		return ;
	}
	else
	{
		MPI_Comm_size ( comm, &mpiSize) ;
		MPI_Comm_rank ( comm, &mpiId ) ;
	}

	// Get source_size from anyone of the source group
	buffer = boost::shared_array<int> ( new int [mpiSize] ) ;
	for ( i = 0 ; i < mpiSize ; i++ )
	{
		buffer[i] =  -1 ;
	}
	MPI_Gather ( &source_size, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	if ( mpiId == 0  )
	{
		for ( i = 0 ; i < mpiSize ; i++)
		{
			if ( buffer[i] != -1 )
			{
				source_size = buffer[i] ;
			}
		}
	}
	MPI_Bcast ( &source_size, 1, MPI_INT, 0, comm );

	// Build the relation from source_rank to recent comm_rank. -1 means the comm_rank is not involving the source data
	source_rank_2_comm_rank = boost::shared_array<int> ( new int [source_size] ) ;
	for ( i = 0 ; i < source_size ; i++ )
	{
		source_rank_2_comm_rank[ i ] = -1 ;
	}
	MPI_Gather ( &source_rank, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	MPI_Bcast ( buffer.get(), mpiSize, MPI_INT, 0, comm );
	for ( i = 0 ; i < mpiSize; i++ )
	{
		if ( buffer[i] != -1 )
		{
			source_rank_2_comm_rank[ buffer[i] ] = i ;
		}
	}

	// Build the relation from target_rank to recent comm_rank. -1 means the comm_rank is not involving the target data
	MPI_Gather ( &target_size, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	if ( mpiId == 0 )
	{
		for ( i = 0 ; i < mpiSize ; i++)
			if ( buffer[i] != -1 )
			{
				target_size = buffer[i] ;
			}
	}
	MPI_Bcast ( &target_size, 1, MPI_INT, 0, comm );

	target_rank_2_comm_rank = boost::shared_array<int> ( new int [target_size] ) ;
	for ( i = 0 ; i < target_size ; i++ )
	{
		target_rank_2_comm_rank[i] = -1;
	}

	MPI_Gather ( &target_rank, 1, MPI_INT, buffer.get(), 1, MPI_INT, 0, comm ) ;
	MPI_Bcast ( buffer.get(), mpiSize, MPI_INT, 0, comm );
	for ( i = 0 ; i < mpiSize; i++ )
	{
		if ( buffer[i] != -1 )
		{
			target_rank_2_comm_rank[ buffer[i] ] = i ;
		}
	}

	cout << source_size << " " << source_rank << " " <<  source_root << endl; 







/*


	MPI_Bcast ( &(domain_target->global_cell_number), 1, MPI_INT, target_rank_2_comm_rank[0], comm );
	MPI_Bcast ( &(domain_source->global_cell_number), 1, MPI_INT, source_rank_2_comm_rank[0], comm );

	target_CellID2MeshOrdering = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
	target_MeshOrdering2CellID = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
	source_CellID2MeshOrdering = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
	source_MeshOrdering2CellID = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
	if ( target_rank == 0 )
	{
		for ( i = 0 ; i < domain_target->global_cell_number ; i++ )
		{
			//target_CellID2MeshOrdering[i] = domain_target->CellID2MeshOrdering[i] ;
			//target_MeshOrdering2CellID[i] = domain_target->MeshOrdering2CellID[i] ;
			target_CellID2MeshOrdering[ i ]	=	domain_target->GlobalCell_MeshCellNo[ i ] ;
			target_MeshOrdering2CellID[ i ]	=	domain_target->MeshCell_GlobalCellNo[ i ] ;
		}
	}
	MPI_Bcast ( target_CellID2MeshOrdering.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm );
	MPI_Bcast ( target_MeshOrdering2CellID.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm );

	if ( source_rank == 0 )
	{
		for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
		{
			//source_CellID2MeshOrdering[i] = domain_source->CellID2MeshOrdering[i] ;
			//source_MeshOrdering2CellID[i] = domain_source->MeshOrdering2CellID[i] ;
			source_CellID2MeshOrdering[ i ] =	domain_source->GlobalCell_MeshCellNo[ i ] ;
			source_MeshOrdering2CellID[ i ] =	domain_source->MeshCell_GlobalCellNo[ i ] ;
		}
	}

	MPI_Bcast ( source_CellID2MeshOrdering.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );
	MPI_Bcast ( source_MeshOrdering2CellID.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );

	source_cell_mpi_id = boost::shared_array< int > ( new int [ domain_source->global_cell_number ] ) ;
	if ( source_rank == 0 )
	{
		for ( i = 0 ; i < domain_source->global_cell_number ; i++ )
		{
			//source_cell_mpi_id[i] = domain_source->global_cell[i].mpi_id ;
			source_cell_mpi_id[ i ]	=	domain_source->Mesh.Cell_ProcessorNo[ domain_source->GlobalCell_MeshCellNo[ i ] ] ;
		}
	}
	MPI_Bcast ( source_cell_mpi_id.get(), domain_source->global_cell_number, MPI_INT, source_rank_2_comm_rank[0], comm );

	target_cell_mpi_id = boost::shared_array< int > ( new int [ domain_target->global_cell_number ] ) ;
	if ( target_rank == 0 )
	{
		for ( i = 0 ; i < domain_target->global_cell_number ; i++ )
		{
			//target_cell_mpi_id[i] = domain_target->global_cell[i].mpi_id  ;
			target_cell_mpi_id[ i ] = domain_target->Mesh.Cell_ProcessorNo[ domain_target->GlobalCell_MeshCellNo[ i ] ]  ;
		}
	}
	MPI_Bcast ( target_cell_mpi_id.get(), domain_target->global_cell_number, MPI_INT, target_rank_2_comm_rank[0], comm ) ;
	// End of communicator construction. Above code give the mpiId, mpiSize and comm. The cpu who is not involving the communicator will not participate the following code.

	MPISEND_number = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_number = boost::shared_array<int> ( new int [mpiSize] );
	MPISEND_index = boost::shared_array<int> ( new int [mpiSize] );
	MPIRECV_index = boost::shared_array<int> ( new int [mpiSize] );

	for ( i = 0 ; i < mpiSize ; i++ )
	{
		MPISEND_number[i] 	= 0;
		MPIRECV_number[i] 	= 0;
		MPISEND_index[i]	= 0;
		MPIRECV_index[i] 	= 0;
	}

	// Ghost face is playing passive role in the program. In the case of Ghost updating, we actually don't need ghost face exchanging. Local face copy is enough. 
	LOCAL_TRANS_number = collectrd_id.size() ;
	LOCAL_TRANS_source_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
	LOCAL_TRANS_target_index = boost::shared_array<int> ( new int [LOCAL_TRANS_number] ) ;
	for ( i = 0 ; i < LOCAL_TRANS_number ; i++ )
	{
		LOCAL_TRANS_source_index[i] = collectrd_id[i] ;
		LOCAL_TRANS_target_index[i] = collectrd_id[i] ;
	}

	*/

}




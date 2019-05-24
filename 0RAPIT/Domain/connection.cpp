#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>
#include <cmath>
#include <mpi.h>
#include "domain.h"
#include "cellbycell_tracking.h"
#include "sys_log.h"
#include "connection.h"


using namespace std ;


/*!

@file connection.cpp Functions for building the relation from given points to a Domain. 

*/


/*! \brief  Building connection from the points 

\param source_number a given number tells how many points we have.
\param position_x	The array of x-coord with source_number elements.
\param position_y	The array of y-coord with source_number elements.
\param position_z	The array of z-coord with source_number elements.
\param connection_number The array that returns how many related cell to each points. It is negative is we cannot find the related cell. 
\param related_local_id The array that returns the related cell local id to each points. Note that the local id migh not belong to the mpi_id it self.
\param related_global_id The array that returns the related cell global id to each points. It is negative is we cannot find the related cell. 
\param related_mpi_id The array that returns the related cell global id to each points. It is negative is we cannot find the related cell. 
\param weighing The array that returns the related cell weighing to each points. Rightnow we build this accroading the inversed distance. 
*/

void build_connection ( int source_number, boost::shared_array<double> position_x, boost::shared_array<double> position_y, boost::shared_array<double> position_z, boost::shared_array<int> connection_number,  boost::shared_array< vector<int>  > related_local_id, boost::shared_array< vector<int> > related_global_id, boost::shared_array< vector<int> > related_mpi_id, boost::shared_array< vector<double> > weighing, Domain *BaseDomain ) 
{
	int 	i, j, k;
	int 	cell_id, closest_node_id, node_id ;
	int 	tag ;
	int 	*CCRT_id ;
	int		mpi_size, mpi_rank ;
	double	*x, *y, *z, *w ;
	double	shortest_dist, dist, lx2, ly2, lz2, px, py, pz ;
	int 	*buffer_ptr ;
	MPI_Status status;
	MPI_Comm comm ;

	boost::shared_array < set<int> > global_id_from ;
	boost::shared_array < set<int> > global_id_to ;
	boost::shared_array < set<int> > local_id_from ;
	boost::shared_array < set<int> > local_id_to ;
	boost::shared_array < int > number_from ;
	boost::shared_array < int > number_to ;

	CellByCellRayTracking CCRT( BaseDomain ) ;

	mpi_size 	=	BaseDomain->comm_size ;
	mpi_rank 	=	BaseDomain->comm_rank ;
	comm 		=	BaseDomain->comm ;

	CCRT_id = new int [ source_number ] ; /**< Closest cell global ID */
	for ( i = 0 ; i < source_number ; i++ )
	{
		CCRT_id[ i ] = -1 ;
	}
	CCRT.PointsMapping( source_number, position_x.get(), position_y.get(), position_z.get(), CCRT_id ) ;

	for ( i = 0 ; i < source_number ; i++ )
	{
		if ( CCRT_id[ i ] < 0 )
		{
			Log().MPITagDump( logWARNING ) << "Cannot find related Cell in Domain(" << BaseDomain->meshfile << "). ID = " << i ;
			
			connection_number[ i ] = -1 ;
		} else 
		{
			// Since we know the closest cell is CCRT_id[ i ]. We need to find the closest Node from this cell.
			// cell_id = CCRT_id[ i ] ;
			// Global_id
			cell_id			=	CCRT_id[ i ] ;
			// Mesh_id
			cell_id 		=	BaseDomain->GlobalCell_MeshCellNo[ cell_id ] ;
			shortest_dist 	=	1.e20 ;
			closest_node_id =	-1 ;
			//for ( j = 0 ; j < BaseDomain->pre_Cell[ cell_id ].NodeNum ; j++ )
			for ( j = 0 ; j < BaseDomain->Mesh.Cell_Node[ cell_id ].size() ; j++ ) 
			{
				//node_id = BaseDomain->pre_Cell[ cell_id ].pNode[ j ]->id ;
				//px = BaseDomain->pre_Node[node_id].x ;
				//py = BaseDomain->pre_Node[node_id].y ;
				//pz = BaseDomain->pre_Node[node_id].z ;
				node_id =	BaseDomain->Mesh.Cell_Node[ cell_id ][ j ] ;
				px 		=	BaseDomain->Mesh.Node_Position[ 0 ][ node_id ] ;
				py 		=	BaseDomain->Mesh.Node_Position[ 1 ][ node_id ] ;
				pz 		=	BaseDomain->Mesh.Node_Position[ 2 ][ node_id ] ;

				lx2 	=	( px - position_x[ i ] ) * ( px - position_x[ i ] ) ;
				ly2 	=	( py - position_y[ i ] ) * ( py - position_y[ i ] ) ;
				lz2 	=	( pz - position_z[ i ] ) * ( pz - position_z[ i ] ) ;

				dist = sqrt ( lx2 + ly2 + lz2 ) ;

				if ( shortest_dist > dist ) 
				{
					shortest_dist 	=	dist ;
					closest_node_id =	node_id ;
				}
			}

			if ( closest_node_id < 0 )
			{
				Log().MPITagDump( logERROR ) << "Trouble to find the relation of point ID " << i ;
			} else 
			{
				//connection_number[ i ] = BaseDomain->pre_Node[closest_node_id].CellRelation.size() ;
				connection_number[ i ]	=	BaseDomain->Mesh.Node_Cell[ closest_node_id ].size() ;

				for ( j = 0 ; j < connection_number[ i ] ; j++ ) 
				{
					//cell_id = BaseDomain->pre_Node[closest_node_id].CellRelation[j] ;
					// Mesh_id
					cell_id =	BaseDomain->Mesh.Node_Cell[ closest_node_id ][ j ] ;
					// Global_id
					cell_id	=	BaseDomain->MeshCell_GlobalCellNo[ cell_id ] ;
					related_global_id[ i ].push_back( cell_id )	;
				}
			}
		}	
	}
	delete [] CCRT_id ;

/*
The above section give us number of related Cell and their global cell ID, for example:

Point ID / Number of related Cell / global cell ID / 
	0	/	3	/	1, 2, 3
	1	/	4	/	2, 7, 9, 10
*/


/*
	The following part it to know the CPU id of each cell (related_mpi_id)

Point ID / Number of related Cell / global cell ID / related_mpi_id
	0	/	3	/	1, 2, 3			/ 0, 0, 0 
	1	/	4	/	2, 7, 9, 10		/ 0, 0, 0, 1

*/
 	for ( i = 0 ; i < source_number ; i++ )
 	{
 		for ( j = 0 ; j < connection_number[ i ] ; j++ ) 
		{
			// Global_id
			cell_id =	related_global_id[ i ][ j ] ;
			// Mesh_id
			cell_id =	BaseDomain->GlobalCell_MeshCellNo[ cell_id ] ;
			//related_mpi_id[i].push_back(  BaseDomain->pre_Cell[ cell_id ].ProcessorNo ) ;
			related_mpi_id[ i ].push_back( BaseDomain->Mesh.Cell_ProcessorNo[ cell_id ] ) ;
		}
 	}

/*
	Keep the list of cell not from local cpu. 

Point ID / Number of related Cell / global cell ID / related_mpi_id
	0	/	3	/	1, 2, 3			/ 0, 0, 0 
	1	/	4	/	2, 7, 9, 10		/ 0, 0, 0, 1
	2	/	4	/	3, 8, 9, 11		/ 0, 0, 0, 2
	3	/	4	/	10, 11, 12, 13	/ 1, 2, 2, 3

If mpi_id is 0, we will collect the list of cell data from CPU1, CPU2, CPU3 ...
The global_id_from[1] ( number_from[1] = 1 )
	10

The global_id_from[2] ( number_from[2] = 2 )
	11
	12

The global_id_from[3] ( number_from[3] = 1 )
	13

*/
	global_id_from = boost::shared_array < set<int> > ( new set<int> [ mpi_size ] ) ;
 	for ( i = 0 ; i < source_number ; i++ )
 	{
 		for ( j = 0 ; j < connection_number[ i ] ; j++ ) 
		{
			if ( related_mpi_id[ i ][ j ] != mpi_rank )
			{
				global_id_from[ related_mpi_id[ i ][ j ] ].insert( related_global_id[ i ][ j ] ) ;
			}
		}
 	}
 	number_from = boost::shared_array < int > ( new int [ mpi_size ] ) ;
 	for ( i = 0 ; i < mpi_size ; i++ )
 	{
 		number_from[ i ]	=	global_id_from[ i ].size() ;
 	}

/*
To send the request to other CPU to tell them the information needed for interpolation, including the number of cells and cell ids.


Point ID / Number of related Cell / global cell ID / related_mpi_id
	0	/	3	/	1, 2, 3			/ 0, 0, 0 
	1	/	4	/	2, 7, 9, 10		/ 0, 0, 0, 1
	2	/	4	/	3, 8, 9, 11		/ 0, 0, 0, 2
	3	/	4	/	10, 11, 12, 13	/ 1, 2, 2, 3

====== CPU[0] =======
The global_id_from[1] ( number_from[1] = 1 )
	10

The global_id_from[2] ( number_from[2] = 2 )
	11
	12

The global_id_from[3] ( number_from[3] = 1 )
	13

====== CPU[1] =======
number_to[0] = 1
global_id_to[0]
	10 

====== CPU[2] =======
number_to[0] = 2 
global_id_to[0]
	11
	12

====== CPU[3] =======
number_to[0] = 1 
global_id_to[0]
	13

*/
	for ( i = 0 ; i < mpi_size ; i++ )
	{
		tag = 100 ;
		if ( i != mpi_rank )
		{
			j =	number_from[ i ] ;
			MPI_Send( &j, 1, MPI_INT, i, tag, comm ) ;
		} else
		{
			for ( j = 0 ; j < mpi_size ; j++ )
			{
				if ( j != mpi_rank )
				{
					MPI_Recv( &k, 1, MPI_INT, j, tag, comm, &status ) ;
					number_to[ j ]	=	k ;
				}
			}
		}
	}

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		tag = 200 ;
		if ( i != mpi_rank )
		{
			j =	number_from[ i ] ;
			buffer_ptr = new int [ j ] ;
			k = 0 ;
			for ( set<int>::iterator it = global_id_from[ i ].begin() ; it != global_id_from[ i ].end() ; it++ )
			{
				buffer_ptr[ k ]	=	*it ;
				k++ ;
			} 
			MPI_Send( buffer_ptr, j, MPI_INT, i, tag, comm ) ;
			delete [] buffer_ptr ;
		} else
		{
			for ( j = 0 ; j < mpi_size ; j++ )
			{
				if ( j != mpi_rank )
				{
					k =	number_to[ j ] ;
					buffer_ptr = new int [ k ] ;
					MPI_Recv( buffer_ptr, k, MPI_INT, j, tag, comm, &status ) ;
					for ( k = 0 ; k < number_to[ j ] ; k++ )
						global_id_to[ j ].insert( buffer_ptr[ k ] ) ; 
					delete [] buffer_ptr ;
				}
			}
		}
	}


/*
Each cpu will have the list. The list provide the exchange information before the interpolation. The local_id_to[other cpu] is the data that need to send to other cpu. And the sending data from others will obtain with the order in global_id_from [other cpu]. Since the global_id_from from other cpu is not belong to local CPU, we have no local id to sign with. 

====== CPU[0] =======
global_id_from[1] ( number_from[1] = 1 )
	10				

global_id_from[2] ( number_from[2] = 2 )
	11				
	12				


global_id_to[1] local_id_to[1] ( number_to[1] = 1 )
	a					a

global_id_to[2] local_id_to[1] ( number_to[2] = 2 )
	b					b
	c					c
*/

	for ( i = 0 ; i < BaseDomain->local_cell_number ; i++ )
	{
		for ( j = 0 ; j < mpi_size ; j++ )
		{
			cell_id =	BaseDomain->cell[ i ].id ;
			if ( global_id_to[ j ].count ( cell_id ) == 1 )
			{
				local_id_to[ j ].insert ( BaseDomain->cell[ i ].local_id ) ;
			}
		}
	}


/*

The weighing coefficient calculation:

Point ID / Number of related Cell / global cell ID / related_mpi_id / weighing 
	0	/	3	/	1, 2, 3			/ 0, 0, 0 		/ 0.3, 0.4, 0.3
	1	/	4	/	2, 7, 9, 10		/ 0, 0, 0, 1	/ 0.2, 0.2, 0.2, 0.4
	2	/	4	/	3, 8, 9, 11		/ 0, 0, 0, 2	/ 0.2, 0.2, 0.2, 0.4
	3	/	4	/	10, 11, 12, 13	/ 1, 2, 2, 3 	/ 0.2, 0.2, 0.2, 0.4
*/

	for ( i = 0 ; i < source_number ; i++ )
	{
		x =	new double [ connection_number[ i ] + 1 ] ;
		y = new double [ connection_number[ i ] + 1 ] ;
		z = new double [ connection_number[ i ] + 1 ] ;
		w = new double [ connection_number[ i ] ] ;

		x[ 0 ]	=	position_x[ i ] ;
		y[ 0 ]	=	position_y[ i ] ;
		z[ 0 ]	=	position_z[ i ] ;

		for ( j = 0 ; j < connection_number[ i ] ; j++ )
		{
			// Global_id
			cell_id 	=	related_global_id[ i ][ j ] ;
			// Mesh_id
			cell_id		=	BaseDomain->GlobalCell_MeshCellNo[ cell_id ] ;

			//x[ j + 1 ] =  BaseDomain->pre_Cell[ cell_id ].x ;
			//y[ j + 1 ] =  BaseDomain->pre_Cell[ cell_id ].y ;
			//z[ j + 1 ] =  BaseDomain->pre_Cell[ cell_id ].z ;
			x[ j + 1 ]	=	BaseDomain->Mesh.Cell_Position[ 0 ][ cell_id ] ;
			y[ j + 1 ]	=	BaseDomain->Mesh.Cell_Position[ 1 ][ cell_id ] ;
			z[ j + 1 ]	=	BaseDomain->Mesh.Cell_Position[ 2 ][ cell_id ] ;
		}

		weighing_inversed_distance ( connection_number[ i ], x, y, z, w ) ;

		delete [] x ;
		delete [] y ;
		delete [] z ;
		delete [] w ;
	}
}


/*
	Calculating weight by inversed distance. 

	x[0], y[0], z[0] is the original point. 
	xyz[1 to (base_numer+1)] are the weighting bases. 

	if xyz[0] is too cloase to one of the base, the weighing will assign to that base point entirely.

*/
void weighing_inversed_distance ( int base_number, double *x, double *y, double *z, double *weight ) 
{
	int 	i ;
	double 	L[ 3 ] ;
	double 	*distance ;
	double 	total_weight ;

	distance =	new double [ base_number ] ;

	// Calculate distance to all the point and store in weight[] array. 
	for ( i = 0 ; i < base_number ; i ++ )	
	{
		L[ 0 ]	=	( x[ 0 ] - x[ i + 1 ] ) * ( x[ 0 ] - x[ i + 1 ] ) ;
		L[ 1 ]	=	( y[ 0 ] - y[ i + 1 ] ) * ( y[ 0 ] - y[ i + 1 ] ) ;
		L[ 2 ]	=	( z[ 0 ] - z[ i + 1 ] ) * ( z[ 0 ] - z[ i + 1 ] ) ;

		distance[ i ]	=	sqrt ( L[ 0 ] + L[ 1 ] + L[ 2 ] ) ;
		weight[ i ]		=	0. ;
	}

	// Check if the point is very close to one of the base point 
	for ( i = 0 ; i < base_number ; i ++ )	
	{
		if ( distance[ i ] < 1.e-25 )
		{
			weight[ i ]	=	1. ;
			i =	123456 ;
		}
	}

	// Inversed distance
	if ( i != 123456 )
	{
		total_weight =	0. ; 
		for ( i = 0 ; i < base_number ; i++ )	
		{
			weight[ i ]		=	1. / distance[ i ] ;
			total_weight	+=	weight[ i ] ;
		}

		for ( i = 0 ; i < base_number ; i++ )	
		{
			weight[ i ]	=	weight[ i ] / total_weight ;
		}		
	}

	delete [] distance ;
}
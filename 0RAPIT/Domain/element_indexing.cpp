#include <mpi.h>
#include "domain.h"
#include "element.h"
#include <iostream>
#include "sys_log.h"

using namespace std ;

void Domain::Element_indexing ()
{
	int 	i, j, k, m, mm ;
	int 	NodeNo, FaceNo, CellNo, MeshNo, LocalNo ;
	double	sum_1_r;
	vector< int > neighbordata ;

	global_node_number	=	Mesh.NodeNum ;
	global_face_number	=	Mesh.FaceNum ;
	global_cell_number	=	Mesh.CellNum ;

	local_node_number	=	Processor_MeshNode[ comm_rank ].size() ;
	local_face_number	=	Processor_MeshFace[ comm_rank ].size() ;
	local_cell_number	=	Processor_MeshCell[ comm_rank ].size() ;

	ghost_node_number_level_1	=	Processor_GhostMeshNode_lv1[ comm_rank ].size() ;
	ghost_node_number 			=	ghost_node_number_level_1 ;

	ghost_face_number_level_1	=	Processor_GhostMeshFace_lv1[ comm_rank ].size() ;
	ghost_face_number 			=	ghost_face_number_level_1 ;

	ghost_cell_number_level_1	=	Processor_GhostMeshCell_lv1[ comm_rank ].size() ;
	ghost_cell_number_level_2	=	Processor_GhostMeshCell_lv2[ comm_rank ].size() ;
	ghost_cell_number_level_3	=	Processor_GhostMeshCell_lv3[ comm_rank ].size() ;
	ghost_cell_number 			=	ghost_cell_number_level_1 + ghost_cell_number_level_2 + ghost_cell_number_level_3 ;

	node = boost::shared_array< Node > ( new Node [ local_node_number + ghost_node_number ] ) ;
	face = boost::shared_array< Face > ( new Face [ local_face_number + ghost_face_number ] ) ;
	cell = boost::shared_array< Cell > ( new Cell [ local_cell_number + ghost_cell_number ] ) ;

	// Node
	for ( i = 0 ; i < local_node_number + ghost_node_number ; i++ )
	{
		MeshNo =	LocalNode_MeshNodeNo[ i ] ;

		// set_id ( global, local )
		node[ i ].set_id ( MeshNo, i ) ;
		node[ i ].set_position( Mesh.Node_Position[ 0 ][ MeshNo ], Mesh.Node_Position[ 1 ][ MeshNo ], Mesh.Node_Position[ 2 ][ MeshNo ] ) ;

		// node-face relationship
		check_neighbor_data( &Mesh.Node_Face[ MeshNo ], MeshFace_LocalFaceNo.get(), &neighbordata ) ;
		node[ i ].set_face_number( neighbordata.size() ) ;
		for ( j = 0 ; j < node[ i ].face_number ; j++ )
		{
			LocalNo =	MeshFace_LocalFaceNo[ neighbordata[ j ] ] ;
			node[ i ].face[ j ]	=	&face[ LocalNo ] ;
		}

		// node-cell relationship
		if ( i < local_node_number )
		{
			node[ i ].set_cell_number( Mesh.Node_Cell[ MeshNo ].size() ) ;

			for ( j = 0 ; j < node[ i ].cell_number ; j++ )
			{
				LocalNo =	MeshCell_LocalCellNo[ Mesh.Node_Cell[ MeshNo ][ j ] ] ;
				node[ i ].cell[ j ]	=	&cell[ LocalNo ] ;
			}
		} else
		{
			check_neighbor_data( &Mesh.Node_Cell[ MeshNo ], MeshCell_LocalCellNo.get(), &neighbordata ) ;
			node[ i ].set_cell_number( neighbordata.size() ) ;

			for ( j = 0 ; j < node[ i ].cell_number ; j++ )
			{
				LocalNo =	MeshCell_LocalCellNo[ neighbordata[ j ] ] ;
				node[ i ].cell[ j ]	=	&cell[ LocalNo ] ;
			}
		}
	}

	// Face
	for ( i = 0 ; i < local_face_number + ghost_face_number ; i++ )
	{
		MeshNo =	LocalFace_MeshFaceNo[ i ] ;

		face[ i ].mpi_id	=	Mesh.Face_ProcessorNo[ MeshNo ] ;

		// set_id ( global, local )
		face[ i ].set_id ( MeshNo, i ) ;
		face[ i ].set_type( Mesh.Face_Type[ MeshNo ], Mesh.Face_Typename[ MeshNo ] ) ;
		face[ i ].set_position( Mesh.Face_Position[ 0 ][ MeshNo ], Mesh.Face_Position[ 1 ][ MeshNo ], Mesh.Face_Position[ 2 ][ MeshNo ] ) ;

		// face-node relationship
		face[ i ].set_node_number( Mesh.Face_Node[ MeshNo ].size() ) ;
		for ( j = 0 ; j < face[ i ].node_number ; j++ )
		{
			LocalNo =	MeshNode_LocalNodeNo[ Mesh.Face_Node[ MeshNo ][ j ] ] ;
			face[ i ].node[ j ]	=	&node[ LocalNo ] ;
		}

		// face-cell relationship
		face[ i ].set_cell_number( Mesh.Face_Cell[ MeshNo ].size() ) ;
		for ( j = 0 ; j < face[ i ].cell_number ; j++ )
		{
			LocalNo =	MeshCell_LocalCellNo[ Mesh.Face_Cell[ MeshNo ][ j ] ] ;
			face[ i ].cell[ j ]	=	&cell[ LocalNo ] ;			
		}
	}

	// cell
	for ( i = 0 ; i < local_cell_number + ghost_cell_number ; i++ )
	{
		MeshNo =	LocalCell_MeshCellNo[ i ] ;

		cell[ i ].mpi_id	=	Mesh.Cell_ProcessorNo[ MeshNo ] ;

		// set_id ( global, local, mesh )
		cell[ i ].set_id ( LocalCell_GlobalCellNo[ i ], i, MeshNo ) ;
		cell[ i ].set_type( Mesh.Cell_Type[ MeshNo ], Mesh.Cell_Typename[ MeshNo ] ) ;
		cell[ i ].set_volume( Mesh.Cell_Volume[ MeshNo ] ) ;
		cell[ i ].set_position( Mesh.Cell_Position[ 0 ][ MeshNo ], Mesh.Cell_Position[ 1 ][ MeshNo ], Mesh.Cell_Position[ 2 ][ MeshNo ] ) ;

		if ( i < local_cell_number + ghost_cell_number_level_1 )
		{
			// cell-node relationship
			cell[ i ].set_node_number( Mesh.Cell_Node[ MeshNo ].size() ) ;
			for ( j = 0 ; j < cell[ i ].node_number ; j++ )
			{
				LocalNo =	MeshNode_LocalNodeNo[ Mesh.Cell_Node[ MeshNo ][ j ] ] ;
				cell[ i ].node[ j ]	=	&node[ LocalNo ] ;
			}

			// cell-face relationship
			cell[ i ].set_face_number( Mesh.Cell_Face[ MeshNo ].size() ) ;
			for ( j = 0 ; j < cell[ i ].face_number ; j++ )
			{
				LocalNo =	MeshFace_LocalFaceNo[ Mesh.Cell_Face[ MeshNo ][ j ] ] ;
				cell[ i ].face[ j ]	=	&face[ LocalNo ] ;
			}

			// cell-cell relationship
			cell[ i ].set_cell_number( Mesh.Cell_Cell[ MeshNo ].size() ) ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				LocalNo =	MeshCell_LocalCellNo[ Mesh.Cell_Cell[ MeshNo ][ j ] ] ;
				cell[ i ].cell[ j ]	=	&cell[ LocalNo ] ;	
			}
		}
	}

	// Node
	for ( i = 0 ; i < local_node_number + ghost_node_number_level_1 ; i++ )
	{
		sum_1_r	=	0. ;
		for ( j = 0 ; j < node[ i ].cell_number ; j++ )
		{
			node[ i ].alpha[ j ]	=	1. / dist( node[ i ].r[ 0 ], node[ i ].r[ 1 ], node[ i ].r[ 2 ], node[ i ].cell[ j ]->r[ 0 ], node[ i ].cell[ j ]->r[ 1 ], node[ i ].cell[ j ]->r[ 2 ] ) ;
			sum_1_r	+=	node[ i ].alpha[ j ] ;
		}

		for ( j = 0 ; j < node[ i ].cell_number ; j++ )
			node[ i ].alpha[ j ]	=	node[ i ].alpha[ j ] / sum_1_r ;;
	}

	// face
	for ( i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
	{
		sum_1_r	=	0. ;
		for ( j = 0 ; j < face[ i ].cell_number ; j++ )
		{
			face[ i ].alpha[ j ] = 1. / dist( face[ i ].r[ 0 ], face[ i ].r[ 1 ], face[ i ].r[ 2 ], face[ i ].cell[ j ]->r[ 0 ], face[ i ].cell[ j ]->r[ 1 ], face[ i ].cell[ j ]->r[ 2 ] ) ;
			sum_1_r +=	face[ i ].alpha[ j ] ;
		}

		for ( j = 0 ; j < face[ i ].cell_number ; j++ )
			face[ i ].alpha[ j ]	=	face[ i ].alpha[ j ] / sum_1_r ;
	}
}

// 
double dist( double x1, double y1, double z1, double x2, double y2, double z2 )
{
	double	dx, dy, dz ;

	dx	=	x1 - x2 ;
	dy	=	y1 - y2 ;
	dz	=	z1 - z2 ;

	return	sqrt( dx * dx + dy * dy + dz * dz ) ;
}

// Return No. is MeshNo.
void check_neighbor_data( vector<int> *meshdata, int *table_mesh_local, vector<int> *neighbordata )
{
	neighbordata->clear() ;

	for ( int i = 0 ; i < meshdata->size() ; i++ )
	{
		if ( table_mesh_local[ (*meshdata)[ i ] ] >= 0 )
		{
			neighbordata->push_back( (*meshdata)[ i ] ) ;
		}
	}
}

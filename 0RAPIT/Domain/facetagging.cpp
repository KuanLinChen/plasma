#include <domain.h>
#include <iostream>
#include <set>
#include "sys_log.h"

using namespace std ;

void Domain::FaceTagging()
{
	bool flg = 0 ;
	int i, j ;

	BCFace_type_number =	Mesh.BC_TypeNum_forFace ;
	BCFace = boost::shared_array< vector<int> >  ( new vector<int> [ BCFace_type_number ] ) ;

	for ( i = 0 ; i < BCFace_type_number ; i++ )
	{
		BCFace_type.push_back( i ) ;
		BCFace_typename.push_back( Mesh.BC_Typename[ i ] ) ;
	}

	for ( i = 0 ; i < local_face_number + ghost_face_number ; i++ ) 
	{
		if ( face[ i ].type == 0 )
		{
			BulkFace.push_back( i ) ;
			NormalFace.push_back ( i ) ;

			if ( face[ i ].cell[ 0 ]->local_id >= local_cell_number || face[ i ].cell[ 1 ]->local_id >= local_cell_number ) 
			{
				MPIFace.push_back( i ) ;
			} 
		} else if ( face[ i ].type == BCFace_type_number )
		{
			InterFace.push_back( i ) ;
			NormalFace.push_back ( i ) ;

			if ( face[ i ].cell[ 0 ]->local_id >= local_cell_number || face[ i ].cell[ 1 ]->local_id >= local_cell_number ) 
			{
				MPIFace.push_back( i ) ;
			} 
		} else
		{
			BCFace[ face[ i ].type ].push_back( i ) ;
		}
	}

	if ( flg && comm_rank == 0 )
	{
		for ( i = 0 ; i < BCFace_type_number ; i++ )
		{
			cout << "Type = " << i << ", name = " << BCFace_typename[ i ] << endl ;
			for ( j = 0 ; j < BCFace[ i ].size() ; j++ )
			{
				cout << "\tF = " << BCFace[ i ][ j ] << ", id = " << face[ BCFace[ i ][ j ] ].id << endl ;
			}
		}
	}
}
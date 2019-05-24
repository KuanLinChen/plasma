#include "domain.h"

using namespace std ;

void Domain::CalculateFaceSign( ) 
{
	int 	i, j , k ;
	double 	r_cf_dot_A ;
	// sign of cell's face
	for ( i = 0 ; i < local_cell_number + ghost_cell_number_level_1 ; i++ )
	{
		for ( j = 0 ; j < cell[ i ].face_number ; j++ )
		{
			if ( cell[ i ].id == cell[ i ].face[ j ]->cell[ 0 ]->id )
			{
				cell[ i ].face_index_c0[ j ]	=	0 ;
				cell[ i ].face_index_c1[ j ]	=	1 ;
			} else
			{
				cell[ i ].face_index_c0[ j ]	=	1 ;
				cell[ i ].face_index_c1[ j ]	=	0 ;
			}

			r_cf_dot_A	=	0. ;

			for ( k = 0 ; k < dimension ; k++ )
			{
				r_cf_dot_A	+=	cell[ i ].face[ j ]->r_cf[ cell[ i ].face_index_c0[ j ] ][ k ] * cell[ i ].face[ j ]->A[ k ] ;
			}

			if ( r_cf_dot_A < 0. )
			{
				cell[ i ].face_sign[ j ]	=	- 1.0 ;
			} else
			{
				cell[ i ].face_sign[ j ]	=	1.0 ;

			}

			cell[ i ].face[ j ]->cell_sign[ cell[ i ].face_index_c0[ j ] ]	=	cell[ i ].face_sign[ j ] ;
			if ( cell[ i ].face[ j ]->cell_number > 1 )
				cell[ i ].face[ j ]->cell_sign[ 1 - cell[ i ].face_index_c0[ j ] ]	=	- cell[ i ].face_sign[ j ] ;
		}
	}
}
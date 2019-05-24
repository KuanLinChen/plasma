#include "domain.h"


void Domain::CalculateCellFaceAreaComponent( )
{
	int	i, j, k, c0_index ;

	for ( i = 0 ; i < local_cell_number + ghost_cell_number_level_1 ; i++ )
	{
		for ( j = 0 ; j < cell[ i ].face_number ; j++ )
		{
			c0_index	=	cell[ i ].face_index_c0[ j ] ;

			for ( k = 0 ; k < dimension ; k++ )
			{
				cell[ i ].A[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->A[ k ] ;
				cell[ i ].nA[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->nA[ k ] ;
			}
		}
	}

/*  // old version
	if ( dimension == 2 )
	{
		for ( i = 0 ; i < local_cell_number + ghost_cell_number_level_1 ; i++ )
		{
			for ( j = 0 ; j < cell[ i ].face_number ; j++ )
			{
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// New version
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				c0_index		=	cell[ i ].face_index_c0[ j ] ;

				for ( k = 0 ; k < dimension ; k++ )
				{
				//	cell[ i ]._xi[ j ][ k ]		=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->_xi[ k ] ;
				//	cell[ i ]._eta[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->_eta[ k ] ;
					cell[ i ].A[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->A[ k ] ;
					cell[ i ].nA[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->nA[ k ] ;
				}
				//cell[ i ].A_dot_xi[ j ]	=	cell[ i ].face_sign[ j ] * pow( -1.0, c0_index ) * cell[ i ].face[ j ]->A_dot_xi ;
				//cell[ i ].Jacobian[ j ]	=	cell[ i ].A_dot_xi[ j ] / cell[ i ].face[ j ]->dA ;
			}
		}
	} else
	{
		for ( i = 0 ; i < local_cell_number + ghost_cell_number_level_1 ; i++ )
		{
			for ( j = 0 ; j < cell[ i ].face_number ; j++ )
			{
				c0_index	=	cell[ i ].face_index_c0[ j ] ;

				for ( k = 0 ; k < dimension ; k++ )
				{
					//cell[ i ]._xi[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->_xi[ k ] ;
					//cell[ i ]._eta[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->_eta[ k ] ;
					//cell[ i ]._zeta[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->_zeta[ k ] ;
					cell[ i ].A[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->A[ k ] ;
					cell[ i ].nA[ j ][ k ]	=	cell[ i ].face_sign[ j ] * cell[ i ].face[ j ]->nA[ k ] ;
				}
				//cell[ i ].A_dot_xi[ j ]	=	cell[ i ].face_sign[ j ] * pow( -1.0, c0_index ) * cell[ i ].face[ j ]->A_dot_xi ;
				//cell[ i ].Jacobian[ j ]	=	cell[ i ].A_dot_xi[ j ] / cell[ i ].face[ j ]->dA ;
			}
		}
	}
*/
}

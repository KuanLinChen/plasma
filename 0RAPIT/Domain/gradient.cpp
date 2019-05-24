#include "domain.h"
#include <cmath>

void Domain::Cell_gradient( double *variable_at_face, double **gradient_variable )
{
	int	i, j, k, faceid ;

	if ( cylindrical_y == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					faceid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ faceid ].x > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ faceid ].x / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].x / 2. / M_PI ) ;
			}

		}
	}  else if ( cylindrical_x == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					faceid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ faceid ].y > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ faceid ].y / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].y / 2. / M_PI ) ;
			}
		}
	} else
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
					gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
				gradient_variable[ j ][ i ]	/=	cell[ i ].volume ;
			}
		}
	}
}

void Domain::Cell_gradient( double *variable_at_face, double *gradient_variable_x, double *gradient_variable_y, double *gradient_variable_z )
{
	int	i, j, k, fid ;

	double * gradient_variable[ 3 ] ;

	gradient_variable[ 0 ]	=	gradient_variable_x ;
	gradient_variable[ 1 ]	=	gradient_variable_y ;
	gradient_variable[ 2 ]	=	gradient_variable_z ;

	if ( cylindrical_y == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					fid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ fid ].x > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ fid ].x / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
					//if (  i == 873 || i == 875 || i == 874 )
					//{
					//	if ( j == 1 )
					//	cout << i << " k " << k << " v " << variable_at_face[ cell[ i ].face[ k ]->local_id ]  << " A " <<  cell[ i ].A[ k ][ j ] << " fid "<< cell[ i ].face[k]->local_id << " dx " << gradient_variable[ j ][ i ]	 << endl;
					//}
				}
				//if (  i == 873 || i == 875 || i == 874 )
				//{
				//	cout  << "b " << gradient_variable[ 0 ][ i ] << "  "<<  cell[ i ].volume / cell[i].x / 2. / M_PI<< endl;
				//	cout  << gradient_variable[ 1 ][ i ] << "  "<<  cell[ i ].volume / cell[i].x / 2. / M_PI<< endl;
				//}

				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].x / 2. / M_PI ) ;
				//if (  i == 873 || i == 875 || i == 874 )
				//{
				//	cout  << gradient_variable[ 0 ][ i ] << "  "<<  cell[ i ].volume / cell[i].x / 2. / M_PI<< endl;
				//	cout  << gradient_variable[ 1 ][ i ] << "  "<<  cell[ i ].volume / cell[i].x / 2. / M_PI<< endl;
				//}
			}

		}
	}  else if ( cylindrical_x == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					fid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ fid ].y > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ fid ].y / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].y / 2. / M_PI ) ;
			}
		}
	} else
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;

					//if ( i == 875 || i == 896 )
					//{
					//	cout << "cell " << i << " f " << j << " face_x "<< cell[ i ].face[ k ]->x << " dv " << gradient_variable[ j ][ i ] << " A " << cell[ i ].A[ k ][ j ] << endl;
					//}
				}
				gradient_variable[ j ][ i ]	/=	cell[ i ].volume ;
			}
		}
	}
}

void Domain::Cell_gradient_value( double *variable, double *variable_at_face, double **gradient_variable, double *min_value, double *max_value )
{
	int	i, j, k, cellid, faceid ;

	if ( cylindrical_y == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					faceid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ faceid ].x > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ faceid ].x / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].x / 2. / M_PI ) ;
			}

			min_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			max_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				cellid	=	cell[ i ].cell[ j ]->local_id ;

				// min
				if ( variable[ cellid ] < min_value[ i ] )
					min_value[ i ]	=	variable[ cellid ] ;
				// max
				if ( variable[ cellid ] > max_value[ i ] )
					max_value[ i ]	=	variable[ cellid ] ;
			}
		}
	}  else if ( cylindrical_x == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					faceid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ faceid ].y > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ faceid ].y / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].y / 2. / M_PI ) ;
			}

			min_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			max_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				cellid	=	cell[ i ].cell[ j ]->local_id ;

				// min
				if ( variable[ cellid ] < min_value[ i ] )
					min_value[ i ]	=	variable[ cellid ] ;
				// max
				if ( variable[ cellid ] > max_value[ i ] )
					max_value[ i ]	=	variable[ cellid ] ;
			}
		}
	} else
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
					gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
				gradient_variable[ j ][ i ]	/=	cell[ i ].volume ;
			}

			min_value[ i ]	=	variable[ i ] ; //variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			max_value[ i ]	=	variable[ i ] ; //variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				cellid	=	cell[ i ].cell[ j ]->local_id ;

				// min
				if ( variable[ cellid ] < min_value[ i ] )
					min_value[ i ]	=	variable[ cellid ] ;
				// max
				if ( variable[ cellid ] > max_value[ i ] )
					max_value[ i ]	=	variable[ cellid ] ;
			}
		}
	}
}

void Domain::Cell_gradient_value( double *variable, double *variable_at_face, double *gradient_variable_x, double *gradient_variable_y, double *gradient_variable_z, double *min_value, double *max_value )
{
	int	i, j, k, cellid, faceid ;
	double	*gradient_variable[ 3 ] ;

	gradient_variable[ 0 ]	=	gradient_variable_x ;
	gradient_variable[ 1 ]	=	gradient_variable_y ;
	gradient_variable[ 2 ]	=	gradient_variable_z ;

	if ( cylindrical_y == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					faceid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ faceid ].x > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ faceid ].x / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].x / 2. / M_PI ) ;
			}

			min_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			max_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				cellid	=	cell[ i ].cell[ j ]->local_id ;

				// min
				if ( variable[ cellid ] < min_value[ i ] )
					min_value[ i ]	=	variable[ cellid ] ;
				// max
				if ( variable[ cellid ] > max_value[ i ] )
					max_value[ i ]	=	variable[ cellid ] ;
			}
		}
	}  else if ( cylindrical_x == 1 )
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
				{
					faceid	=	cell[ i ].face[ k ]->local_id ;
					if ( face[ faceid ].y > 1.e-10  )
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * ( cell[ i ].A[ k ][ j ] / face[ faceid ].y / 2. / M_PI ) ;
					} else
					{
						gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
					}
				}
				gradient_variable[ j ][ i ]	/=	( cell[ i ].volume / cell[ i ].y / 2. / M_PI ) ;
			}

			min_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			max_value[ i ]	=	variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				cellid	=	cell[ i ].cell[ j ]->local_id ;

				// min
				if ( variable[ cellid ] < min_value[ i ] )
					min_value[ i ]	=	variable[ cellid ] ;
				// max
				if ( variable[ cellid ] > max_value[ i ] )
					max_value[ i ]	=	variable[ cellid ] ;
			}
		}
	} else
	{
		for ( i = 0 ; i < local_cell_number ; i++ )
		{
			for ( j = 0 ; j < dimension ; j++ )
			{
				gradient_variable[ j ][ i ]	=	0. ;
				for ( k = 0 ; k < cell[ i ].face_number ; k++ )
					gradient_variable[ j ][ i ]	+=	variable_at_face[ cell[ i ].face[ k ]->local_id ] * cell[ i ].A[ k ][ j ] ;
				gradient_variable[ j ][ i ]	/=	cell[ i ].volume ;
			}

			min_value[ i ]	=	variable[ i ] ; //variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			max_value[ i ]	=	variable[ i ] ; //variable[ cell[ i ].cell[ 0 ]->local_id ] ;
			for ( j = 0 ; j < cell[ i ].cell_number ; j++ )
			{
				cellid	=	cell[ i ].cell[ j ]->local_id ;

				// min
				if ( variable[ cellid ] < min_value[ i ] )
					min_value[ i ]	=	variable[ cellid ] ;
				// max
				if ( variable[ cellid ] > max_value[ i ] )
					max_value[ i ]	=	variable[ cellid ] ;
			}
		}
	}
}
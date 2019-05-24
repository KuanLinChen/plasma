#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include "preprocessing.h"
#include "sys_log.h"

void Structured_mesh( int dimension, double scale, int nx, int ny, int nz, double dx, double dy, double dz, int NodeNum, int CellNum, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, vector<string> *BC_Typename, int *BCFace_Type, string *BCFace_Typename, int *Cell_FaceNum, vector<int> *Cell_Node, vector<int> *Node_Cell, vector<int> *BCFace_Node, vector<int> *Node_BCFace )
{
	int i, j, k, n ;
	int NodeNo, FaceNo[ 2 ], CellNo ;
	vector<string>	_BC_Typename ;
	map< string, int > Typename_Type ;
	map< string, int > check_Typename ;

	_BC_Typename.push_back( "Bulk" ) ;
	check_Typename[ "Bulk" ]	=	1 ;

	(*BC_TypeNum) =	1 ;
	for ( i = 0 ; i < BC_Typename->size() ; i++ )
	{
		if ( check_Typename[ (*BC_Typename)[ i ] ] == 0 )
		{
			_BC_Typename.push_back( (*BC_Typename)[ i ] ) ;

			Typename_Type[ (*BC_Typename)[ i ] ]	=	( *BC_TypeNum ) ;
			check_Typename[ (*BC_Typename)[ i ]	]	=	1 ;

			(*BC_TypeNum)++ ;
		}
	}

	if ( check_Typename[ "Interface" ] == 0 )
	{
		_BC_Typename.push_back( "Interface" ) ;

		Typename_Type[ "Interface" ]	=	( *BC_TypeNum ) ;
		(*BC_TypeNum)++ ;
	}

	if ( dimension == 2 )
	{
		Log().TagDump( logLEVEL4 )  << "\tNx = " << nx << ", Ny = " << ny << ", Dx = " << dx << ", Dy = " << dy ;
		Log().TagDump( logLEVEL4 )  << "\tBoundary:" ;
		for ( i = 0 ; i < 4 ; i++ )
			Log().TagDump( logLEVEL4 ) << "\t\t" << (*BC_Typename)[ i ] ;

		for ( j = 0 ; j < ny ; j++)
		{
			for ( i = 0 ; i < nx ; i++ )
			{
				NodeNo				=	i + j * nx ;

				Node_x[ NodeNo ]	=	double ( i ) * dx * scale ;
				Node_y[ NodeNo ]	=	double ( j ) * dy * scale ;
				Node_z[ NodeNo ]	=	0.0 ;
			}
		}

		for ( j = 0 ; j < ( ny - 1 ) ; j++ )
		{
			for ( i = 0 ; i < ( nx - 1 ) ; i++ )
			{
				CellNo	=	i + j * ( nx -1 ) ;

				Cell_FaceNum[ CellNo ]	=	4 ;

				Cell_Node[ CellNo ].push_back( i + j * nx ) ;
				Cell_Node[ CellNo ].push_back( i + j * nx + 1 ) ;
				Cell_Node[ CellNo ].push_back( i + ( j + 1 ) * nx + 1 ) ;
				Cell_Node[ CellNo ].push_back( i + ( j + 1 ) * nx ) ;

				Cell_x[ CellNo ]	= 0.0 ;
				Cell_y[ CellNo ]	= 0.0 ;
				Cell_z[ CellNo ]	= 0.0 ;
				for ( n = 0 ; n < 4 ; n++ )
				{
					NodeNo =	Cell_Node[ CellNo ][ n ] ;

					Node_Cell[ NodeNo ].push_back( CellNo ) ;

					Cell_x[ CellNo ]	+=	Node_x[ NodeNo ] ;
					Cell_y[ CellNo ]	+=	Node_y[ NodeNo ] ;			
				}

				Cell_x[ CellNo ]	/=	4 ;
				Cell_y[ CellNo ]	/=	4 ;
			}
		}

		FaceNo[ 0 ]	=	0 ;
		// BCFace_1 & BCFace_2, x = 0. & x = xmax
		for ( j = 0 ; j < ( ny - 1 ) ; j++ )
		{
			// BCFace_1, x = 0.0
			BCFace_Typename[ FaceNo[ 0 ] ]	=	(*BC_Typename)[ 0 ] ;
			BCFace_Type[ FaceNo[ 0 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 0 ] ] ] ;

			NodeNo	=	j * nx ;
			BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

			NodeNo	=	( j + 1 ) * nx ;
			BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

			// BCFace_2, x = xmax
			FaceNo[ 1 ]	=	FaceNo[ 0 ] + ( ny - 1 ) ;
			BCFace_Typename[ FaceNo[ 1 ] ]	=	(*BC_Typename)[ 1 ] ;
			BCFace_Type[ FaceNo[ 1 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 1 ] ] ] ;

			NodeNo	=	( j + 1 ) * nx - 1 ;
			BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

			NodeNo	=	( j + 2 ) * nx - 1 ;
			BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

			FaceNo[ 0 ]++ ;
		}

		FaceNo[ 0 ]	=	2 * ( ny - 1 ) ;
		// BCFace_3 & BCFace_4, y = 0. & y = ymax
		for ( i = 0 ; i < ( nx - 1 ) ; i++ )
		{
			// BCFace_3, y = 0.0
			BCFace_Typename[ FaceNo[ 0 ] ]	=	(*BC_Typename)[ 2 ] ;
			BCFace_Type[ FaceNo[ 0 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 0 ] ] ] ;

			NodeNo	=	i ;
			BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

			NodeNo	=	i + 1 ;
			BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

			// BCFace_4, y = ymax
			FaceNo[ 1 ]	=	FaceNo[ 0 ] + nx - 1 ;
			BCFace_Typename[ FaceNo[ 1 ] ]	=	(*BC_Typename)[ 3 ] ;
			BCFace_Type[ FaceNo[ 1 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 1 ] ] ] ;

			NodeNo	=	( ny - 1 ) * nx + i ;
			BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] + nx -1 ) ;

			NodeNo	=	( ny - 1 ) * nx + i + 1 ;
			BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
			Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

			FaceNo[ 0 ]++ ;
		}
	} else
	{
		Log().TagDump( logLEVEL4 )  << "\tNx = " << nx << ", Ny = " << ny << ", Nz = " << nz << ", Dx = " << dx << ", Dy = " << dy << ", Dz = " << dz ;
		Log().TagDump( logLEVEL4 )  << "\tBoundary:" ;
		for ( i = 0 ; i < 6 ; i++ )
			Log().TagDump( logLEVEL4 ) << "\t\t" << (*BC_Typename)[ i ] ;

		for ( k = 0 ; k < nz ; k++ )
		{
			for ( j = 0 ; j < ny ; j++ )
			{
				for ( i = 0 ; i < nx ; i++ )
				{
					NodeNo	=	i + j * nx + k * nx * ny ;

					Node_x[ NodeNo ]	=	double ( i ) * dx * scale ;
					Node_y[ NodeNo ]	=	double ( j ) * dy * scale ;
					Node_z[ NodeNo ]	=	double ( k ) * dz * scale ;
				}
			}
		}

		for ( k = 0 ; k < ( nz - 1 ) ; k++ )
		{
			for ( j = 0 ; j < ( ny - 1 ) ; j++ )
			{
				for ( i = 0 ; i < ( nx - 1 ) ; i++ )
				{
					CellNo	=	i + j * ( nx - 1 ) + k * ( nx - 1 ) * ( ny - 1 ) ;

					Cell_FaceNum[ CellNo ]	=	6 ;

					Cell_Node[ CellNo ].push_back( i + j * nx + k * nx * ny ) ;
					Cell_Node[ CellNo ].push_back( i + j * nx + k * nx * ny + 1 ) ;
					Cell_Node[ CellNo ].push_back( i + ( j + 1 ) * nx + k * nx * ny + 1 ) ;
					Cell_Node[ CellNo ].push_back( i + ( j + 1 ) * nx + k * nx * ny ) ;
					Cell_Node[ CellNo ].push_back( i + j * nx + ( k + 1 ) * nx * ny ) ;
					Cell_Node[ CellNo ].push_back( i + j * nx + ( k + 1 ) * nx * ny + 1 ) ;
					Cell_Node[ CellNo ].push_back( i + ( j + 1 ) * nx + ( k + 1 ) * nx * ny + 1 ) ;
					Cell_Node[ CellNo ].push_back( i + ( j + 1 ) * nx + ( k + 1 ) * nx * ny ) ;

					Cell_x[ CellNo ]	= 0. ;
					Cell_y[ CellNo ]	= 0. ;
					Cell_z[ CellNo ]	= 0. ;
					for ( n = 0 ; n < 8 ; n++ )
					{
						NodeNo =	Cell_Node[ CellNo ][ n ] ;

						Node_Cell[ NodeNo ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo ] ;
						Cell_z[ CellNo ]	+=	Node_z[ NodeNo ] ;
					}

					Cell_x[ CellNo ]	/=	8 ;
					Cell_y[ CellNo ]	/=	8 ;
					Cell_z[ CellNo ]	/=	8 ;					
				}
			}
		}

		// BCFace_1 & BCFace_2, x = 0. & x = xmax
		FaceNo[ 0 ]	=	0 ;
		for ( k = 0 ; k < ( nz - 1 ) ; k++ )
		{
			for ( j = 0 ; j < ( ny - 1 ) ; j++ )
			{
				// BCFace_1, x = 0.0
				BCFace_Typename[ FaceNo[ 0 ] ]	=	(*BC_Typename)[ 0 ] ;
				BCFace_Type[ FaceNo[ 0 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 0 ] ] ] ;

				NodeNo	=	k * nx * ny + j * nx ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	k * nx * ny + ( j + 1 ) * nx ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + ( j + 1 ) * nx ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + j * nx ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				// BCFace_2, x = xmax
				FaceNo[ 1 ]	=	FaceNo[ 0 ] + ( ny -1 ) * ( nz -1 ) ;
				BCFace_Typename[ FaceNo[ 1 ] ]	=	(*BC_Typename)[ 1 ] ;
				BCFace_Type[ FaceNo[ 1 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 1 ] ] ] ;

				NodeNo	=	k * nx * ny + ( j + 1 ) * nx - 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	k * nx * ny + ( j + 2 ) * nx - 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + ( j + 2 ) * nx - 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + ( j + 1 ) * nx - 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				FaceNo[ 0 ]++ ;
			}
		}

		// BCFace_3 & BCFace_4, y = 0. & y = ymax
		FaceNo[ 0 ]	=	2 * ( ny - 1 ) * ( nz -1 ) ;
		for ( k = 0 ; k < ( nz - 1 ) ; k++ )
		{
			for ( i = 0 ; i < ( nx - 1 ) ; i++ )
			{
				// BCFace_3, y = 0.0
				BCFace_Typename[ FaceNo[ 0 ] ]	=	(*BC_Typename)[ 2 ] ;
				BCFace_Type[ FaceNo[ 0 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 0 ] ] ] ;

				NodeNo	=	k * nx * ny + i ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	k * nx * ny + i + 1 ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + i + 1 ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + i ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				// BCFace_2, x = xmax
				FaceNo[ 1 ]	=	FaceNo[ 0 ] + ( nx -1 ) * ( nz -1 ) ;
				BCFace_Typename[ FaceNo[ 1 ] ]	=	(*BC_Typename)[ 3 ] ;
				BCFace_Type[ FaceNo[ 1 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 1 ] ] ] ;

				NodeNo	=	k * nx * ny + ( ny - 1 ) * nx + i ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	k * nx * ny + ( ny - 1 ) * nx + i + 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + (  ny - 1  ) * nx + i + 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( k + 1 ) * nx * ny + ( ny - 1 ) * nx + i ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				FaceNo[ 0 ]++ ;
			}
		}

		// BCFace_5 & BCFace_5, z = 0. & z = zmax
		FaceNo[ 0 ]	=	2 * ( ny - 1 ) * ( nz -1 ) + 2 * ( nx -1 ) * ( nz -1 ) ;
		for ( j = 0 ; j < ( ny - 1 ) ; j++ )
		{
			for ( i = 0 ; i < ( nx - 1 ) ; i++ )
			{				
				// BCFace_5, z = 0.0
				BCFace_Typename[ FaceNo[ 0 ] ]	=	(*BC_Typename)[ 4 ] ;
				BCFace_Type[ FaceNo[ 0 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 0 ] ] ] ;

				NodeNo	=	j * nx + i ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	j * nx + i + 1 ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	( j + 1 ) * nx + i + 1 ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				NodeNo	=	( j + 1 ) * nx + i ;
				BCFace_Node[ FaceNo[ 0 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 0 ] ) ;

				// BCFace_6, z = zmax
				FaceNo[ 1 ]	=	FaceNo[ 0 ] + ( nx -1 ) * ( ny -1 ) ;
				BCFace_Typename[ FaceNo[ 1 ] ]	=	(*BC_Typename)[ 5 ] ;
				BCFace_Type[ FaceNo[ 1 ] ]		=	Typename_Type[ BCFace_Typename[ FaceNo[ 1 ] ] ] ;

				NodeNo	=	( nz - 1 ) * nx * ny + j * nx + i ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( nz - 1 ) * nx * ny + j * nx + i + 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( nz - 1 ) * nx * ny + ( j + 1 ) * nx + i + 1 ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				NodeNo	=	( nz - 1 ) * nx * ny + ( j + 1 ) * nx + i ;
				BCFace_Node[ FaceNo[ 1 ] ].push_back( NodeNo ) ;
				Node_BCFace[ NodeNo ].push_back( FaceNo[ 1 ] ) ;

				FaceNo[ 0 ]++ ;
			}
		}
	}

	(*BC_Typename).clear() ;
	for ( i = 0 ; i < (*BC_TypeNum) ; i++ )
	{
		(*BC_Typename).push_back( _BC_Typename[ i ] ) ;
	}
}

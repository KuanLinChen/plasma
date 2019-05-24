#include <iostream>
#include "domain.h"

bool tol ( double x1, double y1, double z1, double x2, double y2, double z2 ) 
{
	double	L, l1, l2, l3 ;
	l1	=	x1 - x2 ;
	l2	=	y1 - y2 ;
	l3	=	z1 - z2 ;
	L	=	sqrt ( l1 * l1 + l2 * l2 + l3 * l3 ) ;

	if ( L < 1.0e-15 ) 
		return true ;
	else 
		return false ;
}

int main( int argc, char * argv[] )
{
	int	i, j, k ;
	int	mpi_id, mpi_size ;
	double	x, y, z ;

	MPI_Init( &argc, &argv ) ;
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_size ) ;
	MPI_Comm_rank ( MPI_COMM_WORLD, &mpi_id ) ;

	if ( argc != 2 ) 
	{
		cout << "run ./doamin_test <mesh_file>" << endl ;
		return 0 ;
	}

	Domain d ( argv[ 1 ] ) ;

	// Test cell-node
	for ( i = 0 ; i < d.local_cell_number ; i++ )
	{
		//cout << d.cell[ i ].r[ 0 ] << "\t" <<  d.cell[ i ].r[ 1 ] << "\t" <<  d.cell[ i ].r[ 2 ] << endl ;
		x = y = z = 0. ;
		for ( j = 0 ; j < d.cell[ i ].node_number ; j++ )
		{
			x	+=	d.cell[ i ].node[ j ]->x ;
			y	+=	d.cell[ i ].node[ j ]->y ;
			z	+=	d.cell[ i ].node[ j ]->z ;
		}
		x	=	x / d.cell[ i ].node_number ;
		y	=	y / d.cell[ i ].node_number ;
		z	=	z / d.cell[ i ].node_number ;
		if ( !tol ( x, y, z, d.cell[ i ].x, d.cell[ i ].y , d.cell[ i ].z ) )
			cout << "Something thing wrong at cell[" << i << "]" << endl ;
		//cout << d.cell[ i ].r[ 0 ] << "\t" << x << "\t" <<  d.cell[ i ].r[ 1 ] << "\t" << y << "\t" <<  d.cell[ i ].r[ 2 ] << "\t" << z << endl ;
	}
	cout << "[" << mpi_id << "] Pass cell-node relation.." << endl ;

	// Test face-node
	for (  i = 0 ; i < d.local_face_number ; i++ )
	{
		x = y = z = 0. ;
		for ( j = 0 ; j < d.face[ i ].node_number ; j++ )
		{
			x	+=	d.face[ i ].node[ j ]->x ;
			y	+=	d.face[ i ].node[ j ]->y ;
			z	+=	d.face[ i ].node[ j ]->z ;
		}
		x	=	x / d.face[ i ].node_number ;
		y	=	y / d.face[ i ].node_number ;
		z	=	z / d.face[ i ].node_number ;
		if ( !tol ( x, y, z, d.face[ i ].x, d.face[ i ].y, d.face[ i ].z ) )
			cout << "Something thing wrong at face[" << i << "]" << endl ;  
		//cout << d.face [ i ].r[ 0 ] << "\t"  << x << "\t" <<  d.face[ i ].r[ 1 ] << "\t"  << y << "\t" <<  d.face[ i ].r[ 2 ] << "\t" << z << endl ;
	}
	cout << "[" << mpi_id << "] Pass face-node relation.." << endl ;

	// Test cell-face
	for ( i = 0 ; i < d.local_cell_number ; i++ )
	{
		//cout << d.cell[ i ].r[ 0 ] << "\t" <<  d.cell[ i ].r[ 1 ] << "\t" <<  d.cell[ i ].r[ 2 ] << endl ;
		x = y = z = 0. ;
		for ( j = 0 ; j < d.cell[ i ].face_number ; j++ )
		{
			x	+=	d.cell[ i ].face[ j ]->x ;
			y	+=	d.cell[ i ].face[ j ]->y ;
			z	+=	d.cell[ i ].face[ j ]->z ;
		}
		x	=	x / d.cell[ i ].face_number ;
		y	=	y / d.cell[ i ].face_number ;
		z	=	z / d.cell[ i ].face_number ;
		if ( !tol ( x, y, z, d.cell[ i ].x, d.cell[ i ].y, d.cell[ i ].z ) )
			cout << "Something thing wrong at cell[" << i << "]" << endl ;
		//cout << d.cell[ i ].r[ 0 ] << "\t"  << x << "\t" <<  d.cell[ i ].r[ 1 ] << "\t"  << y << "\t" <<  d.cell[ i ].r[ 2 ] << "\t"  << z << endl ;
	}
	cout << "["<< mpi_id << "] Pass cell-face relation.." << endl ;
}

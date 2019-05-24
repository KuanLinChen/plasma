#include <string>
#include "element.h"
//#include "domain.h"

using namespace std ;

/*! @file element.cpp The elements of Domain, including Node, Face, and Cell.

The design of the locker is used for the operator that need no updating. 

*/

/*! \brief Node constructor. 
*/
Node::Node()
{
	id			=	-1 ;
	local_id	=	-1 ;
	type		=	-1 ; 
	
	face_number	=	-1 ;
	cell_number	=	-1 ;
	x			=	0. ;
	y			=	0. ;
	z			=	0. ;

	r[ 0 ]		=	0. ;
	r[ 1 ]		=	0. ;
	r[ 2 ]		=	0. ;
}

/*! \brief Node destructor. 
*/
Node::~Node()
{
}


/*! \brief Set the local and global id. Since we have global ID, we can easily know the position from pre_Node.
*/
void Node::set_id ( int g, int l ) //( global, local, mesh )
{
	id			=	g ;
	local_id	=	l ;
}

void Node::set_position ( double p1, double p2, double p3 ) 
{
	x	=	p1 ;
	y	=	p2 ;
	z	=	p3 ;

	r[ 0 ]	=	p1 ;
	r[ 1 ]	=	p2 ;
	r[ 2 ]	=	p3 ;
}

void Node::set_face_number ( int size )
{
	face_number = size ;
	face = boost::shared_array < Face * > ( new Face* [ size ] ) ;
}

void Node::set_cell_number ( int size )
{
	cell_number = size ;
	cell 	= boost::shared_array < Cell * > ( new Cell* [ size ] ) ;
	alpha 	= boost::shared_array < double > ( new double [ size ] ) ;
	for ( int i = 0 ; i < size ; i++ )
		alpha[ i ] =	0. ;
}

Face::Face()
{
	mpi_id 		=	-1 ;
	id			=	-1 ;
	local_id	=	-1 ; 
	type		=	-1 ;

	x		=	0. ;
	y		=	0. ;
	z		=	0. ; 
	dx		=	0. ;
	dy		=	0. ;
	dz		=	0. ;	

	r[ 0 ]	=	0. ;
	r[ 1 ]	=	0. ;
	r[ 2 ]	=	0. ;
}

Face::~Face()
{
}


/*! \brief Set the local and global id. Since we have global ID, we can easily know the position from pre_Node.
*/
void Face::set_id ( int g, int l ) //( global, local )
{
	id			=	g ;
	local_id	=	l ;
}

void Face::set_position ( double p1, double p2, double p3 ) 
{
	x	=	p1 ;
	y	=	p2 ;
	z	=	p3 ;

	r[ 0 ]	=	p1 ;
	r[ 1 ]	=	p2 ;
	r[ 2 ]	=	p3 ;
}

void Face::set_type ( int t, string tn ) 
{
	type		=	t ;
	Typename	=	tn ; 
}

void Face::set_node_number ( int size )  
{
	node_number	=	size ;
	node		=	boost::shared_array < Node * > ( new Node* [ size ] ) ;
} 

void Face::set_cell_number ( int size ) 
{ 
	cell_number	=	size ;
	cell		=	boost::shared_array < Cell * > ( new Cell* [ size ] ) ; 
	cell_sign	=	boost::shared_array < double > ( new double [ size ] ) ; 
	alpha		=	boost::shared_array < double > ( new double [ size ] ) ; 
	dr			=	boost::shared_array < double > ( new double [ size ] ) ; 
	dr_x		=	boost::shared_array < double > ( new double [ size ] ) ; 
	dr_y		=	boost::shared_array < double > ( new double [ size ] ) ; 

	for ( int i = 0 ; i < size ; i++ )
	{
		cell_sign[ i ]	=	0. ;
		alpha[ i ]		=	0. ;
		dr[ i ]			=	0. ;
		dr_x[ i ]		=	0. ;
		dr_y[ i ]		=	0. ;
	}
}


Cell::Cell()
{
	id			=	-1 ;
	local_id 	=	-1 ;
	type		=	-1 ; 
	mpi_id		=	-1 ;

	x			=	0. ;
	y			=	0. ;
	z			=	0. ;
	volume		=	0. ;
	
	node_number	=	-1 ;
	face_number	=	-1 ;
	cell_number	=	-1 ;
}

Cell::~Cell()
{
}

/*! \brief Connect the Cell object with the system MeshCell.
*/
void Cell::set_id ( int g, int l, int m )  //( global, local, mesh )
{
	id			=	g ;
	local_id	=	l ;
	mesh_id 	=	m ;
}

void Cell::set_type ( int t, string tn ) 
{
	type		=	t ;
	Typename	=	tn ; 
}

void Cell::set_position ( double p1, double p2, double p3 ) 
{
	x	=	p1 ;
	y	=	p2 ;
	z	=	p3 ;

	r[ 0 ]	=	p1 ;
	r[ 1 ]	=	p2 ;
	r[ 2 ]	=	p3 ;
}

void Cell::set_volume ( double v )
{
	volume	=	v ;
}

void Cell::set_node_number ( int size ) 
{
	node_number	=	size ;
	node		=	boost::shared_array < Node * > ( new Node* [ size ] ) ; 
} 

void Cell::set_cell_number ( int size ) 
{ 
	cell_number	=	size ;
	cell		=	boost::shared_array < Cell * > ( new Cell* [ size ] ) ; 
}

void Cell::set_face_number ( int size ) 
{ 
	face_number 	=	size ;
	face			=	boost::shared_array < Face * > ( new Face* [ size ] ) ; 
	dr				=	boost::shared_array < double > ( new double [ size ] ) ; 
	face_sign		=	boost::shared_array < double > ( new double [ size ] ) ;
	face_index_c0	=	boost::shared_array < int > ( new int [ size ] ) ; 
	face_index_c1	=	boost::shared_array < int > ( new int [ size ] ) ; 

	//Jacobian	=	boost::shared_array < double  > ( new double [ size ] ) ; 
	//x_xi		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//y_xi		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//z_x		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//x_eta		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//y_eta		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//z_eta		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//x_zeta	=	boost::shared_array < double > ( new double [ size ] ) ; 
	//y_zeta	=	boost::shared_array < double > ( new double [ size ] ) ; 
	//z_zeta	=	boost::shared_array < double > ( new double [ size ] ) ; 
	Ax			=	boost::shared_array < double > ( new double [ size ] ) ; 
	Ay			=	boost::shared_array < double > ( new double [ size ] ) ; 
	Az			=	boost::shared_array < double > ( new double [ size ] ) ; 
	//A_xi		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//A_eta		=	boost::shared_array < double > ( new double [ size ] ) ; 
	//A_zeta	=	boost::shared_array < double > ( new double [ size ] ) ; 
	nx			=	boost::shared_array < double > ( new double [ size ] ) ; 
	ny			=	boost::shared_array < double > ( new double [ size ] ) ; 
	nz			=	boost::shared_array < double > ( new double [ size ] ) ; 
	//xi_coefficient	=	boost::shared_array < double > ( new double [ size ] ) ; 
	//_xi		=	boost::multi_array < double , 2 > ( (boost::extents[ size ][ 3 ] ) ) ;
	//_eta		=	boost::multi_array < double , 2 > ( (boost::extents[ size ][ 3 ] ) ) ;
	//_zeta		=	boost::multi_array < double , 2 > ( (boost::extents[ size ][ 3 ] ) ) ;
	A			=	boost::shared_array < vector<double> > ( new vector <double> [ size ] ) ;
	nA			=	boost::shared_array < vector<double> > ( new vector <double> [ size ] ) ;

	for ( int i = 0 ; i < size ; i++ )
	{
		A[ i ].reserve( 3 ) ;
		nA[ i ].reserve( 3 ) ;
		for ( int j = 0 ; j < 3 ; j++ )
		{
			A[ i ].push_back ( 0. ) ;
			nA[ i ].push_back ( 0. ) ;
		}

		dr[ i ]				=	0. ;
		face_sign[ i ]		=	0. ;
		face_index_c0[ i ]	=	0 ;
		face_index_c1[ i ]	=	0 ;
	}
	//A_dot_xi	=	boost::shared_array < double > ( new double [ size ] ) ; 
} 
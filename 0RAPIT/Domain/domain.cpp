#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <string.h>
#include <mpi.h>
#include <metis.h>
#include <algorithm>
#include "domain.h"
//#include "tecplot_output.h"
#include "uuid_functions.h"
#include "main.h"
#include "sys_log.h"

using namespace std ;

set<string> Domain::ID_bank ;

/*
string analysor( string s )
{
	s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() ) ;
	return s ;
}
*/
Domain::Domain( )
{
	parent_comm		=	MPI_COMM_WORLD ;
	comm			=	MPI_COMM_WORLD ;
	MPI_Comm_size ( comm, &comm_size ) ;
	MPI_Comm_rank ( comm, &comm_rank ) ;
	world_comm_size		=	comm_size ;
	world_comm_rank		=	comm_rank ;
	ID			=	generate_id() ;
	ID_bank.insert( ID ) ;

	flag_geometry				=	false ;
	flag_boundaryfile			=	false ;
	flag_meshfile				=	false ;
	flag_mesh_scale				=	false ;
	flag_structured_mesh		=	false ;

	dimension 					=	0 ;
	cylindrical_x				=	0 ;
	cylindrical_y				=	0 ;
	global_node_number			=	0 ;
	global_face_number			=	0 ;
	global_cell_number			=	0 ;
	local_node_number			=	0 ;
	local_face_number			=	0 ;
	local_cell_number			=	0 ;
	ghost_node_number			=	0 ;
	ghost_node_number_level_1	=	0 ;
	ghost_face_number			=	0 ;
	ghost_face_number_level_1	=	0 ;
	ghost_cell_number			=	0 ;
	ghost_cell_number_level_1	=	0 ;
	ghost_cell_number_level_2	=	0 ;
}

Domain::~Domain( )
{
	ID_bank.erase( ID ) ;
}

Domain::Domain( MPI_Comm mc )
{
	parent_comm	=	MPI_COMM_WORLD ;
	comm		=	mc ;
	MPI_Comm_size ( comm, &comm_size ) ;
	MPI_Comm_rank ( comm, &comm_rank ) ;
	MPI_Comm_size ( parent_comm, &world_comm_size ) ;
	MPI_Comm_rank ( parent_comm, &world_comm_rank ) ;
	ID		=	generate_id() ;
	ID_bank.insert( ID ) ;

	flag_geometry				=	false ;
	flag_boundaryfile			=	false ;
	flag_meshfile				=	false ;
	flag_mesh_scale				=	false ;
	flag_structured_mesh		=	false ;

	dimension 					=	0 ;
	cylindrical_x				=	0 ;
	cylindrical_y				=	0 ;
	global_node_number			=	0 ;
	global_face_number			=	0 ;
	global_cell_number			=	0 ;
	local_node_number			=	0 ;
	local_face_number			=	0 ;
	local_cell_number			=	0 ;
	ghost_node_number			=	0 ;
	ghost_node_number_level_1	=	0 ;
	ghost_face_number			=	0 ;
	ghost_face_number_level_1	=	0 ;
	ghost_cell_number			=	0 ;
	ghost_cell_number_level_1	=	0 ;
	ghost_cell_number_level_2	=	0 ;
}

Domain::Domain( string input_filename )
{
	parent_comm	=	MPI_COMM_WORLD ;
	comm	=	MPI_COMM_WORLD ;
	MPI_Comm_size ( comm, &comm_size ) ;
	MPI_Comm_rank ( comm, &comm_rank ) ;
	world_comm_size	=	comm_size ;
	world_comm_rank	=	comm_rank ;
	ID		=	generate_id() ;
	ID_bank.insert( ID ) ;

	flag_geometry				=	false ;
	flag_boundaryfile			=	false ;
	flag_meshfile				=	false ;
	flag_mesh_scale				=	false ;
	flag_structured_mesh		=	false ;

	dimension 					=	0 ;
	cylindrical_x				=	0 ;
	cylindrical_y				=	0 ;
	global_node_number			=	0 ;
	global_face_number			=	0 ;
	global_cell_number			=	0 ;
	local_node_number			=	0 ;
	local_face_number			=	0 ;
	local_cell_number			=	0 ;
	ghost_node_number			=	0 ;
	ghost_node_number_level_1	=	0 ;
	ghost_face_number			=	0 ;
	ghost_face_number_level_1	=	0 ;
	ghost_cell_number			=	0 ;
	ghost_cell_number_level_1	=	0 ;
	ghost_cell_number_level_2	=	0 ;

	Init( input_filename ) ;
}

Domain::Domain( string input_filename, int *w )
{
	parent_comm	=	MPI_COMM_WORLD ;
	comm		=	MPI_COMM_WORLD ;
	MPI_Comm_size ( comm, &comm_size ) ;
	MPI_Comm_rank ( comm, &comm_rank ) ;
	world_comm_size	=	comm_size ;
	world_comm_rank	=	comm_rank ;
	ID		=	generate_id() ;
	ID_bank.insert( ID ) ;

	flag_geometry				=	false ;
	flag_boundaryfile			=	false ;
	flag_meshfile				=	false ;
	flag_mesh_scale				=	false ;
	flag_structured_mesh		=	false ;

	dimension					=	0 ;
	cylindrical_x				=	0 ;
	cylindrical_y				=	0 ;
	global_cell_number			=	0 ;
	global_face_number			=	0 ;
	global_node_number			=	0 ;
	local_node_number			=	0 ;
	local_face_number			=	0 ;
	local_cell_number			=	0 ;
	ghost_cell_number			=	0 ;
	ghost_cell_number_level_1	=	0 ;
	ghost_cell_number_level_2	=	0 ;
	ghost_node_number			=	0 ;

	Init( input_filename, w ) ;
}

Domain::Domain( string _input_filename, int size )
{
	MPI_Group	world_group, mpi_group ;
	char		buf[ 300 ] ;
	int			i ;

	parent_comm	=	MPI_COMM_WORLD ;
	MPI_Barrier ( parent_comm )  ;

	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;

	MPI_Comm_size ( parent_comm, &world_comm_size ) ;
	MPI_Comm_rank ( parent_comm, &world_comm_rank ) ;

	flag_geometry			=	false ;
	flag_boundaryfile		=	false ;
	flag_meshfile			=	false ;
	flag_mesh_scale			=	false ;
	flag_structured_mesh	=	false ;

	if ( size > world_comm_size )
	{
		size	=	world_comm_size ;
	}
	group_ranks	=	boost::shared_array<int> ( new int [ size ] ) ;
	for ( i = 0 ; i < size ; i++ )
		group_ranks[ i ]	=	i ;

	MPI_Comm_group( parent_comm , &world_group ) ;
	MPI_Group_incl( world_group, size, group_ranks.get(), &mpi_group ) ;
	MPI_Comm_create_group( parent_comm, mpi_group, 0, &comm ) ;

	//comm = parent_comm ;

	comm_size	=	-1 ;
	comm_rank	=	-1 ;
	if ( comm != MPI_COMM_NULL )
	{
		MPI_Comm_size ( comm, &comm_size ) ;
		MPI_Comm_rank ( comm, &comm_rank ) ;
	}
	MPI_Group_free( &world_group ) ;
	MPI_Group_free( &mpi_group ) ;

	dimension					=	0 ;
	cylindrical_x				=	0 ;
	cylindrical_y				=	0 ;
	global_cell_number			=	0 ;
	global_face_number			=	0 ;
	global_node_number			=	0 ;
	local_node_number			=	0 ;
	local_face_number			=	0 ;
	local_cell_number			=	0 ;
	ghost_cell_number			=	0 ;
	ghost_cell_number_level_1	=	0 ;
	ghost_cell_number_level_2	=	0 ;
	ghost_node_number			=	0 ;

	Init( _input_filename ) ;

	// Only the local communicator know the meshfile name at this moment. We need to broadcast the meshfile name to all the processor, in which the meshfile is needed in the DataManager and DataExchanger class.
	MPI_Barrier( parent_comm )  ;
	if ( world_comm_rank == 0 )
	{
		i	=	meshfile.size() + 1 ;  // +1 to keep the '\0' in the end of the string
		strcpy( buf, meshfile.c_str() ) ;
	}
	MPI_Bcast( &i , 1, MPI_INT, 0, parent_comm ) ;
	MPI_Bcast( buf, i, MPI_CHAR, 0, parent_comm ) ;

	string s( buf ) ;
	meshfile	=	s ;
}

Domain::Domain( string _input_filename, int size, int *w )
{
	MPI_Group	world_group, mpi_group ;
	char		buf[ 300 ] ;
	int			i ;

	parent_comm	=	MPI_COMM_WORLD ;
	MPI_Barrier( parent_comm )  ;

	ID	=	generate_id() ;
	ID_bank.insert( ID ) ;

	MPI_Comm_size( parent_comm, &world_comm_size ) ;
	MPI_Comm_rank( parent_comm, &world_comm_rank ) ;

	flag_geometry		=	false ;
	flag_boundaryfile	=	false ;
	flag_meshfile		=	false ;
	flag_mesh_scale		=	false ;
	flag_structured_mesh	=	false ;

	if ( size > world_comm_size )
	{
		size	=	world_comm_size ;
	}
	group_ranks	=	boost::shared_array<int> ( new int [ size ] ) ;
	for ( i = 0 ; i < size ; i++ )
		group_ranks[ i ]	=	i ;

	MPI_Comm_group( parent_comm, &world_group ) ;
	MPI_Group_incl( world_group, size, group_ranks.get(), &mpi_group ) ;
	MPI_Comm_create_group( parent_comm, mpi_group, 0, &comm ) ;

	//comm = parent_comm;

	comm_size	=	-1 ;
	comm_rank	=	-1 ;
	if (  comm != MPI_COMM_NULL )
	{
		MPI_Comm_size( comm, &comm_size ) ;
		MPI_Comm_rank( comm, &comm_rank ) ;
	}
	MPI_Group_free( &world_group ) ;
	MPI_Group_free( &mpi_group ) ;

	dimension					=	0 ;
	cylindrical_x				=	0 ;
	cylindrical_y				=	0 ;
	global_cell_number			=	0 ;
	global_face_number			=	0 ;
	global_node_number			=	0 ;
	local_node_number			=	0 ;
	local_face_number			=	0 ;
	local_cell_number			=	0 ;
	ghost_cell_number			=	0 ;
	ghost_cell_number_level_1	=	0 ;
	ghost_cell_number_level_2	=	0 ;
	ghost_node_number			=	0 ;

	Init( _input_filename , w ) ;

	// Only the local communicator know the meshfile name at this moment. We need to broadcast the meshfile name to all the processor, in which the meshfile is needed in the DataManager and DataExchanger class.
	MPI_Barrier ( parent_comm )  ;
	if ( world_comm_rank == 0 )
	{
		i	=	meshfile.size() + 1 ;  // +1 to keep the '\0' in the end of the string
		strcpy( buf, meshfile.c_str() ) ;
	}
	MPI_Bcast ( &i , 1, MPI_INT, 0, parent_comm ) ;
	MPI_Bcast ( buf, i, MPI_CHAR, 0, parent_comm ) ;

	string	s( buf ) ;
	meshfile	=	s ;
}

void Domain::set_geometry ( string Geometry )
{
	flag_geometry	=	true ;
	geometry 		=	Geometry ;
	
	dimension 		=	3 ;
	if ( geometry != "3D" )
	{
		dimension 	=	2 ;

		if ( geometry == "2D" )
		{

		} else if ( geometry == "Axisymmetric_X" )
		{
			cylindrical_x	=	1 ;
		} else if ( geometry == "Axisymmetric_Y" )
		{
			cylindrical_y	=	1 ;
		}
	}	
}
void Domain::set_structured_mesh ( int nx, int ny, double dx, double dy, vector<string> *BC_Typename )
{
	set_structured_mesh( nx, ny, 1, dx, dy, dx, BC_Typename ) ;
}

void Domain::set_structured_mesh ( int nx, int ny, int nz, double dx, double dy, double dz, vector<string> *BC_Typename )
{
	Nx =	nx ;
	Ny =	ny ;
	Nz =	nz ;

	Dx =	dx ;
	Dy =	dy ;
	Dz =	dz ;

	for ( int i = 0 ; i < BC_Typename->size() ; i++ )
		Mesh.BC_Typename.push_back( (*BC_Typename)[ i ] ) ;

	flag_structured_mesh	=	true ;
}

void Domain::set_meshfile ( string Meshfile )
{
	meshfile	=	Meshfile ;
	flag_meshfile	=	true ;
}

void Domain::set_scale ( double scale )
{
	mesh_scale		=	scale ;
	flag_mesh_scale	=	true ;
}

void Domain::set_processor_number ( int _n_cpu )
{
	// force the comm_size to the required cpu number
	MPI_Group	world_group, mpi_group ;
	char		buf[ 300 ] ;
	int		i ;
	int		size ;

	size		=	_n_cpu ;
	parent_comm	=	MPI_COMM_WORLD ;

	MPI_Barrier( parent_comm )  ;
	MPI_Comm_size( parent_comm, &world_comm_size ) ;
	MPI_Comm_rank( parent_comm, &world_comm_rank ) ;

	if ( size > world_comm_size )
	{
		cout << "The size of using processor is larger than that system can provide. Using the default number instead = " << world_comm_size << endl;
		size	=	world_comm_size ;
	}

	group_ranks	=	boost::shared_array<int> ( new int [ size ] ) ;
	for ( i = 0 ; i < size ; i++ )
		group_ranks[ i ]	=	i ;

	MPI_Comm_group( parent_comm, &world_group ) ;
	MPI_Group_incl( world_group, size, group_ranks.get(), &mpi_group ) ;
	MPI_Comm_create_group( parent_comm, mpi_group, 0, &comm ) ;

	//comm = parent_comm ;

	comm_size	=	-1 ;
	comm_rank	=	-1 ;
	if (  comm != MPI_COMM_NULL )
	{
		MPI_Comm_size ( comm, &comm_size ) ;
		MPI_Comm_rank ( comm, &comm_rank ) ;
	}
	MPI_Group_free( &world_group ) ;
	MPI_Group_free( &mpi_group ) ;

	MPI_Barrier ( parent_comm )  ;
}

void Domain::set_mpi_comm ( MPI_Comm _comm )
{
	parent_comm	=	MPI_COMM_WORLD ;
	comm		=	_comm ;

	MPI_Comm_size( parent_comm, &world_comm_size ) ;
	MPI_Comm_rank( parent_comm, &world_comm_rank ) ;

	comm_size	=	-1 ;
	comm_rank	=	-1 ;
	if (  comm != MPI_COMM_NULL )
	{
		MPI_Comm_size( comm, &comm_size ) ;
		MPI_Comm_rank( comm, &comm_rank ) ;
	}
}

//****************************************************************************************************************************************************
// New Domain. by MH
//****************************************************************************************************************************************************
void Domain::Init( string s_input_filename )
{
	input_filename	=	s_input_filename ;

	if (  comm != MPI_COMM_NULL )
	{
		Log().TagDump( logLEVEL4 ) << "Preprocessing..." ;
		Preprocessing ( ) ;

		Log().TagDump( logLEVEL4 ) << "Indexing..." ;
		Element_indexing() ;

		Log().TagDump( logLEVEL4 ) << "Calculating coefficients..." ;
		CalculateFacedAdL( ) ;
		CalculateFaceSign( ) ;
		CalculateFaceCoefficient( ) ;
		CalculateCellFaceAreaComponent( ) ;
		FaceTagging() ;
	}
}

void Domain::Init( string s_input_filename , int *weight )
{
	input_filename	=	s_input_filename ;

	if (  comm != MPI_COMM_NULL )
	{
		Log().TagDump( logLEVEL4 ) << "Preprocessing..." ;
		Preprocessing ( weight ) ;

		Log().TagDump( logLEVEL4 ) << "Indexing..." ;
		Element_indexing() ;

		Log().TagDump( logLEVEL4 ) << "Calculating coefficients..." ;
		CalculateFacedAdL( ) ;
		CalculateFaceSign( ) ;
		CalculateFaceCoefficient( ) ;
		CalculateCellFaceAreaComponent( ) ;
		FaceTagging() ;
	}
}

void Domain::Init ( )
{
	int	*NULLINT ;

	NULLINT	=	NULL ;
	Init ( NULLINT ) ;
}

void Domain::Init( int *weight )
{
	int	i , j, k ;

	if (  comm != MPI_COMM_NULL )
	{
		if ( !flag_structured_mesh && !flag_meshfile )
		{
			Log().TagDump( logERROR ) << "Need mesh file or define structured grid." ;
			exit( 1 ) ;
		}

		if ( !flag_geometry )
		{
			Log().TagDump( logERROR ) << "Please define geometry, which should be 2D, 3D, Axisymmetric_Y, or Axisymmetric_X" ;

			exit( 1 ) ;
		}

		if ( !flag_mesh_scale )
		{
			Log().TagDump( logLEVEL4 ) << "Scale undefined in the domain. Use 1 as default." ;
			mesh_scale = 1. ;
		}

		Log().TagDump( logLEVEL4 ) << "Preprocessing..." ;
		Preprocessing ( weight ) ;

		Log().TagDump( logLEVEL4 ) << "Indexing..." ;
		Element_indexing() ;

		Log().TagDump( logLEVEL4 ) << "Calculating coefficients..." ;
		CalculateFacedAdL( ) ;
		CalculateFaceSign( ) ;
		CalculateFaceCoefficient( ) ;
		CalculateCellFaceAreaComponent( ) ;
		FaceTagging() ;
	}
}

void Domain::set_table_processor( int size )
{
	Processor_MeshCell 				=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
	Processor_GhostMeshCell_lv1 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
	Processor_GhostMeshCell_lv2 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
	Processor_GhostMeshCell_lv3 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;

	Processor_MeshFace 				=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
	Processor_GhostMeshFace_lv1 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;

	Processor_MeshNode 				=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
	Processor_GhostMeshNode_lv1 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
}

void Domain::set_table_globalcell( int size )
{
	GlobalCell_MeshCellNo	=	boost::shared_array<int> ( new int [ size ] ) ;
	GlobalCell_LocalCellNo	=	boost::shared_array<int> ( new int [ size ] ) ;

	MeshCell_GlobalCellNo	=	boost::shared_array<int> ( new int [ size ] ) ;
	MeshCell_LocalCellNo	=	boost::shared_array<int> ( new int [ size ] ) ;

	for ( int i = 0 ; i < size ; i++ )
	{
		GlobalCell_MeshCellNo[ i ]	=	-999 ;
		GlobalCell_LocalCellNo[ i ]	=	-999 ;

		MeshCell_GlobalCellNo[ i ]	=	-999 ;
		MeshCell_LocalCellNo[ i ]	=	-999 ;
	}
}

void Domain::set_table_localcell( int size )
{
	LocalCell_GlobalCellNo	=	boost::shared_array<int> ( new int [ size ] ) ;
	LocalCell_MeshCellNo	=	boost::shared_array<int> ( new int [ size ] ) ;

	for ( int i = 0 ; i < size ; i++ )
	{
		LocalCell_GlobalCellNo[ i ]	=	-999 ;
		LocalCell_MeshCellNo[ i ]	=	-999 ;
	}
}
void Domain::set_table_face( int m_size, int l_size )
{
	MeshFace_LocalFaceNo	=	boost::shared_array<int> ( new int [ m_size ] ) ;
	LocalFace_MeshFaceNo	=	boost::shared_array<int> ( new int [ l_size ] ) ;

	for ( int i = 0 ; i < m_size ; i++ )
		MeshFace_LocalFaceNo[ i ]	=	-999 ;

	for ( int i = 0 ; i < l_size ; i++ )
		LocalFace_MeshFaceNo[ i ]	=	-999 ;
}

void Domain::set_table_node( int m_size, int l_size )
{
	MeshNode_LocalNodeNo	=	boost::shared_array<int> ( new int [ m_size ] ) ;
	LocalNode_MeshNodeNo	=	boost::shared_array<int> ( new int [ l_size ] ) ;

	for ( int i = 0 ; i < m_size ; i++ )
		MeshNode_LocalNodeNo[ i ]	=	-999 ;

	for ( int i = 0 ; i < l_size ; i++ )
		LocalNode_MeshNodeNo[ i ]	=	-999 ;
}

//****************************************************************************************************************************************************

/*! \brief Calculation of the face coefficient.

For 2D cylindrical coodinate system, we assume the problem is axisymmetric and the \f$\phi\f$-direction is omitted. Here we use the finite volumn method for diffusion equation in cylindrical coodinate system as an example:
\f[
	\int_V \frac{\partial n }{ \partial t } dV + \int_V \nabla \cdot \left( D \nabla n \right) dV = 0
\f]
\f[
	\int_V \frac{\partial n }{ \partial t } dV + \int_A \left( D \vec{\nabla} n \right) d\vec{A} = 0
\f]
\f[
	\int_V \frac{\partial n }{ \partial t } r dr d\phi dz + \int_A \left( D \vec{\nabla} n \right) ( r d\phi dz \vec{e}_{r} + dr dz \vec{e}_{\phi} + r dr d\phi \vec{e}_{z}  ) = 0
\f]
The volumn integral term is angle independent and the intetion of \f$\phi\f$ is \f$2\pi\f$. In the second term, we know the axisymmetric problem will not have the surface on \f$\vec{e}_{\phi}\f$. The same, the angle independent leads a constant \f$2\pi\f$.
\f[
	2\pi \int_V \frac{\partial n }{ \partial t } r dr dz + 2\pi \int_A \left( D \vec { \nabla } n \right) ( r dz \vec{e}_{r} + r dr \vec{e}_{z}  ) = 0
\f]
For convenient, the constant \f$2\pi\f$ will be ignored in the following derivation.
\f[
	\int_V \frac{\partial n }{ \partial t } r dr dz + \int_A \left( D \vec{\nabla} n \right) ( r dz \vec{e}_{r} + r dr \vec{e}_{z}  ) = 0
\f]

For more detail about the axisymmetric FVM, refer to the paper G. Yu et al, Procedia Computer Science 18 (2013) pp.2117. For the detail of the unstructured FVM, refer to the draft note "Numberical Methods in Heat, Mass, and Momentum Transfer" by Jayathi Y. Murthy from School of Mechanical Engineering Purdue University.

We have to note the difference between the 2D/3D Cartesian and cylindircal coorinates. Here we use 2D as an example. We can defined a face by given 2 points \f$p_{1}\f$ and \f$p_{2}\f$, also a vector lies on this face \f$\vec{v} = \vec{p_{2}p_{1}}\f$. A face normal vector can be obtained by the cross product \f$ \vec{A} = \vec{v} \times \vec{e}_{z} \f$

\f[
	\vec{v} = { p_{2,x} - p_{1,x} ,  p_{2,y} - p_{1,y} , 0 }
\f]
\f[
	\vec{A} = { Ax, Ay, Az } = \vec{v} \times \vec{e}_{z} =  \left|
	\begin{array}{ccc}
	\vec{e}_{x} & \vec{e}_{y} & \vec{e}_{z} \\
	p_{2,x} - p_{1,x}& p_{2,y} - p_{1,y} & 0 \\
	0 & 0 & 1
	\end{array} \right|
\f]
\f[
	Ax = p_{2,y} - p_{1,y}
\f]
\f[
	Ay = - ( p_{2,x} - p_{1,x} )
\f]

For the 2D cylindrical coordinate in r and z, we can define a face by given 2 points \f$p_{1}\f$ and \f$p_{2}\f$, also a vector lies on this face \f$\vec{v} = \vec{p_{2}p_{1}}\f$, which is similar to the 2D Cartesian case. The face normal vector can be obtained by the cross product \f$ \vec{A} = \vec{v} \times \vec{e}_{\phi} \f$

\f[
	\vec{v} = \bar{r} { p_{2,r} - p_{1,r} ,  0 , p_{2, z} - p_{1, z}  }
\f]
, where \f$\bar{r} = 0.5 \times ( r_{1} + r_{2} ) \f$
\f[
	\vec{A} = { Ar, A\phi, Az } = \vec{v} \times \vec{e}_{\phi} =  \left|
	\begin{array}{ccc}
	\vec{e}_{r} & \vec{e}_{\phi} & \vec{e}_{z} \\
	p_{2,r} - p_{1,r} &  0 &  p_{2,z} - p_{1,z} \\
	0 & 1 & 0
	\end{array} \right|
\f]
\f[
	Ar = - ( p_{2,z} - p_{1,z} )
\f]
\f[
	Az =  p_{2,r} - p_{1,r}
\f]

*/

void Domain::cell_to_node ( double *source_at_cell, double *target_at_node )
{
	int	i, j ;

	for ( i = 0 ; i < local_node_number ; i++ )
	{
		target_at_node[ i ]	=	0. ;
		for ( j = 0 ; j < node[ i ].cell_number ; j++ )
		{
			target_at_node[ i ]	+=	node[ i ].alpha[ j ] * source_at_cell[ node[ i ].cell[ j ]->local_id ] ;
		}

	}
}

void Domain::cell_to_face ( double *source_at_cell, double *target_at_face )
{
	int	i, j ;

	for ( i = 0 ; i < local_face_number ; i++ )
	{
		target_at_face[ i ]	=	0. ;
		for ( j = 0 ; j < face[ i ].cell_number ; j++ )
			target_at_face[ i ]	+=	face[ i ].alpha[ j ] * source_at_cell[ face[ i ].cell[ j ]->local_id ] ;
	}
}

void Domain::node_to_face ( double *source_at_node, double *target_at_face )
{
	int	i, j ;

	for ( i = 0 ; i < local_face_number ; i++ )
		target_at_face[ i ]	=	0.5 * ( source_at_node [ face[ i ].node[ 0 ]->local_id ] + source_at_node [ face[ i ].node[ 1 ]->local_id ] ) ;
}

void Domain::Output_MeshInformation( string Meshfile, double scale )
{
	int			j, k, MeshNo ;
	ofstream	Output ;
	string		filename ;
	int			mpi_size, mpi_rank ;

	mpi_size	=	comm_size ;
	mpi_rank	=	comm_rank ;

	if ( mpi_rank == 0 )
	{
		filename	=	"Mesh_" ;
		filename.append( Meshfile.substr( 0, Meshfile.size() - 4 ) ) ;
		filename.append( ".dat" ) ;

		Output.open( filename.c_str(), ios::out | ios::trunc ) ;

		Output << "CellNum = " << global_cell_number << endl ;

		for ( int j = 0 ; j < global_cell_number ; j++ )
		{
			GlobalCell_MeshCellNo[ j ] ;
			Output << setw( 20 ) << Mesh.Cell_Position[ 0 ][ MeshNo ] << setw( 20 ) << Mesh.Cell_Position[ 1 ][ MeshNo ] << setw( 20 ) << Mesh.Cell_Position[ 2 ][ MeshNo ] << setw( 20 ) << Mesh.Cell_Node[ MeshNo ].size() ;

			for ( k = 0 ; k <  Mesh.Cell_Node[ MeshNo ].size()  ; k++ )
				Output << setw( 20 ) << ( Mesh.Cell_Node[ MeshNo ][ k ] + 1 ) ;
			Output << endl ;
		}

		Output << "NodeNum = " << global_node_number << endl ;

		for ( int j = 0 ; j < global_node_number ; j++ )
		{
			Output << setw( 20 ) << Mesh.Node_Position[ 0 ][ j ] << setw( 20 ) << Mesh.Node_Position[ 1 ][ j ] << setw( 20 ) << Mesh.Node_Position[ 2 ][ j ] << setw( 20 )  << Mesh.Node_Cell[ j ].size() ;
			for ( int k = 0 ; k < Mesh.Node_Cell[ j ].size() ; k++ )
				Output << setw( 20 ) << Mesh.Node_Cell[ j ][ k ] ;
			Output << endl ;
		}
		Output.clear() ;
		Output.close() ;
	}
	MPI_Barrier( comm ) ;
}

Domain& Domain::operator<< ( ostream & (*f)( ostream& ) )
{

	int	i, field_number;
	int	mpi_size, mpi_rank ;
	boost::shared_array<string> output_field_name ;
	boost::shared_array<double *> output_field_data ;
	string	output_filename ;

	mpi_size	=	comm_size ;
	mpi_rank	=	comm_rank ;

	if ( mpi_rank == 0 && ( ( output_fields.size() - 1 ) !=  output_data.size() ) )
		cout << "Output field captions and data are not matched!" << endl ;

	// For Tecplot format
	field_number		=	output_fields.size() - 1 ;
	output_field_name	=	boost::shared_array<string> ( new string [ field_number ] ) ;
	for ( i = 0 ; i < field_number ; i++)
	{
		output_field_name[ i ]	=	output_fields[ i + 1 ] ;
	} ;
	output_field_data	=	boost::shared_array<double *>  ( new double * [ field_number ] ) ;
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_data [ i ]	=	output_data[ i ] ;
	} ;
	output_filename	=	output_fields[ 0 ] + ".dat" ;

	if ( dimension == 2 )
	{
		//tecplot_output_2D ( field_number, output_field_name.get() , output_field_data.get() , output_filename , this ) ;
	} else if ( dimension == 3 )
	{
		//tecplot_output_3D ( field_number, output_field_name.get() , output_field_data.get() , output_filename , this ) ;
	} ;
	// End of Tecplot format

	output_fields.clear() ;
	output_data.clear() ;

	return *this ;
}

Domain & Domain::operator<< ( double * v )
{
	output_data.push_back( v ) ;
	return *this ;
}

Domain & Domain::operator<< ( const string& c )
{
	output_fields.push_back( c ) ;
	return *this ;
}


void Domain::tecplot_dump_zone( )
{
	string	s( "empty" ) ;
	tecplot_dump_zone( s, 0 ) ;
}

void Domain::tecplot_dump_zone( string zone_info, int NewFileTag )
{
	int	i, field_number ;
	int	mpi_size, mpi_rank ;
	string	*output_field_name ;
	double	**output_field_data ;
	string	output_filename ;

	mpi_rank	=	comm_rank ;
	mpi_size	=	comm_size ;

	if ( mpi_rank == 0 && ( ( output_fields.size() - 1 ) !=  output_data.size() ) )
		cout << "Output field captions and data are not matched!" << endl ;

	// For Tecplot format
	field_number		=	output_fields.size() - 1 ;
	output_field_name	=	new string [ field_number ] ;
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_name[ i ]	=	output_fields[ i + 1 ] ;
	} ;

	output_field_data	=	new double * [ field_number ] ;
	for( i = 0 ; i < field_number ; i++)
	{
		output_field_data[ i ]	=	output_data[ i ] ;
	} ;

	output_filename	=	output_fields[ 0 ] + ".dat" ;

	if ( dimension == 2 )
	{
		if( NewFileTag > 0 )
		{
			//tecplot_output_zone_2D( zone_info, field_number, output_field_name, output_field_data, output_filename, this ) ;
		} else
		{
			//tecplot_output_2D( zone_info, field_number, output_field_name, output_field_data, output_filename, this ) ;
		} ;
	} else if ( dimension == 3 )
	{
		if ( NewFileTag > 0 )
		{
			//tecplot_output_zone_3D( zone_info, field_number, output_field_name, output_field_data, output_filename, this ) ;
		} else
		{
			///tecplot_output_3D( zone_info, field_number, output_field_name, output_field_data, output_filename, this ) ;
		} ;
	} ;

	delete [] output_field_name ;
	delete [] output_field_data ;

	return ;
} ;

void Domain::tecplot_global_cell_out(string filename, double *data )
{
	int	dataNum = 1 ;
	string	*dataname ;
	double	**data_data ;

	data_data	=	new double*[ dataNum ] ;
	dataname	=	new string[ dataNum ] ;

	dataname[ 0 ]	=	"dataName1" ;
	data_data[ 0 ]	=	data ;

	if ( dimension == 2 )
	{
		//tecplot_output_2D_global( "global", dataNum, dataname, data_data, filename, this ) ;
	} else if ( dimension == 3 )
	{
		//tecplot_output_3D_global( "global", dataNum, dataname, data_data, filename, this ) ;
	} ;

	delete [] dataname ;
	delete [] data_data ;

	return ;
} ;

void Domain::tecplot_local_cell_out( string filename, double *data )
{
	int	dataNum = 1 ;
	string	*dataname ;
	double	**data_data ;

	data_data	=	new double*[ dataNum ] ;
	dataname	=	new string[ dataNum ] ;

	dataname[ 0 ]	=	"dataName1" ;
	data_data[ 0 ]	=	data ;

	if ( dimension == 2 )
	{
		//tecplot_output_2D( "local", dataNum, dataname, data_data, filename, this ) ;
	} else if ( dimension == 3 )
	{
		//tecplot_output_3D( "local", dataNum, dataname, data_data, filename, this ) ;
	} ;

	delete [] dataname ;
	delete [] data_data ;

	return ;
} ;

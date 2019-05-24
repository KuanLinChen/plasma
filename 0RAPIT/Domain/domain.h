#include <fstream>
#include <set>
#include <mpi.h>
#include "json.hpp"
#include "preprocessing.h"
#include "element.h"

using namespace std;
using nlohmann::json;

#if !defined(__DOMAIN_H)
#define __DOMAIN_H

//#define FormatFIELDVIEW 0
//#define FormatCGNS 1
//#define FormatMSH 2

/*! \brief class Domain

\image latex ArrangementOfCells.eps "The arrangement of the non-orthogonal mesh." width=10cm

*/

class Domain
{
public:


	Domain( ) ;
	Domain( MPI_Comm ) ;
	// create a new communicator is the
	Domain( MPI_Comm, string mesh_filename, int new_size ) ;
	Domain( MPI_Comm, string mesh_filename, int new_size, int *weighting ) ;
	~Domain( ) ;
	Domain( string ) ;
	Domain( string , int *) ; // mesh_file, weighting
	Domain( string , int ) ; // mesh_file, cpu number
	Domain( string , int , int * ) ; // mesh_file, cpu number, weighing

	MPI_Comm	parent_comm ;
	MPI_Comm	comm ;

	//MPI_Group	mpi_group ;
	int 	comm_size, comm_rank ;
	int 	world_comm_size, world_comm_rank ;
	boost::shared_array<int> group_ranks ;

	string ID ;
	static set<string> ID_bank ;

	string	input_filename ;
	string	geometry, boundaryfile, meshfile, restart_meshfile ;
	bool    flag_geometry, flag_boundaryfile, flag_meshfile, flag_mesh_scale ;
	double	min_cell_length ;

	int	dimension ;
	int	cylindrical_x, cylindrical_y ;
	int	global_cell_number ;
	int	global_face_number ;
	int	global_node_number ;

	int	local_node_number ;
	int	local_face_number ;
	int	local_cell_number ;

	int	ghost_cell_number ;
	int	ghost_cell_number_level_1, ghost_cell_number_level_2, ghost_cell_number_level_3 ;
	int	ghost_node_number, ghost_node_number_level_1 ; // This is for FEM, but useless at this moment
	int	ghost_face_number, ghost_face_number_level_1, ghost_face_number_level_2 ;

	boost::shared_array<int>  CellID2MeshOrdering ;
	boost::shared_array<int>  MeshOrdering2CellID ;

	int BCFace_type_number ;
	boost::shared_array< vector<int> > BCFace ;
	vector<int> BCFace_type ;
	vector<string> BCFace_typename ;
	vector<int> InterFace, BulkFace, NormalFace, MPIFace ;

	int BCNum, BCNum_forFace, BCNum_forCell ;

	double mesh_scale ;

	// Collecting all the boundary faces in to some classes. We use 2 dirchlet BCs and 2 neumann BCs in a mesh as an example.
	// classified_face_set_number: How many different boundaries in the domain? classified_face_set_number = 4
	// classified_face_set: if first dirchlet BC is defined as -11 in the mesh file, classified_face_set[0] = 11, which here we only take positive value
	// classified_face_local_number: the number of face of sets. ex. If a boundary is defined as -11 in the meshfile, classified_face_set[0] = 11, then classified_face_local_number [ 11 ]  is the number of -11 faces.
	// classified_face: if there are 100 faces of type -11, the classified_face[11][0] ~ classified_face[11][99] are the pointer to the BC -11 faces.
	//Face	**classified_face [ 200 ] ;
	//int	classified_face_local_number [ 200 ] ;
	//int 	number_of_BC_face ;
	//int 	*classified_face_set ;
	//int 	classified_face_set_number ;

//	void	set_structured_mesh( int nx, int ny, double dx, double dy ) ;
	int 	Nx, Ny, Nz ;
	double 	Dx, Dy, Dz ;
	bool	flag_structured_mesh ;

	void	set_structured_mesh( int nx, int ny, double dx, double dy, vector<string> *BC_Typename ) ;
	void	set_structured_mesh( int nx, int ny, int nz, double dx, double dy, double dz, vector<string> *BC_Typename ) ;

	void 	set_meshfile ( string Meshfile ) ;
	void 	set_geometry ( string Geometry ) ;
	void 	set_scale ( double scale ) ;
	void 	set_processor_number ( int ) ;	// if force the cpu number, the system will generate a mpi communicator, which the use of set_mpi_comm() is not necessary.
	void	set_mpi_comm ( MPI_Comm ) ;		// If force to use the communicator, the use of set_processor_number() is not necessary.

	//*************************************************************************************************************************
	// New version with Mesh library
	//*************************************************************************************************************************
	boost::shared_array< Node > node ;
	boost::shared_array< Face > face ;
	boost::shared_array< Cell > cell ;

	void 	Init ( string ) ;
	void	Init ( string, int * ) ;
	void 	Init (  ) ;
	void 	Init ( int * ) ;

	void	Preprocessing( ) ;
	void	Preprocessing( int * ) ;
	void 	Element_indexing () ;
	void	CalculateFacedAdL( ) ;
	void	CalculateFaceSign() ;
	void	CalculateFaceCoefficient( ) ;
	void	CalculateCellFaceAreaComponent( ) ;
	void	FaceTagging() ;

	MeshData	Mesh ;
	void 	GetMeshInformation( ) ;
	void	DomainPartition( int * ) ;
	void	CreateCellInformation( ) ;
	void	CreateFaceInformation( ) ;

	//-------------------------------------------------------------------------------------------------------------
	// Processor
	//-------------------------------------------------------------------------------------------------------------
	// Processor_MeshNode:			the relation between a processor and its corresponding MeshNodes
	// Processor_GhostMeshFace_lv1:
	// Processor_MeshFace:			the relation between a processor and its corresponding MeshFaces
	// Processor_GhostMeshFace_lv1:
	// Processor_MeshCell:			the relation between a processor and its corresponding MeshCells
	// Processor_GhostMeshCell_lv1:
	// Processor_GhostMeshCell_lv2:
	boost::shared_array< vector< int > >	Processor_MeshNode ;
	boost::shared_array< vector< int > >	Processor_GhostMeshNode_lv1 ;

	boost::shared_array< vector< int > >	Processor_MeshFace ;
	boost::shared_array< vector< int > >	Processor_GhostMeshFace_lv1 ;

	boost::shared_array< vector< int > >	Processor_MeshCell ;
	boost::shared_array< vector< int > >	Processor_GhostMeshCell_lv1 ;
	boost::shared_array< vector< int > >	Processor_GhostMeshCell_lv2 ;
	boost::shared_array< vector< int > >	Processor_GhostMeshCell_lv3 ;

	//-------------------------------------------------------------------------------------------------------------
	// Table: the relation between MeshData, GlobalData & LocalData
	//-------------------------------------------------------------------------------------------------------------
	boost::shared_array<int>  MeshNode_LocalNodeNo ;
	boost::shared_array<int>  LocalNode_MeshNodeNo ;

	boost::shared_array<int>  MeshFace_LocalFaceNo ;
	boost::shared_array<int>  LocalFace_MeshFaceNo ;

	boost::shared_array<int>  GlobalCell_MeshCellNo ;
	boost::shared_array<int>  MeshCell_GlobalCellNo ;

	boost::shared_array<int>  LocalCell_MeshCellNo ;
	boost::shared_array<int>  MeshCell_LocalCellNo ;

	boost::shared_array<int>  GlobalCell_LocalCellNo ;
	boost::shared_array<int>  LocalCell_GlobalCellNo ;

	void	set_table_processor( int size ) ;
	void	set_table_globalcell( int size ) ;
	void	set_table_localcell( int size ) ;
	void	set_table_face( int g_size, int l_size ) ;
	void	set_table_node( int g_size, int l_size ) ;
	//*************************************************************************************************************************

	//void	OutputResultTec_node( double *variable, string output_filename ) ;
	//void	OutputResultTec( double *variable, string output_filename ) ;
	//void	OutputResultTec( string output_filename ) ;
	void	cell_to_node ( double *, double * ) ;
	void	cell_to_face ( double *, double * ) ;
	void	node_to_face ( double *, double * ) ;

	void	Cell_gradient( double *, double ** ) ; // ( du, vec_u )
	void	Cell_gradient( double *, double *, double *, double * ) ; // (du, ux, uy, uz )
	void	Cell_gradient_value( double *, double *, double **, double *, double * ) ;
	void	Cell_gradient_value( double *, double *, double *, double *, double *, double *, double * ) ;

	void	Output_MeshInformation( string filename, double scale ) ;
	void	MappingNodeInCell( Domain *GridToBeMapped ) ;

	// For mesh information
	//double	*node_position[ 3 ] ;
	//int		*cell_nodeNum, **cell_node ;

	boost::shared_array<int>	processor_local_cell_number ;

	// For dumping data, Tecplot version for this moment
	vector<string>   output_fields ;
	vector<double *> output_data ;

	Domain & operator << ( ostream & (*f)( ostream& ) ) ;
	Domain & operator << ( double * ) ;
	Domain & operator << ( const string &  ) ;

	//int iNumOut, out_flag, NumOfShareZone;
	//ofstream *debug_out;

	void tecplot_dump_zone( ) ;
	void tecplot_dump_zone( string zonename, int NewFileTag ) ;
	void tecplot_global_cell_out( string filename, double *data ) ;
	void tecplot_local_cell_out( string filename, double *data ) ;
} ;


//------------------------------------------------------------------
/*
class auto_Hex_grid
{
public:
	bool make_Hex_grid(json *auto_mesh_json, json *grid_json);
};

class HexCellDomain
{
public:

	int iNumNx, iNumNy, iNumNz;
	int iNumCx, iNumCy, iNumCz;
	int iNumCell, iNumNode, iNumFace, iNumBndyFace;
	int bndy_type_X1, bndy_type_X2, bndy_type_Y1, bndy_type_Y2, bndy_type_Z1, bndy_type_Z2;
	int iNumBCs;
	int *BC_TypeID;
	string *BC_Name;

	double dLength[3],dDx[3];

	Cell *HexCell;
	Node *HexNode;
	Face *HexFace;

	ifstream inpFile;
	ofstream *debugfile;

	void initial(string);
	void SetDebugFile(ofstream *);
	void set_bndy(void);
	void write_uns(ofstream *);

protected:
private:

};
*/






#endif

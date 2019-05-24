#include <string>
#include <vector>
#include <set>
#include <boost/shared_array.hpp>

using namespace std ; 

#if !defined(__PREPROCESSING_H)
#define __PREPROCESSING_H

class MeshData ;

class MeshData
{
	public:
		int NodeNum, FaceNum, CellNum, BC_TypeNum, BC_TypeNum_forFace, BCFaceNum, InnerFaceNum, InterFaceNum ;
		//***********************************************************************************************
		// Boundary Condition
		//***********************************************************************************************	
		// BCFace_Type:	
		// BCFace_Typename: 	the typenames of boundary faces	
		// BCFace_Node:		the relation between a boundary face and its vertices
		vector<string>							BC_Typename ;
		boost::shared_array< int > 				BCFace_Type ;
		boost::shared_array< string > 			BCFace_Typename ;
		boost::shared_array< vector< int > >	BCFace_Node ;

		void set_datasize_BCFaceNum( int size )
		{
			BCFace_Type		=	boost::shared_array< int > ( new int [ size ] ) ;
			BCFace_Typename	=	boost::shared_array< string > ( new string [ size ] ) ;
			for ( int i = 0 ; i < size ; i++ )
			{
				BCFace_Type[ i ]		=	0 ;
				BCFace_Typename[ i ]	=	"Bulk" ;
			}

			BCFace_Node 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
		}

		//***********************************************************************************************
		// Node
		//***********************************************************************************************	
		// Node_Position:	the positions of nodes
		// Node_BCFace:		the relation between a node and its corresponding boundary facess
		// Node_Face:		the relation between a node and its corresponding facess
		// Node_Cell:		the relation between a node and its corresponding Cells
		boost::shared_array< double >			Node_Position[ 3 ] ;
		boost::shared_array< vector< int > >	Node_BCFace ;
		boost::shared_array< vector< int > >	Node_Face ;
		boost::shared_array< vector< int > >	Node_Cell ;

		void set_datasize_NodeNum( int size )
		{
			for ( int i = 0 ; i < 3 ; i++ )
			{
				Node_Position[ i ]	=	boost::shared_array< double > ( new double [ size ] ) ;
				for ( int j = 0 ; j < size ; j++ )
					Node_Position[ i ][ j ] =	0. ;
			}

			Node_BCFace 	=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
			Node_Face 		=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
			Node_Cell 		=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
		}

		//***********************************************************************************************
		// Face
		//***********************************************************************************************
		// Face_BCNo:		
		// Face_Node: 	the relation between a face and its vertices
		// Face_Type:		
		// Face_Typename:	
		// Face_Cell:	the relation between a face and its corresponding Cells
		boost::shared_array< double >			Face_Position[ 3 ] ;
		boost::shared_array< int > 				Face_Type ;
		boost::shared_array< string > 			Face_Typename ;
		boost::shared_array< int > 				Face_ProcessorNo ;
		boost::shared_array< vector< int > >	Face_Node ;
		boost::shared_array< vector< int > >	Face_Cell ;

		void set_datasize_FaceNum( int size )
		{
			for ( int i = 0 ; i < 3 ; i++ )
			{
				Face_Position[ i ]	=	boost::shared_array< double > ( new double [ size ] ) ;
				for ( int j = 0 ; j < size ; j++ )
					Face_Position[ i ][ j ]	=	0. ;
			}

			Face_ProcessorNo	=	boost::shared_array< int > ( new int [ size ] ) ;
			Face_Type			=	boost::shared_array< int > ( new int [ size ] ) ;
			Face_Typename 		=	boost::shared_array< string > ( new string [ size ] ) ;
			for ( int i = 0 ; i < size ; i++ )
			{
				Face_Type[ i ]		=	0 ;
				Face_Typename[ i ]	=	"Bulk" ;
			}

			Face_Node		=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
			Face_Cell		=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
		}

		//***********************************************************************************************
		// Cell
		//***********************************************************************************************
		// Cell_Form: 				
		// Cell_BCNo:			
		// Cell_FaceNum:	
		// Cell_ProcessorNo:
		// Cell_Volume:					
		// Cell_Node: 			the relation between a Cell and its vertices
		// Cell_Face: 			the relation between a Cell and its faces
		// Cell_Cell: 			the relation between a Cell and its neighbor Cells
		boost::shared_array< double > 			Cell_Position[ 3 ] ;
		boost::shared_array< int > 				Cell_Form ;
		boost::shared_array< int > 				Cell_Type ;
		boost::shared_array< string > 			Cell_Typename ;		
		boost::shared_array< int > 				Cell_FaceNum ;
		boost::shared_array< int > 				Cell_ProcessorNo ;
		boost::shared_array< double > 			Cell_Volume ;
		boost::shared_array< vector< int > >	Cell_Node ;
		boost::shared_array< vector< int > >	Cell_Face ;
		boost::shared_array< vector< int > >	Cell_Cell ;

		void set_datasize_CellNum( int size )
		{
			for ( int i = 0 ; i < 3 ; i++ )
			{
				Cell_Position[ i ]	=	boost::shared_array< double > ( new double [ size ] ) ;
				for ( int j = 0 ; j < size ; j++ )
					Cell_Position[ i ][ j ]	=	0. ;
			}

			Cell_Form			=	boost::shared_array< int > ( new int [ size ] ) ;
			Cell_Type			=	boost::shared_array< int > ( new int [ size ] ) ;
			Cell_Typename		=	boost::shared_array< string > ( new string [ size ] ) ;
			Cell_FaceNum 		=	boost::shared_array< int > ( new int [ size ] ) ;
			Cell_ProcessorNo 	=	boost::shared_array< int > ( new int [ size ] ) ;
			Cell_Volume		 	=	boost::shared_array< double > ( new double [ size ] ) ;
			for ( int i = 0 ; i < size ; i++ )
			{
				Cell_Form[ i ]			=	-999 ;
				Cell_Type[ i ]			=	0 ;
				Cell_Typename[ i ]		=	"Fluid" ;
				Cell_FaceNum[ i 	]	=	0 ;
				Cell_ProcessorNo[ i ]	=	0 ;
				Cell_Volume[ i ]		=	0. ;
			}

			Cell_Node 			=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
			Cell_Face 			=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
			Cell_Cell 			=	boost::shared_array< vector< int > > ( new vector< int > [ size ] ) ;
		}
} ;

class CellMapping
{
	public:
		int	NodeNum[ 6 ], SurfaceNum[ 6 ], TetraCellNum[ 6 ] ;
		int	SurfaceNodeNum[ 6 ][ 6 ], Node[ 6 ][ 6 ][ 4 ], TetraCell[ 6 ][ 5 ][ 4 ] ;

		// SurfaceNodeNum[ Cell_Form ][ Face_index ]
		// Node[ Cell_Form ][ Face_index ][ Node_index ]
		// TetraCell[ Cell_Form ][ Face_index ][ Node_index ]

		CellMapping()
		{
			for ( int i = 0 ; i < 6 ; i++ )
			{
				for ( int j = 0 ; j < 6 ; j++ )
				{
					SurfaceNodeNum[ i ][ j ]	=	0 ;

					for ( int k = 0 ; k < 4 ; k++ )
						Node[ i ][ j ][ k ]	=	0 ;
				}

				for ( int j = 0 ; j < 5 ; j++ )
				{
					for ( int k = 0 ; k < 4 ; k++ )
						TetraCell[ i ][ j ][ k ]	=	0 ;
				}
			}

			// For Tetra Cell Type
			NodeNum[ 0 ]				=	4 ;
			SurfaceNum[ 0 ]				=	4 ;
			TetraCellNum[ 0 ]			=	1 ;

			SurfaceNodeNum[ 0 ][ 0 ]	=	3 ;
			SurfaceNodeNum[ 0 ][ 1 ]	=	3 ;
			SurfaceNodeNum[ 0 ][ 2 ]	=	3 ;
			SurfaceNodeNum[ 0 ][ 3 ]	=	3 ;

			Node[ 0 ][ 0 ][ 0 ]			=	0 ;
			Node[ 0 ][ 0 ][ 1 ]			=	1 ;
			Node[ 0 ][ 0 ][ 2 ]			=	3 ;

			Node[ 0 ][ 1 ][ 0 ]			=	1 ;
			Node[ 0 ][ 1 ][ 1 ]			=	2 ;
			Node[ 0 ][ 1 ][ 2 ]			=	3 ;

			Node[ 0 ][ 2 ][ 0 ]			=	2 ;
			Node[ 0 ][ 2 ][ 1 ]			=	0 ;
			Node[ 0 ][ 2 ][ 2 ]			=	3 ;

			Node[ 0 ][ 3 ][ 0 ]			=	2 ;
			Node[ 0 ][ 3 ][ 1 ]			=	1 ;
			Node[ 0 ][ 3 ][ 2 ]			=	0 ;

			TetraCell[ 0 ][ 0 ][ 0 ]	=	0 ;
			TetraCell[ 0 ][ 0 ][ 1 ]	=	1 ;
			TetraCell[ 0 ][ 0 ][ 2 ]	=	2 ;
			TetraCell[ 0 ][ 0 ][ 3 ]	=	3 ;
			// End (tetra cell)

			// For Hex. Cell Type	
			NodeNum[ 1 ]				=	8 ;
			SurfaceNum[ 1 ]				=	6 ;
			TetraCellNum[ 1 ]			=	5 ;

			SurfaceNodeNum[ 1 ][ 0 ]	=	4 ;
			SurfaceNodeNum[ 1 ][ 1 ]	=	4 ;
			SurfaceNodeNum[ 1 ][ 2 ]	=	4 ;
			SurfaceNodeNum[ 1 ][ 3 ]	=	4 ;
			SurfaceNodeNum[ 1 ][ 4 ]	=	4 ;
			SurfaceNodeNum[ 1 ][ 5 ]	=	4 ;

			Node[ 1 ][ 0 ][ 0 ]			=	0 ;
			Node[ 1 ][ 0 ][ 1 ]			=	1 ;
			Node[ 1 ][ 0 ][ 2 ]			=	5 ;
			Node[ 1 ][ 0 ][ 3 ]			=	4 ;

			Node[ 1 ][ 1 ][ 0 ]			=	1 ;
			Node[ 1 ][ 1 ][ 1 ]			=	2 ;
			Node[ 1 ][ 1 ][ 2 ]			=	6 ;
			Node[ 1 ][ 1 ][ 3 ]			=	5 ;

			Node[ 1 ][ 2 ][ 0 ]			=	2 ;
			Node[ 1 ][ 2 ][ 1 ]			=	3 ;
			Node[ 1 ][ 2 ][ 2 ]			=	7 ;
			Node[ 1 ][ 2 ][ 3 ]			=	6 ;

			Node[ 1 ][ 3 ][ 0 ]			=	3 ;
			Node[ 1 ][ 3 ][ 1 ]			=	0 ;
			Node[ 1 ][ 3 ][ 2 ]			=	4 ;
			Node[ 1 ][ 3 ][ 3 ]			=	7 ;

			Node[ 1 ][ 4 ][ 0 ]			=	2 ;
			Node[ 1 ][ 4 ][ 1 ]			=	1 ;
			Node[ 1 ][ 4 ][ 2 ]			=	0 ;
			Node[ 1 ][ 4 ][ 3 ]			=	3 ;

			Node[ 1 ][ 5 ][ 0 ]			=	5 ;
			Node[ 1 ][ 5 ][ 1 ]			=	6 ;
			Node[ 1 ][ 5 ][ 2 ]			=	7 ;
			Node[ 1 ][ 5 ][ 3 ]			=	4 ;

			TetraCell[ 1 ][ 0 ][ 0 ]	=	0 ;
			TetraCell[ 1 ][ 0 ][ 1 ]	=	1 ;
			TetraCell[ 1 ][ 0 ][ 2 ]	=	2 ;
			TetraCell[ 1 ][ 0 ][ 3 ]	=	5 ;

			TetraCell[ 1 ][ 1 ][ 0 ]	=	0 ;
			TetraCell[ 1 ][ 1 ][ 1 ]	=	2 ;
			TetraCell[ 1 ][ 1 ][ 2 ]	=	3 ;
			TetraCell[ 1 ][ 1 ][ 3 ]	=	7 ;

			TetraCell[ 1 ][ 2 ][ 0 ]	=	0 ;
			TetraCell[ 1 ][ 2 ][ 1 ]	=	2 ;
			TetraCell[ 1 ][ 2 ][ 2 ]	=	5 ;
			TetraCell[ 1 ][ 2 ][ 3 ]	=	7 ;

			TetraCell[ 1 ][ 3 ][ 0 ]	=	2 ;
			TetraCell[ 1 ][ 3 ][ 1 ]	=	5 ;
			TetraCell[ 1 ][ 3 ][ 2 ]	=	6 ;
			TetraCell[ 1 ][ 3 ][ 3 ]	=	7 ;

			TetraCell[ 1 ][ 4 ][ 0 ]	=	0 ;
			TetraCell[ 1 ][ 4 ][ 1 ]	=	4 ;
			TetraCell[ 1 ][ 4 ][ 2 ]	=	5 ;
			TetraCell[ 1 ][ 4 ][ 3 ]	=	7 ;
			// End (hex. cell)

			// For Pyramid Cell Type
			NodeNum[ 2 ]				=	5 ;
			SurfaceNum[ 2 ]				=	5 ;
			TetraCellNum[ 2 ]			=	2 ;

			SurfaceNodeNum[ 2 ][ 0 ]	=	3 ;
			SurfaceNodeNum[ 2 ][ 1 ]	=	3 ;
			SurfaceNodeNum[ 2 ][ 2 ]	=	3 ;
			SurfaceNodeNum[ 2 ][ 3 ]	=	3 ;
			SurfaceNodeNum[ 2 ][ 4 ]	=	4 ;

			Node[ 2 ][ 0 ][ 0 ]			=	0 ;
			Node[ 2 ][ 0 ][ 1 ]			=	1 ;
			Node[ 2 ][ 0 ][ 2 ]			=	4 ;

			Node[ 2 ][ 1 ][ 0 ]			=	1 ;
			Node[ 2 ][ 1 ][ 1 ]			=	2 ;
			Node[ 2 ][ 1 ][ 2 ]			=	4 ;

			Node[ 2 ][ 2 ][ 0 ]			=	2 ;
			Node[ 2 ][ 2 ][ 1 ]			=	3 ;
			Node[ 2 ][ 2 ][ 2 ]			=	4 ;

			Node[ 2 ][ 3 ][ 0 ]			=	3 ;
			Node[ 2 ][ 3 ][ 1 ]			=	0 ;
			Node[ 2 ][ 3 ][ 2 ]			=	4 ;

			Node[ 2 ][ 4 ][ 0 ]			=	2 ;
			Node[ 2 ][ 4 ][ 1 ]			=	1 ;
			Node[ 2 ][ 4 ][ 2 ]			=	0 ;
			Node[ 2 ][ 4 ][ 3 ]			=	3 ;

			TetraCell[ 2 ][ 0 ][ 0 ]	=	0 ;
			TetraCell[ 2 ][ 0 ][ 1 ]	=	1 ;
			TetraCell[ 2 ][ 0 ][ 2 ]	=	3 ;
			TetraCell[ 2 ][ 0 ][ 3 ]	=	4 ;

			TetraCell[ 2 ][ 1 ][ 0 ]	=	1 ;
			TetraCell[ 2 ][ 1 ][ 1 ]	=	2 ;
			TetraCell[ 2 ][ 1 ][ 2 ]	=	3 ;
			TetraCell[ 2 ][ 1 ][ 3 ]	=	4 ;		
			// End (pyramid cell)

			// For Prism Cell Type
			NodeNum[ 3 ]				=	6 ;
			SurfaceNum[ 3 ]				=	5 ;
			TetraCellNum[ 3 ]			=	3 ;

			SurfaceNodeNum[ 3 ][ 0 ]	=	4 ;
			SurfaceNodeNum[ 3 ][ 1 ]	=	3 ;
			SurfaceNodeNum[ 3 ][ 2 ]	=	4 ;
			SurfaceNodeNum[ 3 ][ 3 ]	=	3 ;
			SurfaceNodeNum[ 3 ][ 4 ]	=	4 ;

			Node[ 3 ][ 0 ][ 0 ]			=	0 ;
			Node[ 3 ][ 0 ][ 1 ]			=	1 ;
			Node[ 3 ][ 0 ][ 2 ]			=	5 ;
			Node[ 3 ][ 0 ][ 3 ]			=	4 ;

			Node[ 3 ][ 1 ][ 0 ]			=	1 ;
			Node[ 3 ][ 1 ][ 1 ]			=	2 ;
			Node[ 3 ][ 1 ][ 2 ]			=	5 ;

			Node[ 3 ][ 2 ][ 0 ]			=	2 ;
			Node[ 3 ][ 2 ][ 1 ]			=	3 ;
			Node[ 3 ][ 2 ][ 2 ]			=	4 ;
			Node[ 3 ][ 2 ][ 3 ]			=	5 ;

			Node[ 3 ][ 3 ][ 0 ]			=	3 ;
			Node[ 3 ][ 3 ][ 1 ]			=	0 ;
			Node[ 3 ][ 3 ][ 2 ]			=	4 ;

			Node[ 3 ][ 4 ][ 0 ]			=	2 ;
			Node[ 3 ][ 4 ][ 1 ]			=	1 ;
			Node[ 3 ][ 4 ][ 2 ]			=	0 ;
			Node[ 3 ][ 4 ][ 3 ]			=	3 ;

			TetraCell[ 3 ][ 0 ][ 0 ]	=	0 ;
			TetraCell[ 3 ][ 0 ][ 1 ]	=	3 ;
			TetraCell[ 3 ][ 0 ][ 2 ]	=	4 ;
			TetraCell[ 3 ][ 0 ][ 3 ]	=	1 ;

			TetraCell[ 3 ][ 1 ][ 0 ]	=	3 ;
			TetraCell[ 3 ][ 1 ][ 1 ]	=	2 ;
			TetraCell[ 3 ][ 1 ][ 2 ]	=	4 ;
			TetraCell[ 3 ][ 1 ][ 3 ]	=	1 ;

			TetraCell[ 3 ][ 2 ][ 0 ]	=	1 ;
			TetraCell[ 3 ][ 2 ][ 1 ]	=	2 ;
			TetraCell[ 3 ][ 2 ][ 2 ]	=	4 ;
			TetraCell[ 3 ][ 2 ][ 3 ]	=	5 ;
			// End (prism cell)

			// For 2-D Triangular Cell Type
			NodeNum[ 4 ]				=	3 ;
			SurfaceNum[ 4 ]				=	3 ;
			TetraCellNum[ 4 ]			=	0 ;

			SurfaceNodeNum[ 4 ][ 0 ]	=	2 ;
			SurfaceNodeNum[ 4 ][ 1 ]	=	2 ;
			SurfaceNodeNum[ 4 ][ 2 ]	=	2 ;

			Node[ 4 ][ 0 ][ 0 ]			=	0 ;
			Node[ 4 ][ 0 ][ 1 ]			=	1 ;

			Node[ 4 ][ 1 ][ 0 ]			=	1 ;
			Node[ 4 ][ 1 ][ 1 ]			=	2 ;

			Node[ 4 ][ 2 ][ 0 ]			=	2 ;
			Node[ 4 ][ 2 ][ 1 ]			=	0 ;
			// End (2-d triangular cell)

			// For 2-D Quadrilateral Cell Type
			NodeNum[ 5 ]				=	4 ;
			SurfaceNum[ 5 ]				=	4 ;
			TetraCellNum[ 5 ]			=	0 ;

			SurfaceNodeNum[ 5 ][ 0 ]	=	2 ;
			SurfaceNodeNum[ 5 ][ 1 ]	=	2 ;
			SurfaceNodeNum[ 5 ][ 2 ]	=	2 ;
			SurfaceNodeNum[ 5 ][ 3 ]	=	2 ;

			Node[ 5 ][ 0 ][ 0 ]			=	0 ;
			Node[ 5 ][ 0 ][ 1 ]			=	1 ;

			Node[ 5 ][ 1 ][ 0 ]			=	1 ;
			Node[ 5 ][ 1 ][ 1 ]			=	2 ;

			Node[ 5 ][ 2 ][ 0 ]			=	2 ;
			Node[ 5 ][ 2 ][ 1 ]			=	3 ;

			Node[ 5 ][ 3 ][ 0 ]			=	3 ;
			Node[ 5 ][ 3 ][ 1 ]			=	0 ;
			// End (2-d quadrilateral cell)
		}
} ;

void Structured_mesh( int dimension, double scale, int nx, int ny, int nz, double dx, double dy, double dz, int NodeNum, int CellNum, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, vector<string> *BC_Typename, int *BCFace_Type, string *BCFace_Typename, int *Cell_FaceNum, vector<int> *Cell_Node, vector<int> *Node_Cell, vector<int> *BCFace_Node, vector<int> *Node_BCFace ) ;

#endif

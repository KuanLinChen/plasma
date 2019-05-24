#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>
#include "preprocessingtool.h"
#include "sys_log.h"
#include "pti.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
using namespace std ;

/*!

@file msh.cpp Functions for reading GMSH mesh.

*/


/*! \brief Read the information from mesh in PTI format

To obtain the numbers from mesh in PTI format, including NodeNum, CellNum, BC_TypeNum, and BCFaceNum.

*/

void GetNumbers_PTI ( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, int *MeshType_Num, string Filename )
{
	int 		i, buffer_int ;
	ifstream	input ;
	string 		buffer ;

	int 		physical_Dim, element_Num ;
	string 		physical_Name, element_Type ;

	input.open ( Filename.c_str ( ) ) ;
	if ( !input )
	{
		cerr << "\n";
		cerr << "PTI_DATA_READ - Fatal error!\n";
		cerr << "  Could not open input file \"" << Filename << "\"\n";
		exit ( 1 ) ;
	}

	// Initial value to the variable
	*NodeNum	=	0 ;
	*CellNum 	=	0 ;
	*BC_TypeNum	=	0 ;
	*BCFaceNum 	=	0 ;

	if ( dimension == 2 )
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$PhysicalName" ) )
			{
				getline( input, buffer ) ;
				(*BC_TypeNum)	=	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < (*BC_TypeNum) ; i++ )
				{
					input >> physical_Dim >> physical_Name ;

					if ( physical_Dim == 1 )
					{
						Log().TagDump( logLEVEL4 )  << "BCFace_Type " << physical_Name ;
					} else if ( physical_Dim == 2 )
					{
						Log().TagDump( logLEVEL4 )  << "Cell_Type " << physical_Name  ;
					}
				}
				//Log().TagDump( logLEVEL4 )  << "BCTypeNum: " << *BC_TypeNum ;
			}

			if ( Compare_string( buffer, "$InfoElement" ) )
			{
				getline( input, buffer ) ;
				buffer_int =	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < buffer_int ; i++ )
				{
					input >> element_Type >> element_Num ;
					if ( element_Type == "Node" )
					{
						(*NodeNum )	=	element_Num ;
					} else if ( element_Type == "Line" )
					{
						(*BCFaceNum)		=	element_Num ;
						MeshType_Num[ 6 ]	=	element_Num ;
					} else if ( element_Type == "Triangular" )
					{
						(*CellNum)			+=	element_Num ;
						MeshType_Num[ 4 ]	=	element_Num ;
					} else if ( element_Type == "Quadrilateral" )
					{
						(*CellNum)			+=	element_Num ;
						MeshType_Num[ 5 ]	=	element_Num ;
					}
				}
			}
		}
	} else
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$PhysicalName" ) )
			{
				getline( input, buffer ) ;
				(*BC_TypeNum)	=	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < (*BC_TypeNum) ; i++ )
				{
					input >> physical_Dim >> physical_Name ;

					if ( physical_Dim == 2 )
					{
						Log().TagDump( logLEVEL4 )  << "BCFace_Type " << physical_Name ;
					} else if ( physical_Dim == 3 )
					{
						Log().TagDump( logLEVEL4 )  << "Cell_Type " << physical_Name  ;
					}
				}
				//Log().TagDump( logLEVEL4 )  << "BCTypeNum: " << *BC_TypeNum ;
			}

			if ( Compare_string( buffer, "$InfoElement" ) )
			{
				getline( input, buffer ) ;
				buffer_int =	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < buffer_int ; i++ )
				{
					input >> element_Type >> element_Num ;

					if ( element_Type == "Node" )
					{
						(*NodeNum )	=	element_Num ;
					} else if ( element_Type == "Triangular" )
					{
						(*BCFaceNum)		+=	element_Num ;
						MeshType_Num[ 4 ]	=	element_Num ;
					} else if ( element_Type == "Quadrilateral" )
					{
						(*BCFaceNum)		+=	element_Num ;
						MeshType_Num[ 5 ]	=	element_Num ;
					} else if ( element_Type == "Tetra")
					{
						(*CellNum)			+=	element_Num ;
						MeshType_Num[ 0 ]	=	element_Num ;
					} else if ( element_Type == "Hexahedron")
					{
						(*CellNum)			+=	element_Num ;
						MeshType_Num[ 1 ]	=	element_Num ;
					} else if ( element_Type == "Prism")
					{
						(*CellNum)			+=	element_Num ;
						MeshType_Num[ 2 ]	=	element_Num ;
					} else if ( element_Type == "Pyramid")
					{
						(*CellNum)			+=	element_Num ;
						MeshType_Num[ 3 ]	=	element_Num ;
					}
				}
			}
		}
	}
	input.clear() ;
	input.close() ;

	Log().TagDump( logLEVEL4 )  << "NodeNum " << (*NodeNum) ;
	Log().TagDump( logLEVEL4 )  << "CellNum " << (*CellNum) ;
	Log().TagDump( logLEVEL4 )  << "BCFaceNum " << (*BCFaceNum)  ;

	(*BC_TypeNum)++ ;
}

void ReadGrid_PTI ( int dimension, double scale, int NodeNum, int *MeshType_Num, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int *Cell_Form, int *Cell_FaceNum, int *Cell_Type, string *Cell_Typename, vector<int> *Cell_Node, vector<int> *Node_Cell, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, vector<int> *Node_BCFace, string Filename )
{
	int 		i, j, No, buffer_int, FaceNo, CellNo ;
	string 		buffer ;
	ifstream	input ;
	int 		physical_Dim, element_PhysGrp, element_NumTag ;
	string 		physical_Name ;
	double 		position[ 3 ] ;

	boost::shared_array< int >		NodeNo ;
	NodeNo =	boost::shared_array< int > ( new int [ 8 ] ) ;
	for ( int i = 0 ; i < 8 ; i++ )
		NodeNo[ i ] = -999 ;

	boost::shared_array< string >	Typename ;
	Typename 		=	boost::shared_array< string > ( new string [ *BC_TypeNum ] ) ;
	Typename[ 0 ]	=	"Bulk" ;

	map<string, int>	check_Typename, Typename_Type ;
	Typename_Type[ "Bulk" ]		=	0 ;
	check_Typename[ "Bulk" ]	=	1 ;

	(*BC_Typename).push_back( "Bulk" ) ;
	(*BC_TypeNum)	=	1 ;

	input.open ( Filename.c_str ( ) ) ;
	while ( getline( input, buffer ) )
	{
		if ( Compare_string( buffer, "$PhysicalName" ) )
		{
			getline( input, buffer ) ;
			buffer_int	=	atoi( buffer.c_str() ) ;

			for ( i = 1 ; i <= buffer_int ; i++ )
			{
				input >> physical_Dim >> physical_Name ;

				Typename[ i ]	=	physical_Name ;

				if ( physical_Dim == ( dimension - 1 ) )
				{
					if ( check_Typename[ physical_Name ] == 0 && physical_Name != "Unspecified" )
					{
						(*BC_Typename).push_back( physical_Name ) ;

						Typename_Type[ physical_Name ]	=	(*BC_TypeNum) ;
						check_Typename[ physical_Name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
				}
			}
		}
	}

	if ( check_Typename[ "Interface" ] == 0 )
	{
		(*BC_Typename).push_back( "Interface" ) ;

		Typename_Type[ "Interface" ]	=	(*BC_TypeNum) ;
		check_Typename[ "Interface" ]	=	1 ;
		(*BC_TypeNum)++ ;
	}
	(*BC_TypeNum_forFace)	=	(*BC_TypeNum) ;	
	input.clear() ;
	input.close() ;

	input.open ( Filename.c_str ( ) ) ;
	while ( getline( input, buffer ) )
	{
		if ( Compare_string( buffer, "$PhysicalName" ) )
		{
			getline( input, buffer ) ;
			buffer_int	=	atoi( buffer.c_str() ) ;

			for ( i = 1 ; i <= buffer_int ; i++ )
			{
				input >> physical_Dim >> physical_Name ;

				Typename[ i ]	=	physical_Name ;

				if ( physical_Dim == dimension )
				{
					if ( check_Typename[ physical_Name ] == 0 && physical_Name != "Unspecified" )
					{
						(*BC_Typename).push_back( physical_Name ) ;

						Typename_Type[ physical_Name ]	=	(*BC_TypeNum) ;
						check_Typename[ physical_Name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
				}
			}
		}
	}
	input.clear() ;
	input.close() ;

	FaceNo 	=	0 ;
	CellNo 	=	0 ;
	input.open ( Filename.c_str ( ) ) ;
	if ( dimension == 3 )
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$Node" ) )
			{
				for ( i = 0 ; i < NodeNum ; i++ )
				{
					input >> position[ 0 ] >> position[ 1 ] >> position[ 2 ] ;

					Node_x[ i ]	=	position[ 0 ] * scale ;
					Node_y[ i ]	=	position[ 1 ] * scale ;
					Node_z[ i ]	=	position[ 2 ] * scale ;
				}
			}

			// BCFace Information
			if ( Compare_string( buffer, "$Triangular" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 4 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;
					
					No =	3 ;
					BCFace_Typename[ FaceNo ]	=	physical_Name ;
					BCFace_Type[ FaceNo ]		=	Typename_Type[ BCFace_Typename[ FaceNo ] ] ;

					for ( j = 0 ; j < No ; j++ )
					{
						input >> NodeNo[ j ] ;

						BCFace_Node[ FaceNo ].push_back( NodeNo[ j ] ) ;
						Node_BCFace[ NodeNo[ j ] ].push_back( FaceNo ) ;
					}
					FaceNo++ ;
				}
			}

			// BCFace Information
			if ( Compare_string( buffer, "$Quadrilateral" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 5 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					No =	4 ;
					BCFace_Typename[ FaceNo ]	=	physical_Name ;
					BCFace_Type[ FaceNo ]		=	Typename_Type[ BCFace_Typename[ FaceNo ] ] ;

					for ( j = 0 ; j < No ; j++ )
					{
						input >> NodeNo[ j ] ;
						BCFace_Node[ FaceNo ].push_back( NodeNo[ j ] ) ;
						Node_BCFace[ NodeNo[ j ] ].push_back( FaceNo ) ;
					}
					FaceNo++ ;
				}
			}

			// Cell Information
			if ( Compare_string( buffer, "$Tetra" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 0 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					Cell_Typename[ CellNo ]	=	physical_Name ;
					Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

					No 						=	4 ;
					Cell_Form[ CellNo ]		=	0 ;
					Cell_FaceNum[ CellNo ]	=	4 ;

					for ( j = 0 ; j < No ; j++ )
						input >> NodeNo[ j ] ;

					Cell_x[ CellNo ]	=	0.0 ;
					Cell_y[ CellNo ]	=	0.0 ;
					Cell_z[ CellNo ]	=	0.0 ;
					for ( j = 0 ; j < No ; j++ )
					{
						Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
						Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
					}
					Cell_x[ CellNo ]	/=	No ;
					Cell_y[ CellNo ]	/=	No ;
					Cell_z[ CellNo ]	/=	No ;

					CellNo++ ;
				}
			}

			if ( Compare_string( buffer, "$Hexahedron" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 1 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					Cell_Typename[ CellNo ]	=	physical_Name ;
					Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

					No 						=	8 ;
					Cell_Form[ CellNo ]		=	1 ;
					Cell_FaceNum[ CellNo ]	=	6 ;

					for ( j = 0 ; j < No ; j++ )
						input >> NodeNo[ j ] ;

					Cell_x[ CellNo ]	=	0.0 ;
					Cell_y[ CellNo ]	=	0.0 ;
					Cell_z[ CellNo ]	=	0.0 ;
					for ( j = 0 ; j < No ; j++ )
					{
						Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
						Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
					}
					Cell_x[ CellNo ]	/=	No ;
					Cell_y[ CellNo ]	/=	No ;
					Cell_z[ CellNo ]	/=	No ;

					CellNo++ ;
				}
			}

			if ( Compare_string( buffer, "$Prism" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 2 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					Cell_Typename[ CellNo ]	=	physical_Name ;
					Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

					No 						=	6 ;
					Cell_Form[ CellNo ]		=	2 ;
					Cell_FaceNum[ CellNo ]	=	5 ;

					for ( j = 0 ; j < No ; j++ )
						input >> NodeNo[ j ] ;

					Cell_x[ CellNo ]	=	0.0 ;
					Cell_y[ CellNo ]	=	0.0 ;
					Cell_z[ CellNo ]	=	0.0 ;
					for ( j = 0 ; j < No ; j++ )
					{
						Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
						Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
					}
					Cell_x[ CellNo ]	/=	No ;
					Cell_y[ CellNo ]	/=	No ;
					Cell_z[ CellNo ]	/=	No ;

					CellNo++ ;
				}
			}

			if ( Compare_string( buffer, "$Pyramid" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 3 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					Cell_Typename[ CellNo ]	=	physical_Name ;
					Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

					No 						=	5 ;
					Cell_Form[ CellNo ]		=	3 ;
					Cell_FaceNum[ CellNo ]	=	5 ;

					for ( j = 0 ; j < No ; j++ )
						input >> NodeNo[ j ] ;

					Cell_x[ CellNo ]	=	0.0 ;
					Cell_y[ CellNo ]	=	0.0 ;
					Cell_z[ CellNo ]	=	0.0 ;
					for ( j = 0 ; j < No ; j++ )
					{
						Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
						Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
					}
					Cell_x[ CellNo ]	/=	No ;
					Cell_y[ CellNo ]	/=	No ;
					Cell_z[ CellNo ]	/=	No ;

					CellNo++ ;
				}
			}
		}
	} else
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$Node" ) )
			{
				for ( i = 0 ; i < NodeNum ; i++ )
				{
					input >> position[ 0 ] >> position[ 1 ] >> position[ 2 ] ;

					Node_x[ i ]	=	position[ 0 ] * scale ;
					Node_y[ i ]	=	position[ 1 ] * scale ;
					Node_z[ i ]	=	position[ 2 ] * scale ;
				}
			}

			// BCFace Information
			if ( Compare_string( buffer, "$Line" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 6 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					No =	2 ;
					BCFace_Typename[ FaceNo ]	=	physical_Name ;
					BCFace_Type[ FaceNo ]		=	Typename_Type[ BCFace_Typename[ FaceNo ] ] ;

					for ( j = 0 ; j < No ; j++ )
					{
						input >> NodeNo[ j ] ;
						BCFace_Node[ FaceNo ].push_back( NodeNo[ j ] ) ;
						Node_BCFace[ NodeNo[ j ] ].push_back( FaceNo ) ;
					}
					FaceNo++ ;
				}
			}

			// Cell Information
			if ( Compare_string( buffer, "$Triangular" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 4 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					Cell_Typename[ CellNo ]	=	physical_Name ;
					Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

					No 						=	3 ;
					Cell_Form[ CellNo ]		=	4 ;
					Cell_FaceNum[ CellNo ]	=	3 ;

					for ( j = 0 ; j < No ; j++ )
						input >> NodeNo[ j ] ;

					CellNodeCCW( No, Node_x, Node_y, NodeNo.get() ) ;

					Cell_x[ CellNo ]	=	0.0 ;
					Cell_y[ CellNo ]	=	0.0 ;
					Cell_z[ CellNo ]	=	0.0 ;
					for ( j = 0 ; j < No ; j++ )
					{
						Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
						Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
					}
					Cell_x[ CellNo ]	/=	No ;
					Cell_y[ CellNo ]	/=	No ;

					CellNo++ ;
				}
			}

			if ( Compare_string( buffer, "$Quadrilateral" ) )
			{
				for ( i = 0 ; i < MeshType_Num[ 5 ] ; i++ )
				{
					input >> element_PhysGrp >> element_NumTag ;

					physical_Name	=	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_Name	=	"Bulk" ;

					Cell_Typename[ CellNo ]	=	physical_Name ;
					Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

					No 						=	4 ;
					Cell_Form[ CellNo ]		=	5 ;
					Cell_FaceNum[ CellNo ]	=	4 ;

					for ( j = 0 ; j < No ; j++ )
						input >> NodeNo[ j ] ;

					CellNodeCCW( No, Node_x, Node_y, NodeNo.get() ) ;

					Cell_x[ CellNo ]	=	0.0 ;
					Cell_y[ CellNo ]	=	0.0 ;
					Cell_z[ CellNo ]	=	0.0 ;
					for ( j = 0 ; j < No ; j++ )
					{
						Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
						Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

						Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
						Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
					}
					Cell_x[ CellNo ]	/=	No ;
					Cell_y[ CellNo ]	/=	No ;

					CellNo++ ;
				}
			}
		}
	}
	input.clear() ;
	input.close() ;
}

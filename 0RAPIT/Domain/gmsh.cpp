#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>
#include "preprocessingtool.h"
#include "sys_log.h"
#include "gmsh.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
using namespace std ;

/*!

@file msh.cpp Functions for reading GMSH mesh.

*/


/*! \brief Read the information from mesh in GMSH format

To obtain the numbers from mesh in GMSH format, including NodeNum, CellNum, BC_TypeNum, and BCFaceNum.

*/

void GetNumbers_GMSH ( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, string Filename )
{
	int 		i, buffer_int, dummy ;

	ifstream	input ;
	string 		buffer ;

	string 		version ;
	int 		msh_filetype, msh_datasize ;

	int 		physical_Dim, physical_ID ;
	string 		physical_Name ;
	stringstream Elements ;
	vector<int> element_data ;

	input.open ( Filename.c_str ( ) ) ;
	if ( !input )
	{
		cerr << "\n";
		cerr << "GMSH_DATA_READ - Fatal error!\n";
		cerr << "  Could not open input file \"" << Filename << "\"\n";
		exit ( 1 ) ;
	}

	// Initial value to the variable 
	*NodeNum	=	0 ;
	*CellNum 	=	0 ;
	*BC_TypeNum	=	0 ;
	*BCFaceNum 	=	0 ;

	buffer_int = 0 ;
	while ( getline( input, buffer ) )
	{
		if ( Find_string( buffer, "$MeshFormat" ) )
		{
			buffer_int++ ;
			input >>  version >>  msh_filetype >>  msh_datasize ;
			Log().TagDump( logLEVEL4 ) << "\tGMSH file:	" << Filename ;
			Log().TagDump( logLEVEL4 ) << "\tver: " << version << ", msh type: " << msh_filetype << ", msh data size: "<< msh_datasize ;
			break ;
		}
	}

	if ( buffer_int == 0 )
	{
		cout << Filename << ": This is not a MSH mesh!" << endl ;
		exit( -1 ) ;
	}

	if ( dimension == 2 )
	{
		while ( getline( input, buffer ) )
		{

			if ( Compare_string( buffer, "$PhysicalNames" ) )
			{
				getline( input, buffer ) ;
				(*BC_TypeNum)	=	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < (*BC_TypeNum) ; i++ )
				{
					input >> physical_Dim >> physical_ID >> physical_Name ;
					
					if ( physical_Dim == 1 )
					{
						Log().TagDump( logLEVEL4 )  << "\t\tBCFace_Type:	" << physical_Name ;
					} else if ( physical_Dim == 2 )
					{
						Log().TagDump( logLEVEL4 )  << "\t\tCell_Type:	" << physical_Name  ;
					}
				}
				//Log().TagDump( logLEVEL4 )  << "\tBCTypeNum: " << *BC_TypeNum ;
			}

			if ( Compare_string( buffer, "$Nodes" ) )
			{
				getline( input, buffer ) ;
				( *NodeNum )	=	atoi( buffer.c_str() ) ;
				//Log().TagDump( logLEVEL4 )  <<  "\tNodeNum: "  << *NodeNum ;
			}

			if ( Compare_string( buffer, "$Elements" ) )
			{
				getline( input, buffer ) ;
				buffer_int =	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < buffer_int ; i++ )
				{
					Elements.clear() ;
					element_data.clear() ;

					getline( input, buffer ) ;
					Elements << buffer ;

					while( Elements >> dummy )
					{
						element_data.push_back( dummy ) ;
					}

					if ( element_data[ 1 ] == 1 )
					{
						(*BCFaceNum)++ ;
					} else //if ( element_data[ 1 ] == 2 || element_data[ 1 ] == 3 )
					{
						(*CellNum)++ ;
					}
				}
			}
		}
	} else
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$PhysicalNames" ) )
			{
				getline( input, buffer ) ;
				(*BC_TypeNum)	=	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < (*BC_TypeNum) ; i++ )
				{
					input >> physical_Dim >> physical_ID >> physical_Name ;
					
					if ( physical_Dim == 2 )
					{
						Log().TagDump( logLEVEL4 )  << "\t\tBCFace_Type:	" << physical_Name ;
					} else if ( physical_Dim == 3 )
					{
						Log().TagDump( logLEVEL4 )  << "\t\tCell_Type:	" << physical_Name  ;
					}
				}
				//Log().TagDump( logLEVEL4 )  << "\tBCTypeNum: " << *BC_TypeNum ;
			}

			if ( Compare_string( buffer, "$Nodes" ) )
			{
				getline( input, buffer ) ;
				( *NodeNum )	=	atoi( buffer.c_str() ) ;
				//Log().TagDump( logLEVEL4 )  <<  "\tNodeNum: "  << *NodeNum ;
			}

			if ( Compare_string( buffer, "$Elements" ) )
			{
				getline( input, buffer ) ;
				buffer_int =	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < buffer_int ; i++ )
				{
					Elements.clear() ;
					element_data.clear() ;

					getline( input, buffer ) ;
					Elements << buffer ;

					while( Elements >> dummy )
					{
						element_data.push_back( dummy ) ;
					}

					if ( element_data[ 1 ] == 2 || element_data[ 1 ] == 3 )
					{
						(*BCFaceNum)++ ;
					} else if ( element_data[ 1 ] != 1 )
					{
						(*CellNum)++ ;
					}
				}
			}
		}		
	}
	input.clear() ;
	input.close() ;

	Log().TagDump( logLEVEL4 )  << "\tCellNum:	" << (*CellNum) ;
	Log().TagDump( logLEVEL4 )  << "\tNodeNum:	" << (*NodeNum) ;
	Log().TagDump( logLEVEL4 )  << "\tBCFaceNum:	" << (*BCFaceNum)  ;

	(*BC_TypeNum)++ ;
}

void ReadGrid_GMSH ( int dimension, double scale, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int *Cell_Form, int *Cell_FaceNum, int *Cell_Type, string *Cell_Typename, vector<int> *Cell_Node, vector<int> *Node_Cell, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, vector<int> *Node_BCFace, string Filename ) 
{
	int 		i, j, No, buffer_int, FaceNo, CellNo ;
	string 		buffer ;
	ifstream	input ;
	int 		physical_dimenstion, physical_number ;
	int 		element_ID, element_Type, element_NumTags, element_PhysGrp, element_ElemGrp ;
	string 		physical_name ; 
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
		if ( Compare_string( buffer, "$PhysicalNames" ) )
		{
			getline( input, buffer ) ;
			buffer_int	=	atoi( buffer.c_str() ) ;

			for ( i = 1 ; i <= buffer_int ; i++ )
			{
				input >> physical_dimenstion >> physical_number >> physical_name ;
				// remove " 
				physical_name.erase ( std::remove( physical_name.begin(), physical_name.end(), '"'), physical_name.end() ) ;
				Typename[ physical_number ]	=	physical_name ;

				if ( physical_dimenstion == ( dimension - 1 ) )
				{
					if ( check_Typename[ physical_name ] == 0 && physical_name != "Unspecified" )
					{
						(*BC_Typename).push_back( physical_name ) ;

						Typename_Type[ physical_name ]	=	(*BC_TypeNum) ;
						check_Typename[ physical_name ]	=	1 ;
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
		if ( Compare_string( buffer, "$PhysicalNames" ) )
		{
			getline( input, buffer ) ;
			buffer_int	=	atoi( buffer.c_str() ) ;

			for ( i = 1 ; i <= buffer_int ; i++ )
			{
				input >> physical_dimenstion >> physical_number >> physical_name ;
				// remove " 
				physical_name.erase ( std::remove( physical_name.begin(), physical_name.end(), '"'), physical_name.end() ) ;

				if ( physical_dimenstion == dimension )
				{				
					if ( check_Typename[ physical_name ] == 0 && physical_name != "Unspecified" )
					{
						(*BC_Typename).push_back( physical_name ) ;

						Typename_Type[ physical_name ]	=	(*BC_TypeNum) ;
						check_Typename[ physical_name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
				}
			}
		}
	}
	input.clear() ;
	input.close() ;

	input.open ( Filename.c_str ( ) ) ;
	if ( dimension == 2 )
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$Nodes" ) )
			{
				getline( input, buffer ) ;
				buffer_int	=	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < buffer_int ; i++ )
				{
					input >> No >> position[ 0 ] >> position[ 1 ] >> position[ 2 ] ;
					No-- ;

					Node_x[ No ]	=	position[ 0 ] * scale ;
					Node_y[ No ]	=	position[ 1 ] * scale ;
					Node_z[ No ]	=	position[ 2 ] * scale ;
				}
			}

			// example:	1	1	2	3	1	2	5
			if ( Compare_string( buffer, "$Elements" ) )
			{
				getline( input, buffer ) ;
				buffer_int	=	atoi( buffer.c_str() ) ;

				FaceNo	=	0 ;
				CellNo 	=	0 ;
				for ( i = 0 ; i < buffer_int ; i++ )
				{
					input >> element_ID >> element_Type >> element_NumTags >> element_PhysGrp >> element_ElemGrp ;

					physical_name =	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_name =	"Bulk" ;

					if ( element_Type == 1 )	// Bar
					{
						BCFace_Typename[ FaceNo ]	=	physical_name ;
						BCFace_Type[ FaceNo ]		=	Typename_Type[ BCFace_Typename[ FaceNo ] ] ;

						for ( j = 0 ; j < 2 ; j++ )
						{
							input >> No ;
							No-- ;
							BCFace_Node[ FaceNo ].push_back( No ) ;
							Node_BCFace[ No ].push_back( FaceNo ) ;
						}
						FaceNo++ ;
					} else //if ( element_Type == 2 || element_Type == 3 )	// 2 = Triangular, 3 = Quadrilateral
					{
						Cell_Typename[ CellNo ]	=	physical_name ;
						Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

						if ( element_Type == 2 )
						{
							No =	3 ;
							Cell_Form[ CellNo ]		=	4 ;
							Cell_FaceNum[ CellNo ]	=	3 ;
						} else
						{
							No =	4 ;
							Cell_Form[ CellNo ]		=	5 ;
							Cell_FaceNum[ CellNo ]	=	4 ;
						}

						for ( j = 0 ; j < No ; j++ )
						{
							input >> NodeNo[ j ] ;
							NodeNo[ j ]-- ;
						}

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
	} else
	{
		while ( getline( input, buffer ) )
		{
			if ( Compare_string( buffer, "$Nodes" ) )
			{
				getline( input, buffer ) ;
				buffer_int	=	atoi( buffer.c_str() ) ;

				for ( i = 0 ; i < buffer_int ; i++ )
				{
					input >> No >> position[ 0 ] >> position[ 1 ] >> position[ 2 ] ;
					No-- ;

					Node_x[ No ]	=	position[ 0 ] * scale ;
					Node_y[ No ]	=	position[ 1 ] * scale ;
					Node_z[ No ]	=	position[ 2 ] * scale ;
				}
			}

			// example:	1	1	2	3	1	2	5
			if ( Compare_string( buffer, "$Elements" ) )
			{
				getline( input, buffer ) ;
				buffer_int	=	atoi( buffer.c_str() ) ;

				FaceNo	=	0 ;
				CellNo 	=	0 ;
				for ( i = 0 ; i < buffer_int ; i++ )
				{
					input >> element_ID >> element_Type >> element_NumTags >> element_PhysGrp >> element_ElemGrp ;

					physical_name =	Typename[ element_PhysGrp ] ;
					if ( Typename[ element_PhysGrp ] == "Unspecified" )
						physical_name =	"Bulk" ;

					if ( element_Type == 2 || element_Type == 3 ) // 2 = Triangular, 3 = Quadrilateral
					{
						BCFace_Typename[ FaceNo ]	=	physical_name ;
						BCFace_Type[ FaceNo ]		=	Typename_Type[ BCFace_Typename[ FaceNo ] ] ;

						for ( j = 0 ; j < element_Type + 1 ; j++ )
						{
							input >> NodeNo[ j ] ;
							NodeNo[ j ]-- ;

							BCFace_Node[ FaceNo ].push_back( NodeNo[ j ] ) ;
							Node_BCFace[ NodeNo[ j ] ].push_back( FaceNo ) ;					
						}						
						FaceNo++ ;
					} else if ( element_Type != 1 )
					{
						Cell_Typename[ CellNo ]	=	physical_name ;
						Cell_Type[ CellNo ]		=	Typename_Type[ Cell_Typename[ CellNo ] ] ;

						if ( element_Type == 4 )		// Tetra
						{
							No =	4 ;
							Cell_Form[ CellNo ]		=	0 ;
							Cell_FaceNum[ CellNo ]	=	4 ;
						} else if ( element_Type == 5 )	// Hexahedron
						{
							No =	8 ;
							Cell_Form[ CellNo ]		=	1 ;
							Cell_FaceNum[ CellNo ]	=	6 ;
						} else if ( element_Type == 6 )	//	Prism
						{
							No =	6 ;
							Cell_Form[ CellNo ]		=	2 ;
							Cell_FaceNum[ CellNo ]	=	5 ;	
						} else if ( element_Type == 7 )	//	Pyramid
						{
							No =	5 ;
							Cell_Form[ CellNo ]		=	3 ;
							Cell_FaceNum[ CellNo ]	=	5 ;													
						}
	
						Cell_x[ CellNo ]	=	0.0 ;
						Cell_y[ CellNo ]	=	0.0 ;
						Cell_z[ CellNo ]	=	0.0 ;
						for ( j = 0 ; j < No ; j++ )
						{
							input >> NodeNo[ j ] ;
							NodeNo[ j ]-- ;

							Cell_Node[ CellNo ].push_back( NodeNo[ j ] ) ;
							Node_Cell[ NodeNo[ j ] ].push_back( CellNo ) ;

							Cell_x[ CellNo ]	+=	Node_x[ NodeNo[ j ] ] ;
							Cell_y[ CellNo ]	+=	Node_y[ NodeNo[ j ] ] ;
							Cell_z[ CellNo ]	+=	Node_z[ NodeNo[ j ] ] ;	
				
						}
						Cell_x[ CellNo ]	/=	No ;
						Cell_y[ CellNo ]	/=	No ;
						Cell_z[ CellNo ]	/=	No ;

						CellNo++ ;		
					}
				}
			}
		}		
	}
	input.clear() ;
	input.close() ;
}

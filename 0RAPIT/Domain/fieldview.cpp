#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>
#include "fieldview.h"
#include "preprocessingtool.h"
//#include "sys_log.h"

using namespace std ;

/*!

@file fieldview.cpp Functions for mesh in ViewView format.

*/


/*! \brief Read the information from mesh in FieldView format

To obtain the numbers from mesh in FieldView format, including node, cell, BC type, and BC face.

*/

void GetNumbers_FieldView ( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, string Filename )
{
	string		buffer ;
	ifstream	Input ;

	// Initial value to the variable 
	( *NodeNum )	=	0 ;
	( *CellNum )	=	0 ;
	( *BC_TypeNum )	=	0 ;
	( *BCFaceNum )	=	0 ;

	Input.open( Filename.c_str(), ios::in ) ;

	getline( Input, buffer ) ;
	if ( !Find_string( buffer, "FIELDVIEW" ) )
	{
		cout << Filename << ": This is not a FIELDVIEW mesh!" << endl ;
		exit( -1 ) ;
	} else
	{
		while ( getline( Input, buffer ) )
		{
			if ( Compare_string( buffer, "Boundary Table" ) )
			{
				getline( Input, buffer ) ;
				( *BC_TypeNum )	=	atoi( buffer.c_str() ) ;
			}

			if ( Compare_string( buffer, "Nodes" ) )
			{
				getline( Input, buffer ) ;
				( *NodeNum )	=	atoi( buffer.c_str() ) ;
			}

			if ( Compare_string( buffer, "Boundary Faces" ) )
			{
				getline( Input, buffer ) ;
				( *BCFaceNum )	=	atoi( buffer.c_str() ) ;
			}

			// Elements does not provide the total number of elment. Simply, tracing line by line and count. 
			if ( Compare_string( buffer, "Elements" ) )
			{
				while ( getline( Input, buffer ) )
				{
					if ( !Compare_string( buffer, "Variables" ) )
						( *CellNum )++ ;
					else
						break ;
				}
			}
		}
	}

	Input.clear() ;
	Input.close() ;
	if ( dimension == 2 )
		(*NodeNum)	/=	2 ;
	(*BC_TypeNum)++ ;
}

void ReadGrid_FieldView ( int dimension, double scale, int NodeNum, int CellNum, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int *BCFaceNum, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, int *Cell_Form, int *Cell_Type, string *Cell_Typename, int *Cell_FaceNum, vector<int> *Cell_Node, vector<int> *Node_Cell, vector<int> *Node_BCFace, string Filename )
{
	int			tmp, BCNo, BCFace_NodeNum, Type, Num ;
	string		buffer, name ;
	ifstream	Input ;

	boost::shared_array< int >		NodeNo ;
	NodeNo =	boost::shared_array< int > ( new int [ 8 ] ) ;
	for ( int i = 0 ; i < 8 ; i++ )
		NodeNo[ i ] = -999 ;

	boost::shared_array< string >	Typename ;
	Typename =	boost::shared_array< string > ( new string [ *BC_TypeNum ] ) ;
	Typename[ 0 ]	=	"Bulk" ;

	Input.open( Filename.c_str(), ios::in ) ;
	while ( getline( Input, buffer ) )
	{
		if ( Compare_string( buffer, "Boundary Table" ) )
		{
			getline( Input, buffer ) ;
			for ( int i = 1 ; i < *BC_TypeNum ; i++ )
				Input >> buffer >> Typename[ i ] ;
		}
	}
	Input.clear() ;
	Input.close() ;

	map< string, int > 	check_Typename, Typename_Type ;
	Typename_Type[ "Bulk" ]		=	0 ;
	check_Typename[ "Bulk" ]	=	1 ;

	BC_Typename->push_back( "Bulk" ) ;
	(*BC_TypeNum)	=	1 ;

	if ( dimension == 2 )
	{
		int 	Node_NewNode[ NodeNum * 2 ], NewNo = 0, NewBCFaceNum = 0 ;
		double 	z, _Node_z[ NodeNum * 2 ], position[ 3 ] ;

		for ( int i = 0 ; i < NodeNum * 2 ; i++ )
		{
			Node_NewNode[ i ]	=	0 ;
			_Node_z[ i ]		=	0. ;
		}

		for ( int i = 0 ; i < 3  ; i++ )
			position[ i ]	=	0. ;

		Input.open( Filename.c_str(), ios::in ) ;
		while ( getline( Input, buffer ) )
		{
			if ( Compare_string( buffer, "Nodes" ) )
			{
				getline( Input, buffer ) ;
				for ( int i = 0 ; i < ( NodeNum * 2 ) ; i++ )
				{
					Input >> position[ 0 ] >> position[ 1 ] >> position[ 2 ] ;
					
					_Node_z[ i ]	=	position[ 2 ] ;				

					if ( _Node_z[ i ] == 0. )
					{
						Node_NewNode[ i ]	=	NewNo ;

						Node_x[ NewNo ]	=	position[ 0 ] * scale ;
						Node_y[ NewNo ]	=	position[ 1 ] * scale ;
						Node_z[ NewNo ]	=	0. ;
							
						NewNo++ ;
					}
				}
			}

			NewNo =	0 ;
			// example:	1  4  10303  10304  10203  10202
			if ( Compare_string( buffer, "Boundary Faces" ) )
			{
				getline( Input, buffer ) ;
				for ( int i = 0 ; i < *BCFaceNum ; i++ )
				{
					Input >> BCNo >> BCFace_NodeNum ;

					z =	0. ;
					tmp = 0 ;
					for ( int j = 0 ; j < BCFace_NodeNum ; j++ )
					{
						Input >> NodeNo[ j ] ;
						z	+=	_Node_z[ NodeNo[ j ] - 1 ] ;
						tmp +=	Node_NewNode[ NodeNo[ j ] - 1 ] ;
					}

					if ( z != 0. && tmp != 0 )	// BCFace
					{
						name =	Typename[ BCNo ] ;
						BCFace_Typename[ NewBCFaceNum ]	=	name ;

						if ( check_Typename[ name ] == 0 )
						{
							BC_Typename->push_back( name ) ;

							Typename_Type[ name ]	=	(*BC_TypeNum) ;
							check_Typename[ name ]	=	1 ;
							(*BC_TypeNum)++ ;
						}
						BCFace_Type[ NewBCFaceNum ] =	Typename_Type[ name ] ;				

						for ( int j = 0 ; j < BCFace_NodeNum ; j++ )
						{
							if ( _Node_z[ NodeNo[ j ] - 1 ] == 0 )
							{
								BCFace_Node[ NewBCFaceNum ].push_back( Node_NewNode[ NodeNo[ j ] - 1 ] ) ;
								Node_BCFace[ Node_NewNode[ NodeNo[ j ] - 1 ] ].push_back( NewBCFaceNum ) ;
							}
						}
						NewBCFaceNum++ ;
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
		Input.clear() ;
		Input.close() ;

		Input.open( Filename.c_str(), ios::in ) ;
		while ( getline( Input, buffer ) )
		{
			NewNo =	0 ;
			// example:	1  4  10303  10304  10203  10202
			if ( Compare_string( buffer, "Boundary Faces" ) )
			{
				getline( Input, buffer ) ;
				for ( int i = 0 ; i < *BCFaceNum ; i++ )
				{
					Input >> BCNo >> BCFace_NodeNum ;

					z =	0. ;
					tmp = 0 ;
					for ( int j = 0 ; j < BCFace_NodeNum ; j++ )
					{
						Input >> NodeNo[ j ] ;
						z	+=	_Node_z[ NodeNo[ j ] - 1 ] ;
						tmp +=	Node_NewNode[ NodeNo[ j ] - 1 ] ;
					}

					if ( z == 0. )	// Cell
					{
						for ( int j = 0 ; j < BCFace_NodeNum ; j++ )
							NodeNo[ j ]	=	Node_NewNode[ NodeNo[ j ] - 1 ] ;

						CellNodeCCW( BCFace_NodeNum, Node_x, Node_y, NodeNo.get() ) ;

						name =	Typename[ BCNo ] ;

						Cell_Typename[ NewNo ]	=	name ;

						if ( check_Typename[ name ] == 0 )
						{
							BC_Typename->push_back( name ) ;

							Typename_Type[ name ]	=	(*BC_TypeNum) ;
							check_Typename[ name ]	=	1 ;
							(*BC_TypeNum)++ ;
						}
						Cell_Type[ NewNo ]	=	Typename_Type[ name ] ;

						Cell_FaceNum[ NewNo ]	=	BCFace_NodeNum ;
						if ( Cell_FaceNum[ NewNo ] == 3 )
						{
							Cell_Form[ NewNo ] 	=	4 ;
						} else
						{
							Cell_Form[ NewNo ] 	=	5 ;
						}

						Cell_x[ NewNo ]	=	0. ;
						Cell_y[ NewNo ]	=	0. ;
						Cell_z[ NewNo ]	=	0. ;
						for ( int j = 0 ; j < BCFace_NodeNum ; j++ )
						{
							Cell_Node[ NewNo ].push_back( NodeNo[ j ] ) ;
							Node_Cell[ NodeNo[ j ] ].push_back( NewNo ) ;

							Cell_x[ NewNo ]	+=	Node_x[ NodeNo[ j ] ] ;
							Cell_y[ NewNo ]	+=	Node_y[ NodeNo[ j ] ] ;
						}
						Cell_x[ NewNo ]	/=	BCFace_NodeNum ;
						Cell_y[ NewNo ]	/=	BCFace_NodeNum ;

						NewNo++ ;
					}
				}			
			}
		}
		Input.clear() ;
		Input.close() ;

		*BCFaceNum 	=	NewBCFaceNum ;
	} else if ( dimension == 3 )
	{
		Input.open( Filename.c_str(), ios::in ) ;
		while ( getline( Input, buffer ) )
		{
			if ( Compare_string( buffer, "Nodes" ) )
			{
				getline( Input, buffer ) ;
				for ( int i = 0 ; i < NodeNum ; i++ )
				{
					Input >> Node_x[ i ] >> Node_y[ i ] >> Node_z[ i ] ;
					Node_x[ i ]	*=	scale ;
					Node_y[ i ]	*=	scale ;
					Node_z[ i ]	*=	scale ;
				}
			}

			// example:	1  4  10303  10304  10203  10202
			if ( Compare_string( buffer, "Boundary Faces" ) )
			{
				getline( Input, buffer ) ;
				for ( int i = 0 ; i < *BCFaceNum ; i++ )
				{
					Input >> BCFace_Type[ i ] >> BCFace_NodeNum ;

					BCFace_Typename[ i ]	=	Typename[ BCFace_Type[ i ] ] ;

					for ( int j = 0 ; j < BCFace_NodeNum ; j++ )
					{
						Input >> NodeNo[ j ] ;
						BCFace_Node[ i ].push_back( NodeNo[ j ] - 1 ) ;
						Node_BCFace[ NodeNo[ j ] - 1 ].push_back( i ) ;
					}
				}
			}
	
			// example:	2  1  1  2  102  103  10202  10203  10303  10304
			if ( Compare_string( buffer, "Elements" ) )
			{
				for ( int i = 0 ; i < CellNum ; i++ )
				{
					// Cell[ i ].form = Grid_type, first "NodeNo" is nothing
					Input >> Cell_Form[ i ] >> NodeNo[ 0 ] ;
					Cell_Form[ i ]-- ;

					Cell_Type[ i ]		=	0 ;
					Cell_Typename[ i ]	=	Typename[ Cell_Type[ i ] ] ;

					// For Tetra Cell Type
					if ( Cell_Form[ i ] == 0 )
					{
						Num =	4 ;
						for ( int j = 0 ; j < Num ; j++ )
							Input >> NodeNo[ j ] ;

						tmp			=	NodeNo[ 1 ] ;
						NodeNo[ 1 ]	=	NodeNo[ 2 ] ;
						NodeNo[ 2 ]	=	tmp ;

						Cell_FaceNum[ i ]	=	4 ;
					// For Hex. Cell Type
					} else if ( Cell_Form[ i ] == 1 )
					{
						Num = 8 ;
						for ( int j = 0 ; j < Num ; j++ )
							Input >> NodeNo[ j ] ;

						tmp			=	NodeNo[ 2 ] ;
						NodeNo[ 2 ]	=	NodeNo[ 3 ] ;
						NodeNo[ 3 ]	=	tmp ;

						tmp			=	NodeNo[ 6 ] ;
						NodeNo[ 6 ]	=	NodeNo[ 7 ] ;
						NodeNo[ 7 ]	=	tmp ;

						Cell_FaceNum[ i ]	=	6 ;
					// For Prism Cell Type
					} else if ( Cell_Form[ i ] == 2 )
					{
						Num =	6 ;
						for ( int j = 0 ; j < Num ; j++ )
							Input >> NodeNo[ j ] ;

						tmp			=	NodeNo[ 4 ] ;
						NodeNo[ 4 ]	=	NodeNo[ 5 ] ;
						NodeNo[ 5 ]	=	tmp ;

						Cell_FaceNum[ i ]	=	5 ;
					} else if ( Cell_Form[ i ] == 3 )
					{
						Num =	5 ;
						for ( int j = 0 ; j < Num ; j++ )
							Input >> NodeNo[ j ] ;

						Cell_FaceNum[ i ]	=	5 ;
					}

					Cell_x[ i ]	=	0. ;
					Cell_y[ i ]	=	0. ;
					Cell_z[ i ]	=	0. ;
					for ( int j = 0 ; j < Num ; j++ )
					{
						Cell_Node[ i ].push_back( NodeNo[ j ] - 1 ) ;
						Node_Cell[ NodeNo[ j ] - 1 ].push_back( i ) ;

						Cell_x[ i ]	+=	Node_x[ NodeNo[ j ] - 1 ] ;
						Cell_y[ i ]	+=	Node_y[ NodeNo[ j ] - 1 ] ;
						Cell_z[ i ]	+=	Node_z[ NodeNo[ j ] - 1 ] ;
					}
					Cell_x[ i ]	/=	Num ;
					Cell_y[ i ]	/=	Num ;
					Cell_z[ i ]	/=	Num ;
				}
				break ;
			}
		}
		Input.clear() ;
		Input.close() ;
	}
}
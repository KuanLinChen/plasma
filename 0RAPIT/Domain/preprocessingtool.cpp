#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include "sys_log.h"
#include "preprocessingtool.h"

using namespace std ;

int Get_filetype ( string _n )
{
	if ( _n.substr(  _n.find_last_of( '.' ) + 1 ) == "uns" )
	{
		return FormatFIELDVIEW ;
	} else if ( _n.substr(  _n.find_last_of( '.' ) + 1 ) == "cgns" )
	{
		return FormatCGNS ;
	} else if ( _n.substr(  _n.find_last_of( '.' ) + 1 ) == "msh" )
	{
		return FormatMSH ;
	} else if ( _n.substr(  _n.find_last_of( '.' ) + 1 ) == "pti" )
	{
		return FormatPTI ;
	}
	return -999 ;
}

bool Compare_string ( string s1, string s2 )
{
	bool	flag = true ;

	for ( int i = 0 ; i < s2.size() ; i++ )
	{
		if ( s1[ i ] != s2[ i ] )
		{
			flag	=	false ;
			i		=	s2.size( ) ;
		}
	}

	return	flag ;
}


bool Find_string ( string s1, string s2 )
{
	size_t	found = s1.find( s2 ) ;

	if ( found != std::string::npos )
	{
    	return true ;
	} else
	{
		return false ;
	}
}

void CellNodeCCW( int NodeNum, double *Node_x, double *Node_y, int *NodeNo )
{
	int		NodeNoMaxx, NodeBuffer, No ;
	double	x0, y0, x1, y1, a, b, c ;

	boost::shared_array<double>	Angle ;
	Angle =	boost::shared_array<double> ( new double [ 4 ] ) ;
	for ( int i = 0 ; i < 4 ; i++ )
		Angle[ i ]	=	0. ;

	NodeNoMaxx	=	0 ;

	for ( int i = 1 ; i < NodeNum ; i++ )
	{
		if ( Node_x[ NodeNo[ i ] ] > Node_x[ NodeNo[ NodeNoMaxx ] ] )
		{
			NodeNoMaxx	=	i ;
		}
	}

	NodeBuffer				=	NodeNo[ 0 ] ;
	NodeNo[ 0 ]				=	NodeNo[ NodeNoMaxx ] ;
	NodeNo[ NodeNoMaxx ]	=	NodeBuffer ;

	x0	=	Node_x[ NodeNo[ 0 ] ] ;
	y0	=	Node_y[ NodeNo[ 0 ] ] ;
	for ( int i = 1 ; i < NodeNum ; i++ )
	{
		No 	=	NodeNo[ i ] ;
		x1	=	Node_x[ No ] ;
		y1	=	Node_y[ No ] ;
		a	=	x1 - x0 ;
		b	=	y1 - y0 ;
		c	=	sqrt ( a * a + b * b ) ;

		Angle[ i ]	=	CalculateAngleDeg( a, b, c ) ;
	}

	// Sort the angle from the smallest to the biggest.
	BubbleSort( &Angle[ 1 ], &NodeNo[ 1 ], NodeNum - 1 ) ;
}

double CalculateAngleDeg( double a, double b, double c )
{
	double	pi = 3.141592654 ;
	double	angle ;

	angle	=	asin( b / c ) ;
	angle	=	angle / pi * 180. ;
	if ( a < 0. )
		angle	=	180. - angle ;
	else if ( angle < 0. )
		angle	=	360. + angle ;

	return	angle ;
}

template <class TYPE>
void BubbleSort( TYPE *a, int *Index, int Num )
{
	for ( int i = 0 ; i < ( Num - 1 ) ; i++ )
	{
		for ( int j = 0 ; j < ( Num - ( 1 + i ) ) ; j++ )
		{
			if ( a[ j ] > a[ j + 1 ] )
			{
				Swap( &a[ j ], &a[ j + 1 ] ) ;
				Swap( &Index[ j ], &Index[ j + 1 ] ) ;
			}
		}
	}
}

template <class TYPE>
void Swap( TYPE *a, TYPE *b )
{
	TYPE	buffer ;
	buffer	=	( *a ) ;
	( *a )	=	( *b ) ;
	( *b )	=	buffer ;
}

// For Cell-Cell relation
bool CheckNeighbor( int TypeA, int TypeB, vector<int> *_NodeA, vector<int> *_NodeB, int SurfaceA, int SurfaceB, CellMapping *pMapping )
{
	bool	Check	= false ;
	int 	Num = 0, NodeA[ 4 ], NodeB ;

	for ( int i = 0 ; i < 4 ; i++ )
		NodeA[ i ]	=	-999 ;

	for ( int i = 0 ; i < pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] ; i++ )
	{
		NodeA[ i ]	=	pMapping->Node[ TypeA ][ SurfaceA ][ i ] ;
		NodeA[ i ]	=	(*_NodeA)[ NodeA[ i ] ] ;

		for ( int j = 0 ; j < pMapping->SurfaceNodeNum[ TypeB ][ SurfaceB ] ; j++ )
		{
			NodeB	=	pMapping->Node[ TypeB ][ SurfaceB ][ j ] ;
			NodeB	=	(*_NodeB)[ NodeB ] ;

			if ( NodeA[ i ] == NodeB )
			{
				Num++ ;
				break ;
			}
		}
	}

	if ( Num == pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] )
	{
		Check = true ;
	}

	return	Check ;
}

bool CheckNeighbor( int TypeA, int TypeB, vector<int> *_NodeA, vector<int> *_NodeB, int SurfaceA, int SurfaceB, vector<int> *_Face_Node, CellMapping *pMapping )
{
	bool	Check	= false ;
	int 	Num = 0, NodeA[ 4 ], NodeB ;

	for ( int i = 0 ; i < 4 ; i++ )
		NodeA[ i ]	=	-999 ;

	for ( int i = 0 ; i < pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] ; i++ )
	{
		NodeA[ i ]	=	pMapping->Node[ TypeA ][ SurfaceA ][ i ] ;
		NodeA[ i ]	=	(*_NodeA)[ NodeA[ i ] ] ;

		for ( int j = 0 ; j < pMapping->SurfaceNodeNum[ TypeB ][ SurfaceB ] ; j++ )
		{
			NodeB	=	pMapping->Node[ TypeB ][ SurfaceB ][ j ] ;
			NodeB	=	(*_NodeB)[ NodeB ] ;

			if ( NodeA[ i ] == NodeB )
			{
				Num++ ;
				break ;
			}
		}
	}

	if ( Num == pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] )
	{
		for ( int i = 0 ; i < Num ; i++ )
			_Face_Node->push_back( NodeA[ i ] ) ;

		Check = true ;
	}

	return	Check ;
}

// InnerFace
bool CheckNeighbor( int TypeA, vector<int> *_NodeA, int SurfaceA, vector<int> *_NodeB, CellMapping *pMapping )
{
	bool	Check	= false ;
	int 	Num = 0, NodeA[ 4 ], NodeB ;

	for ( int i = 0 ; i < 4 ; i++ )
		NodeA[ i ]	=	-999 ;

	for ( int i = 0 ; i < pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] ; i++ )
	{
		NodeA[ i ]	=	pMapping->Node[ TypeA ][ SurfaceA ][ i ] ;
		NodeA[ i ]	=	(*_NodeA)[ NodeA[ i ] ] ;

		for ( int j = 0 ; j < (*_NodeB).size() ; j++ )
		{
			NodeB	=	(*_NodeB)[ j ] ;

			if ( NodeA[ i ] == NodeB )
			{
				Num++ ;
				break ;
			}
		}
	}

	if ( Num == pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] )
	{
		Check = true ;
	}

	return	Check ;
}


bool CheckNeighbor( int TypeA, vector<int> *_NodeA, int SurfaceA, vector<int> *_NodeB, vector<int> *_Face_Node, CellMapping *pMapping )
{
	bool	Check	= false ;
	int 	Num = 0, NodeA[ 4 ], NodeB ;

	for ( int i = 0 ; i < 4 ; i++ )
		NodeA[ i ]	=	-999 ;

	for ( int i = 0 ; i < pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] ; i++ )
	{
		NodeA[ i ]	=	pMapping->Node[ TypeA ][ SurfaceA ][ i ] ;
		NodeA[ i ]	=	(*_NodeA)[ NodeA[ i ] ] ;

		for ( int j = 0 ; j < (*_NodeB).size() ; j++ )
		{
			NodeB	=	(*_NodeB)[ j ] ;

			if ( NodeA[ i ] == NodeB )
			{
				Num++ ;
				break ;
			}
		}
	}

	if ( Num == pMapping->SurfaceNodeNum[ TypeA ][ SurfaceA ] )
	{
		for ( int i = 0 ; i < Num ; i++ )
		{
			//_Face_Node->push_back( NodeA[ i ] ) ;
			_Face_Node->push_back( NodeA[ i ] ) ;
		}

		Check = true ;
	}

	return	Check ;
}
/*
void CreateCellCellRelation( int dimension, int NodeNum, int CellNum, int *InnerFaceNum, int *InterFaceNum, vector<int> *Node_Cell, int *Cell_Type, vector<int> *Cell_Node, vector<int> *Cell_Cell )
{
	int No, k, Num = 0, CheckNum ;
	map<int, int> C ; // binary search tree

	if ( dimension == 2 )
	{
		CheckNum = 2 ;
	} else
	{
		CheckNum = 4 ;
	}

	(*InnerFaceNum)	=	0 ;
	(*InterFaceNum)	=	0 ;
	for ( int i = 0 ; i < CellNum ; i++ )
	{
		for ( vector<int>::iterator j = Cell_Node[ i ].begin() ; j != Cell_Node[ i ].end(); j++ )
		{
			No =	*j ;
			for ( vector<int>::iterator k = Node_Cell[ No ].begin() ; k != Node_Cell[ No ].end(); k++ )
			{
				C[ *k ]++ ;
			}
		}

		for ( auto k: C )
		{
			if ( k.second >= CheckNum && k.first != i )
			{
				Cell_Cell[ i ].push_back( k.first ) ;
				if ( Cell_Type[ i ] != Cell_Type[ k.first ] )
					(*InterFaceNum)++ ;
			}
		}
		C.clear() ;

		*InnerFaceNum +=	Cell_Cell[ i ].size() ;
	}

	(*InnerFaceNum)	/=	2 ;
	(*InterFaceNum) /=	2 ;
}
*/

void CreateCellCellRelation( int NodeNum, int CellNum, int *InnerFaceNum, int *InterFaceNum, int *FaceNum, vector<int> *Node_Cell, int *Cell_Form, int *Cell_Type, vector<int> *Cell_Node, vector<int> *Cell_Cell )
{
	int	i, j, c1, c2, j1, j2, k1, k2, Type1, Type2 ;

	/*--- Modifyed by Kuan-Lin ---*/
	//int Check[ CellNum ][ 6 ] ;
	int **Check = new int *[ CellNum ] ;
	for ( i=0; i < CellNum ; i++ ) Check[ i ] = new int [6] ;
	/*----------------------------*/
		
	CellMapping	pMapping ;
	for ( i = 0 ; i < CellNum ; i++ )
	{
		for ( j = 0 ; j < 6 ; j++ )
		{
			Check[ i ][ j ]	=	-999 ;
		}
	}

	//cout<<"FFFFUUUUUCCCCCCKKKKK"<<endl;
	(*InnerFaceNum)	=	0 ;
	(*InterFaceNum)	=	0 ;
	for ( i = 0 ; i < NodeNum ; i++ )
	{
		for ( j1 = 0 ; j1 < Node_Cell[ i ].size() ; j1++ )
		{
			c1	=	Node_Cell[ i ][ j1 ] ;

			Type1	=	Cell_Form[ c1 ] ;
			for ( k1 = 0 ; k1 < pMapping.SurfaceNum[ Type1 ] ; k1++ )
			{
				if ( Check[ c1 ][ k1 ] == -999 )
				{
					j2	= 0 ;
					while ( j2 < Node_Cell[ i ].size() )
					{
						if ( j1 != j2  )
						{
							c2	=	Node_Cell[ i ][ j2 ] ;

							Type2	=	Cell_Form[ c2 ] ;
							for ( k2 = 0 ; k2 < pMapping.SurfaceNum[ Type2 ] ; k2++ )
							{
								if ( CheckNeighbor( Type1, Type2, &Cell_Node[ c1 ], &Cell_Node[ c2 ], k1, k2, &pMapping ) )
								{
									Cell_Cell[ c1 ].push_back( c2 ) ;
									Cell_Cell[ c2 ].push_back( c1 ) ;

									(*InnerFaceNum)++ ;

									if ( Cell_Type[ c1 ] != Cell_Type[ c2 ] )
										(*InterFaceNum)++ ;

									Check[ c1 ][ k1 ]	=	1 ;
									Check[ c2 ][ k2 ]	=	1 ;

									j2	=	Node_Cell[ i ].size() ;
									break ;
								}
							}
						}

						j2++ ;
					}
				}
			}
		}
	}
	//cout<<"&FFFFUUUUUCCCCCCKKKKK"<<endl<<endl;
	(*FaceNum)	=	(*InnerFaceNum) ;
	for ( i = 0 ; i < CellNum ; i++ )
	{
		Type1 =	Cell_Form[ i ] ;
		for ( j = 0 ; j < pMapping.SurfaceNum[ Type1 ] ; j++ )
		{
			if ( Check[ i ][ j ] == -999 )
			{
				(*FaceNum)++ ;
			}
		}
	}
}

void CalculateCellVolume( string geometry, int CellNum, int *Cell_Form, double *min_cell_length, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, vector<int> *Cell_Node, double *Cell_Volume, CellMapping *pMapping )
{
	int		Type, NodeNo[ 2 ] ;
	double	x[ 3 ], y[ 3 ], z[ 3 ], cell_length ;

	(*min_cell_length)	=	1.0e+10 ;

	if ( geometry == "2D" )
	{
		for ( int i = 0 ; i < CellNum ; i++ )
		{
			Type 				=	Cell_Form[ i ] ;
			Cell_Volume[ i ]	=	0. ;

			for ( int j = 0 ; j < pMapping->NodeNum[ Type ] ; j++ )
			{
				NodeNo[ 0 ] =	Cell_Node[ i ][ j ] ;
				x[ 0 ]		= 	Node_x[ NodeNo[ 0 ] ] ;
				y[ 0 ]		= 	Node_y[ NodeNo[ 0 ] ] ;

				if ( j == ( pMapping->NodeNum[ Type ] - 1 ) )
				{
					NodeNo[ 1 ] =	Cell_Node[ i ][ 0 ] ;
				} else
				{
					NodeNo[ 1 ] = 	Cell_Node[ i ][ j + 1 ] ;
				}

				x[ 1 ]	= 	Node_x[ NodeNo[ 1 ] ] ;
				y[ 1 ]	= 	Node_y[ NodeNo[ 1 ] ] ;

				Cell_Volume[ i ]	+= ( ( x[ 0 ] * y[ 1 ] ) - ( y[ 0 ] * x[ 1 ] ) ) ;
			}
			Cell_Volume[ i ]	/=	2.0 ;

			//cell_length =	pow( Cell_Volume[ i ], 1.0 / 2.0 ) ;
			cell_length =	sqrt( Cell_Volume[ i ] ) ;
			if ( cell_length < (*min_cell_length) )
				(*min_cell_length)	=	cell_length ;
		}
	} else if ( geometry == "3D" )
	{
		for ( int i = 0 ; i < CellNum ; i++ )
		{
			Type 				=	Cell_Form[ i ] ;
			Cell_Volume[ i ]	=	0. ;
			for ( int j = 0 ; j < pMapping->TetraCellNum[ Type ] ; j++ )
			{
				NodeNo[ 0 ]	=	pMapping->TetraCell[ Type ][ j ][ 0 ] ;
				NodeNo[ 0 ]	=	Cell_Node[ i ][ NodeNo[ 0 ] ] ;

				for ( int k = 0 ; k < 3 ; k++ )
				{
					NodeNo[ 1 ]	=	pMapping->TetraCell[ Type ][ j ][ k + 1 ] ;
					NodeNo[ 1 ]	=	Cell_Node[ i ][ NodeNo[ 1 ] ] ;

					x[ k ]	=	Node_x[ NodeNo[ 1 ] ] - Node_x[ NodeNo[ 0 ] ] ;
					y[ k ]	=	Node_y[ NodeNo[ 1 ] ] - Node_y[ NodeNo[ 0 ] ] ;
					z[ k ]	=	Node_z[ NodeNo[ 1 ] ] - Node_z[ NodeNo[ 0 ] ] ;
				}

				Cell_Volume[ i ]	+=	fabs( ( x[ 0 ] * ( y[ 1 ] * z[ 2 ] - z[ 1 ] * y[ 2 ] )
									+	y[ 0 ] * ( z[ 1 ] * x[ 2 ] - x[ 1 ] * z[ 2 ] )
									+	z[ 0 ] * ( x[ 1 ] * y[ 2 ] - y[ 1 ] * x[ 2 ] ) ) / 6. ) ;

				cell_length =	pow( Cell_Volume[ i ], 1.0 / 3.0 ) ;
				if ( cell_length < (*min_cell_length) )
					(*min_cell_length)	=	cell_length ;
			}
		}
	} else if ( geometry == "Axisymmetric_X" )
	{
		for ( int i = 0 ; i < CellNum ; i++ )
		{
			Type 				=	Cell_Form[ i ] ;
			Cell_Volume[ i ]	=	0. ;
			for ( int j = 0 ; j < pMapping->NodeNum[ Type ] ; j++ )
			{
				NodeNo[ 0 ] =	Cell_Node[ i ][ j ] ;
				x[ 0 ]		= Node_x[ NodeNo[ 0 ] ] ;
				y[ 0 ]		= Node_y[ NodeNo[ 0 ] ] ;

				if ( j == ( pMapping->NodeNum[ Type ] - 1 ) )
				{
					NodeNo[ 1 ]	=	Cell_Node[ i ][ 0 ] ;
				} else
				{
					NodeNo[ 1 ]	=	Cell_Node[ i ][ j + 1 ] ;
				}
				x[ 1 ]	= Node_x[ NodeNo[ 1 ] ] ;
				y[ 1 ]	= Node_y[ NodeNo[ 1 ] ] ;

				Cell_Volume[ i ]	+= ( ( x[ 0 ] * y[ 1 ] ) - ( y[ 0 ] * x[ 1 ] ) ) ;
			}
			Cell_Volume[ i ]	*=	M_PI * Cell_y[ i ] ;

			cell_length =	pow( Cell_Volume[ i ] / ( 2.0 * M_PI * Cell_y[ i ] ), 0.5 ) ;
			if ( cell_length < (*min_cell_length) )
				(*min_cell_length)	=	cell_length ;
		}
	} else if ( geometry == "Axisymmetric_Y" )
	{
		for ( int i = 0 ; i < CellNum ; i++ )
		{
			Type 				=	Cell_Form[ i ] ;
			Cell_Volume[ i ]	=	0. ;
			for ( int j = 0 ; j < pMapping->NodeNum[ Type ] ; j++ )
			{
				NodeNo[ 0 ] =	Cell_Node[ i ][ j ] ;
				x[ 0 ]		= Node_x[ NodeNo[ 0 ] ] ;
				y[ 0 ]		= Node_y[ NodeNo[ 0 ] ] ;

				if ( j == ( pMapping->NodeNum[ Type ] - 1 ) )
				{
					NodeNo[ 1 ]	=	Cell_Node[ i ][ 0 ] ;
				} else
				{
					NodeNo[ 1 ]	=	Cell_Node[ i ][ j + 1 ] ;
				}
				x[ 1 ]	= Node_x[ NodeNo[ 1 ] ] ;
				y[ 1 ]	= Node_y[ NodeNo[ 1 ] ] ;

				Cell_Volume[ i ]	+= ( ( x[ 0 ] * y[ 1 ] ) - ( y[ 0 ] * x[ 1 ] ) ) ;
			}
			Cell_Volume[ i ]	*=	M_PI * Cell_x[ i ] ;

			cell_length =	pow( Cell_Volume[ i ] / ( 2.0 * M_PI * Cell_x[ i ] ), 0.5 ) ;
			if ( cell_length < (*min_cell_length) )
				(*min_cell_length)	=	cell_length ;
		}
	}
}

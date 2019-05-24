#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cgnslib.h>
#include <map>
#include "cgns.h"
#include "preprocessingtool.h"
//#include "sys_log.h"

using namespace std ;

/*!

@file cgns.cpp Functions for mesh in CGNS format.

*/


/*! \brief Read the information from mesh in CGNS format

To obtain the numbers from mesh in CGNS format, including NodeNum, CellNum, BC_TypeNum, and BCFaceNum.


*/
void GetNumbers_CGNS ( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, int *MeshType_Num, string Filename )
{
	bool			flg_debug = 0 ; 

	cgsize_t 		isize[ 1 ][ 3 ] ;
	int 			index_file ;
	char 			zonename[ 100 ],  sectionname[ 100 ] ;	
	ElementType_t	itype ;
	int				nsections, index_sect, nbndry, iparent_flag ;
	cgsize_t		istart, iend ;

	for ( int i = 0 ; i < 3 ; i++ )
		isize[ 0 ][ i ]	=	0 ;

	// Open CGNS file for read-only
	if ( cg_open( Filename.c_str(), CG_MODE_READ, &index_file ) )
		cg_error_exit() ;

	cg_zone_read( index_file, 1, 1, zonename, isize[ 0 ] ) ;

	cg_nsections( index_file, 1, 1, &nsections ) ;

	for ( index_sect = 1 ; index_sect <= nsections ; index_sect++ )
	{
		cg_section_read( index_file, 1, 1, index_sect, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag ) ;

		if ( flg_debug ) 
		{
			cout << "\n\tSection = " << index_sect << endl ;
			cout << "\t\tsection name = " << sectionname << endl ;
			cout << "\t\tsection type = " << ElementTypeName[ itype ] << endl ;
			cout << "\t\tistart = " << ( int ) istart << ", iend = " << ( int ) iend << endl ;
		}

		if ( itype == BAR_2 )
		{
			MeshType_Num[ 6 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else if ( itype == TRI_3 )
		{
			MeshType_Num[ 4 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else if ( itype == QUAD_4 )
		{
			MeshType_Num[ 5 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else if ( itype == TETRA_4 )
		{
			MeshType_Num[ 0 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else if ( itype == PYRA_5 )
		{
			MeshType_Num[ 3 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else if ( itype == PENTA_6 )
		{
			MeshType_Num[ 2 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else if ( itype == HEXA_8 )
		{
			MeshType_Num[ 1 ] +=	( int ) iend - ( int ) istart + 1 ;
		} else
		{
			cout << "\t\tNot reading element data for this element type" << endl ;
		}
	}

	if ( dimension == 2 )
	{
		if ( ( MeshType_Num[ 0 ] + MeshType_Num[ 1 ] + MeshType_Num[ 2 ] + MeshType_Num[ 3 ] ) != 0 )
		{
			cout << "***   This a 3D problem. Please check dimension setting   ***" << endl ;
			exit( -1 ) ;
		}

		(*BCFaceNum)	=	MeshType_Num[ 6 ] ;
	} else
	{
		if ( ( MeshType_Num[ 0 ] + MeshType_Num[ 1 ] + MeshType_Num[ 2 ] + MeshType_Num[ 3 ] ) == 0 )
		{
			cout << "***   This a 2D problem. Please check dimension setting   ***" << endl ;
			exit( -1 ) ;
		}

		(*BCFaceNum)		=	MeshType_Num[ 4 ] + MeshType_Num[ 5 ] ;
		MeshType_Num[ 6 ]	=	0 ;
	}

	// find out number of BCs that exist under this zone
	cg_nbocos( index_file, 1, 1, BC_TypeNum ) ;

	cg_close( index_file ) ;

	(*NodeNum) =	isize[ 0 ][ 0 ] ;
	(*CellNum) =	isize[ 0 ][ 1 ] ;
	(*BC_TypeNum)++ ;
}

void ReadGrid_CGNS ( int dimension, double scale, int NodeNum, int CellNum, int *MeshType_Num, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypenNum_forFace, vector<string> *BC_Typename, int BCFaceNum, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, int *Cell_Form, int *Cell_Type, string *Cell_Typename, int *Cell_FaceNum, vector<int> *Cell_Node, vector<int> *Node_Cell, vector<int> *Node_BCFace, string Filename )
{
	bool			flg_debug = 0 ;

	int				index_file ;
	cgsize_t		irmin, irmax, istart, iend ;
	int				nsections, index_sect, nbndry, iparent_flag ;
	cgsize_t		iparentdata ;
	char			sectionname[ 100 ] ;
	ElementType_t	itype ;
	int				nbocos, ibc ;
	int				normalindex[ 3 ], ndataset ;
	int				normallist ;
	char			boconame[ 100 ] ;
	BCType_t		ibocotype ;
	PointSetType_t	iptsettype ;
	DataType_t		normaldatatype ;
	GridLocation_t	igr ;
	cgsize_t		npts, normallistflag ;
	cgsize_t		*ipnts ;

	int i, j, tmp[ 3 ] ;
	map< int, int > BCNo ;

	for ( i = 0 ; i < 3 ; i++ )
		tmp[ i ]	=	0 ;

	boost::shared_array< string >	Typename ;
	Typename =	boost::shared_array< string > ( new string [ (*BC_TypeNum) ] ) ;
	Typename[ 0 ]	=	"Bulk" ;

	// Open CGNS file for read-only
	cg_open( Filename.c_str(), CG_MODE_READ, &index_file ) ;

	// lower & upper range index of vertices
	irmin	=	1 ;
	irmax	=	NodeNum ;

	// Read grid coordinates
	cg_coord_read( index_file, 1, 1, "CoordinateX", RealDouble, &irmin, &irmax, Node_x ) ;
	cg_coord_read( index_file, 1, 1, "CoordinateY", RealDouble, &irmin, &irmax, Node_y ) ;
	cg_coord_read( index_file, 1, 1, "CoordinateZ", RealDouble, &irmin, &irmax, Node_z ) ;

	for ( i = 0 ; i < NodeNum ; i++ )
	{
		Node_x[ i ] *=	scale ;
		Node_y[ i ] *=	scale ;
		Node_z[ i ] *=	scale ;
	}

	// Get the boundary conditions
	for ( ibc = 1 ; ibc < (*BC_TypeNum) ; ibc++ )
	{
		// find out what BC grid location is ( expecting FaceCenter )
		cg_goto( index_file, 1, "Zone_t", 1, "ZoneBC_t", 1, "BC_t", ibc, "end" ) ;
		cg_gridlocation_read( &igr ) ;

		// get BC info
		cg_boco_info( index_file, 1, 1, ibc, boconame, &ibocotype, &iptsettype, &npts, normalindex, &normallistflag, &normaldatatype, &ndataset ) ;
		
		// read point list in here
		ipnts	=	new cgsize_t[ npts ] ;

		cg_boco_read( index_file, 1, 1, ibc, ipnts, &normallist ) ;

		// Create the mesh & BCs relationship
		for ( i = 0 ; i < npts ; i++ )
			BCNo[ ( int ) ipnts[ i ] ]	=	ibc ;

		Typename[ ibc ]	=	boconame ;

		if ( flg_debug )
		{
			cout << "\nBC number: " << ibc
				 << "\n\tname = " << boconame
				 << "\n\tNo. of elements = " << ( int ) npts << endl ;
			for ( int i = 0 ; i < npts ; i++ )
				cout << "\t\tipnts[ " << i << " ] = " << ( int ) ipnts[ i ] << endl ;
		}

		delete [] ipnts ;
	}

	// declaration for element connectivity 
	boost::shared_array<int>	_Node ;
	_Node	=	boost::shared_array<int> ( new int [ 4 ] ) ;

	for( int i = 0 ; i < 4 ; i++ )
		_Node[ i ]	=	-999 ;

	int Num[ 7 ], F_Num = 0, C_Num = 0, NodeNo ;

	for ( i = 0 ; i < 7 ; i++ )
		Num[ i ] =	0 ;

	cgsize_t	ielem_TETRA_4[ MeshType_Num[ 0 ] ][ 4 ] ;
	cgsize_t	ielem_HEXA_8[ MeshType_Num[ 1 ] ][ 8 ] ;
	cgsize_t	ielem_PENTA_6[ MeshType_Num[ 2 ] ][ 6 ] ;
	cgsize_t	ielem_PYRA_5[ MeshType_Num[ 3 ] ][ 5 ] ;
	cgsize_t	ielem_TRI_3[ MeshType_Num[ 4 ] ][ 3 ] ;
	cgsize_t	ielem_QUAD_4[ MeshType_Num[ 5 ] ][ 4 ] ;
	cgsize_t	ielem_BAR_2[ MeshType_Num[ 6 ] ][ 2 ] ;

	string name ;
	map< string, int > 	check_Typename, Typename_Type ;

	Typename_Type[ "Bulk" ]		=	0 ;
	check_Typename[ "Bulk" ]	=	1 ;

	BC_Typename->push_back( "Bulk" ) ;
	(*BC_TypeNum)	=	1 ;

	// find out how many sections
	cg_nsections( index_file, 1, 1, &nsections ) ;

	// read element connectivity
	if ( dimension == 2 )
	{	
		for ( index_sect = 1 ; index_sect <= nsections ; index_sect++ )
		{			
			cg_section_read( index_file, 1, 1, index_sect, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag ) ;

			tmp[ 0 ]	=	( int ) iend - ( int ) istart + 1 ;
			tmp[ 1 ]	=	( int ) istart ;

			if ( itype == BAR_2 )
			{			
				cg_elements_read( index_file, 1, 1, index_sect, ielem_BAR_2[ Num[ 6 ] ], &iparentdata ) ;

				for ( i = Num[ 6 ] ; i < Num[ 6 ] + tmp[ 0 ] ; i++ )
				{
					name 						=	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					BCFace_Typename[ F_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;

						Typename_Type[ name ]	=	(*BC_TypeNum);
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					BCFace_Type[ F_Num ] =	Typename_Type[ name ] ;	

					for ( j = 0 ; j < 2 ; j++ )
					{
						NodeNo =	ielem_BAR_2[ i ][ j ] - 1 ;
						BCFace_Node[ F_Num ].push_back( NodeNo ) ;
						Node_BCFace[ NodeNo ].push_back( F_Num ) ;
					}
					tmp[ 1 ]++ ;
					F_Num++ ;
				}

				Num[ 6 ]	+= tmp[ 0 ] ;
			}
		}

		if ( check_Typename[ "Interface" ] == 0 )
		{
			BC_Typename->push_back( "Interface" ) ;

			Typename_Type[ "Interface" ]	=	(*BC_TypeNum);
			check_Typename[ "Interface" ]	=	1 ;
			(*BC_TypeNum)++ ;
		}
		(*BC_TypenNum_forFace)	=	(*BC_TypeNum) ;

		for ( index_sect = 1 ; index_sect <= nsections ; index_sect++ )
		{
			cg_section_read( index_file, 1, 1, index_sect, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag ) ;

			tmp[ 0 ]	=	( int ) iend - ( int ) istart + 1 ;
			tmp[ 1 ]	=	( int ) istart ;

			if ( itype == TRI_3 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_TRI_3[ Num[ 4 ] ], &iparentdata ) ;

				for ( i = Num[ 4 ] ; i < Num[ 4 ] + tmp[ 0 ] ; i++ )
				{
					Cell_Form[ C_Num ]		=	4 ;
					Cell_FaceNum[ C_Num ]	=	3 ;

					name 					=	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					Cell_Typename[ C_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;

						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					Cell_Type[ C_Num ]	=	Typename_Type[ name ] ;

					for ( j = 0 ; j < 3 ; j++ )
						_Node[ j ] =	ielem_TRI_3[ i ][ j ] - 1 ;

					CellNodeCCW( 3, Node_x, Node_y, _Node.get() ) ;

					Cell_x[ C_Num ]	=	0. ;
					Cell_y[ C_Num ]	=	0. ;
					Cell_z[ C_Num ]	=	0. ;
					for ( j = 0 ; j < 3 ; j++ )
					{
						NodeNo =	_Node[ j ] ;
						Cell_Node[ C_Num ].push_back( NodeNo ) ;
						Node_Cell[ NodeNo ].push_back( C_Num ) ;

						Cell_x[ C_Num ] +=	Node_x[ NodeNo ] ;
						Cell_y[ C_Num ] +=	Node_y[ NodeNo ] ;
					}
					Cell_x[ C_Num ]	/= 3 ;
					Cell_y[ C_Num ]	/= 3 ;

					tmp[ 1 ]++ ;
					C_Num++ ;
				}
				Num[ 4 ]	+=	tmp[ 0 ] ;
			} else if ( itype == QUAD_4 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_QUAD_4[ Num[ 5 ] ], &iparentdata ) ;

				for ( i = Num[ 5 ] ; i < Num[ 5 ] + tmp[ 0 ] ; i++ )
				{
					Cell_Form[ C_Num ]		=	5 ;
					Cell_FaceNum[ C_Num ]	=	4 ;

					name 					=	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					Cell_Typename[ C_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;
								
						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					Cell_Type[ C_Num ]	=	Typename_Type[ name ] ;					

					for ( j = 0 ; j < 4 ; j++ )
						_Node[ j ] =	ielem_QUAD_4[ i ][ j ] - 1 ;

					CellNodeCCW( 4, Node_x, Node_y, _Node.get() ) ;

					Cell_x[ C_Num ]	=	0. ;
					Cell_y[ C_Num ]	=	0. ;
					Cell_z[ C_Num ]	=	0. ;
					for ( j = 0 ; j < 4 ; j++ )
					{
						NodeNo =	_Node[ j ] ;
						Cell_Node[ C_Num ].push_back( NodeNo ) ;
						Node_Cell[ NodeNo ].push_back( C_Num ) ;

						Cell_x[ C_Num ] +=	Node_x[ NodeNo ] ;
						Cell_y[ C_Num ] +=	Node_y[ NodeNo ] ;		
					}
					Cell_x[ C_Num ]	/= 4 ;
					Cell_y[ C_Num ]	/= 4 ;

					tmp[ 1 ]++ ;
					C_Num++ ;
				}
				Num[ 5 ]	+=	tmp[ 0 ] ;
			}
		}
	} else if ( dimension == 3 )
	{
		for ( index_sect = 1 ; index_sect <= nsections ; index_sect++ )
		{
			cg_section_read( index_file, 1, 1, index_sect, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag ) ;

			tmp[ 0 ]	=	( int ) iend - ( int ) istart + 1 ;
			tmp[ 1 ]	=	( int ) istart ;

			if ( itype == TRI_3 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_TRI_3[ Num[ 4 ] ], &iparentdata ) ;

				for ( i = Num[ 4 ] ; i < Num[ 4 ] + tmp[ 0 ] ; i++ )
				{
					name =	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					BCFace_Typename[ F_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;

						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					BCFace_Type[ F_Num ] =	Typename_Type[ name ] ;

					for ( j = 0 ; j < 3 ; j++ )
					{
						NodeNo =	ielem_TRI_3[ i ][ j ] - 1 ;
						BCFace_Node[ F_Num ].push_back( NodeNo ) ;
						Node_BCFace[ NodeNo ].push_back( F_Num ) ;
					}
					tmp[ 1 ]++ ;
					F_Num++ ;
				}
				Num[ 4 ]	+=	tmp[ 0 ] ;
			} else if ( itype == QUAD_4 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_QUAD_4[ Num[ 5 ] ], &iparentdata ) ;

				for ( i = Num[ 5 ] ; i < Num[ 5 ] + tmp[ 0 ] ; i++ )
				{
					name =	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					BCFace_Typename[ F_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;
						
						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					BCFace_Type[ F_Num ] =	Typename_Type[ name ] ;

					for ( j = 0 ; j < 4 ; j++ )
					{
						NodeNo =	ielem_QUAD_4[ i ][ j ] - 1 ;
						BCFace_Node[ F_Num ].push_back( NodeNo ) ;
						Node_BCFace[ NodeNo ].push_back( F_Num ) ;
					}
					tmp[ 1 ]++ ;
					F_Num++ ;
				}
				Num[ 5 ]	+=	tmp[ 0 ] ;							
			}
		}

		if ( check_Typename[ "Interface" ] == 0 )
		{
			BC_Typename->push_back( "Interface" ) ;

			Typename_Type[ "Interface" ]	=	(*BC_TypeNum) ;
			check_Typename[ "Interface" ]	=	1 ;
			(*BC_TypeNum)++ ;
		}
		(*BC_TypenNum_forFace)	=	(*BC_TypeNum) ;

		for ( index_sect = 1 ; index_sect <= nsections ; index_sect++ )
		{
			cg_section_read( index_file, 1, 1, index_sect, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag ) ;

			tmp[ 0 ]	=	( int ) iend - ( int ) istart + 1 ;
			tmp[ 1 ]	=	( int ) istart ;

			if ( itype == TETRA_4 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_TETRA_4[ Num[ 0 ] ], &iparentdata ) ;

				for ( i = Num[ 0 ] ; i < Num[ 0 ] + tmp[ 0 ] ; i++ )
				{
					Cell_Form[ C_Num ]		=	0 ;
					Cell_FaceNum[ C_Num ]	=	4 ;

					name =	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					Cell_Typename[ C_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;

						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					Cell_Type[ C_Num ]	=	Typename_Type[ name ] ;					

					Cell_x[ C_Num ]	=	0. ;
					Cell_y[ C_Num ]	=	0. ;
					Cell_z[ C_Num ]	=	0. ;
					for ( j = 0 ; j < 4 ; j++ )
					{
						NodeNo =	ielem_TETRA_4[ i ][ j ] - 1 ;
						Cell_Node[ C_Num ].push_back( NodeNo ) ;
						Node_Cell[ NodeNo ].push_back( C_Num ) ;

						Cell_x[ C_Num ] +=	Node_x[ NodeNo ] ;
						Cell_y[ C_Num ] +=	Node_y[ NodeNo ] ;
						Cell_z[ C_Num ] +=	Node_z[ NodeNo ] ;
					}
					Cell_x[ C_Num ]	/= 4 ;
					Cell_y[ C_Num ]	/= 4 ;
					Cell_z[ C_Num ]	/= 4 ;

					tmp[ 1 ]++ ;
					C_Num++ ;
				}
				Num[ 0 ]	+=	tmp[ 0 ] ;
			} else if ( itype == PYRA_5 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_PYRA_5[ Num[ 3 ] ], &iparentdata ) ;

				for ( i = Num[ 3 ] ; i < Num[ 3 ] + tmp[ 0 ] ; i++ )
				{
					Cell_Form[ C_Num ]		=	3 ;
					Cell_FaceNum[ C_Num ]	=	5 ;

					name =	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					Cell_Typename[ C_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;
								
						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					Cell_Type[ C_Num ]	=	Typename_Type[ name ] ;					

					Cell_x[ C_Num ]	=	0. ;
					Cell_y[ C_Num ]	=	0. ;
					Cell_z[ C_Num ]	=	0. ;
					for ( j = 0 ; j < 5 ; j++ )
					{
						NodeNo =	ielem_PYRA_5[ i ][ j ] - 1 ;
						Cell_Node[ C_Num ].push_back( NodeNo ) ;
						Node_Cell[ NodeNo ].push_back( C_Num ) ;

						Cell_x[ C_Num ] +=	Node_x[ NodeNo ] ;
						Cell_y[ C_Num ] +=	Node_y[ NodeNo ] ;
						Cell_z[ C_Num ] +=	Node_z[ NodeNo ] ;
					}
					Cell_x[ C_Num ]	/= 5 ;
					Cell_y[ C_Num ]	/= 5 ;
					Cell_z[ C_Num ]	/= 5 ;

					tmp[ 1 ]++ ;
					C_Num++ ;
				}
				Num[ 3 ]	+=	tmp[ 0 ] ;
			} else if ( itype == PENTA_6 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_PENTA_6[ Num[ 2 ] ], &iparentdata ) ;

				for ( i = Num[ 2 ] ; i < Num[ 2 ] + tmp[ 0 ] ; i++ )
				{
					Cell_Form[ C_Num ]		=	2 ;
					Cell_FaceNum[ C_Num ]	=	5 ;

					name =	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					Cell_Typename[ C_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;

						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					Cell_Type[ C_Num ]	=	Typename_Type[ name ] ;

					Cell_x[ C_Num ]	=	0. ;
					Cell_y[ C_Num ]	=	0. ;
					Cell_z[ C_Num ]	=	0. ;
					for ( j = 0 ; j < 6 ; j++ )
					{
						NodeNo =	ielem_PENTA_6[ i ][ j ] - 1 ;
						Cell_Node[ C_Num ].push_back( NodeNo ) ;
						Node_Cell[ NodeNo ].push_back( C_Num ) ;

						Cell_x[ C_Num ] +=	Node_x[ NodeNo ] ;
						Cell_y[ C_Num ] +=	Node_y[ NodeNo ] ;
						Cell_z[ C_Num ] +=	Node_z[ NodeNo ] ;
					}
					Cell_x[ C_Num ]	/= 6 ;
					Cell_y[ C_Num ]	/= 6 ;
					Cell_z[ C_Num ]	/= 6 ;

					tmp[ 1 ]++ ;
					C_Num++ ;
				}
				Num[ 2 ]	+=	tmp[ 0 ] ;
			} else if ( itype == HEXA_8 )
			{
				cg_elements_read( index_file, 1, 1, index_sect, ielem_HEXA_8[ Num[ 1 ] ], &iparentdata ) ;

				for ( i = Num[ 1 ] ; i < Num[ 1 ] + tmp[ 0 ] ; i++ )
				{					
					Cell_Form[ C_Num ]		=	1 ;
					Cell_FaceNum[ C_Num ]	=	6 ;

					name =	Typename[ BCNo[ tmp[ 1 ] ] ] ;
					Cell_Typename[ C_Num ]	=	name ;

					if ( check_Typename[ name ] == 0 )
					{
						BC_Typename->push_back( name ) ;

						Typename_Type[ name ]	=	(*BC_TypeNum) ;
						check_Typename[ name ]	=	1 ;
						(*BC_TypeNum)++ ;
					}
					Cell_Type[ C_Num ]	=	Typename_Type[ name ] ;

					Cell_x[ C_Num ]	=	0. ;
					Cell_y[ C_Num ]	=	0. ;
					Cell_z[ C_Num ]	=	0. ;
					for ( j = 0 ; j < 8 ; j++ )
					{
						NodeNo =	ielem_HEXA_8[ i ][ j ] - 1 ;
						Cell_Node[ C_Num ].push_back( NodeNo ) ;
						Node_Cell[ NodeNo ].push_back( C_Num ) ;

						Cell_x[ C_Num ] +=	Node_x[ NodeNo ] ;
						Cell_y[ C_Num ] +=	Node_y[ NodeNo ] ;
						Cell_z[ C_Num ] +=	Node_z[ NodeNo ] ;						
					}
					Cell_x[ C_Num ]	/= 8 ;
					Cell_y[ C_Num ]	/= 8 ;
					Cell_z[ C_Num ]	/= 8 ;

					tmp[ 1 ]++ ;	
					C_Num++	;
				}
				Num[ 1 ]	+=	tmp[ 0 ] ;
			}
		}
	}

	// close CGNS file
	cg_close( index_file ) ;
}

#include <mpi.h>
#include <metis.h>
#include <iostream>
#include <cgnslib.h>
#include "domain.h"
#include "preprocessingtool.h"
#include "preprocessing.h"
#include "fieldview.h"
#include "cgns.h"
#include "gmsh.h"
#include "pti.h"

using namespace std ;

void Domain::Preprocessing (  )
{
	int	*NULLINT ;

	NULLINT	=	NULL ;
	Preprocessing ( NULLINT ) ;
}

void Domain::Preprocessing ( int *CellWgt )
{
	//cout<<"GetMeshInformation( ) ;"<<endl;
	GetMeshInformation( ) ;
	//cout<<"&GetMeshInformation( ) ;"<<endl<<endl;

	//cout<<"DomainPartition( ) ;"<<endl;
	DomainPartition( CellWgt ) ;
	//cout<<"&DomainPartition( ) ;"<<endl<<endl;

	//cout<<"CreateCellInformation( ) ;"<<endl;
	CreateCellInformation( ) ;
	//cout<<"&CreateCellInformation( ) ;"<<endl<<endl;

	//cout<<"CreateFaceInformation( ) ;"<<endl;
	CreateFaceInformation( ) ;
	//cout<<"&CreateFaceInformation( ) ;"<<endl<<endl;
}

void Domain::GetMeshInformation( )
{
	// 0 = FieldView; 1 = CGNS; 2 = MSH
	//cout<<"Domain::GetMeshInformation( )"<<endl;
	if ( flag_structured_mesh )
	{	
		cout<<"Domain::GetMeshInformation( ), flag_structured_mesh, block by Kuan-Lin"<<endl;exit(1) ;
		// Mesh.NodeNum =	Nx * Ny * Nz ;

		// if ( dimension == 2 )
		// {
		// 	Mesh.CellNum 	=	( Nx -1 ) * ( Ny - 1 ) ;
		// 	Mesh.BCFaceNum 	=	2 * ( Nx - 1 ) + 2 * ( Ny - 1 ) ;
		// } else
		// {
		// 	Mesh.CellNum 	=	( Nx -1 ) * ( Ny - 1 ) * ( Nz - 1 ) ;
		// 	Mesh.BCFaceNum 	=	2 * ( ( Nx - 1 ) * ( Ny - 1 ) + ( Ny - 1 ) * ( Nz - 1 ) + ( Nz -1 ) * ( Nx - 1 ) ) ;
		// }

		// Mesh.set_datasize_NodeNum( Mesh.NodeNum ) ;
		// Mesh.set_datasize_CellNum( Mesh.CellNum ) ;
		// Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;

		// Structured_mesh( dimension, mesh_scale, Nx, Ny, Nz, Dx, Dy, Dz, Mesh.NodeNum, Mesh.CellNum, Mesh.Node_Position[ 0 ].get(), Mesh.Node_Position[ 1 ].get(), Mesh.Node_Position[ 2 ].get(), Mesh.Cell_Position[ 0 ].get(), Mesh.Cell_Position[ 1 ].get(), Mesh.Cell_Position[ 2 ].get(), &Mesh.BC_TypeNum, &Mesh.BC_Typename, Mesh.BCFace_Type.get(), Mesh.BCFace_Typename.get(), Mesh.Cell_FaceNum.get(), Mesh.Cell_Node.get(), Mesh.Node_Cell.get(), Mesh.BCFace_Node.get(), Mesh.Node_BCFace.get() ) ;
 
		// Mesh.BC_TypeNum_forFace	=	Mesh.BC_TypeNum ;

 	// 	if ( dimension == 2 )
 	// 	{
		// 	for ( int i = 0 ; i < Mesh.CellNum ; i++ )
		// 	{
		// 		Mesh.Cell_Form[ i ]		=	5 ;
		// 		Mesh.Cell_Type[ i ]		=	0 ;
		// 		Mesh.Cell_Typename[ i ]	=	"Bulk" ;
		// 	}
		// } else
		// {
		// 	for ( int i = 0 ; i < Mesh.CellNum ; i++ )
		// 	{
		// 		Mesh.Cell_Form[ i ]		=	1 ;
		// 		Mesh.Cell_Type[ i ]		=	0 ;
		// 		Mesh.Cell_Typename[ i ]	=	"Bulk" ;				
		// 	}
		// }
	} else if ( Get_filetype( meshfile ) == 0 )	// format: uns
	{
		cout<<"Domain::GetMeshInformation( ), format: uns, block by Kuan-Lin"<<endl;exit(1) ;
		// // To get the array size from grids information.
		// GetNumbers_FieldView ( dimension, &Mesh.NodeNum, &Mesh.CellNum, &Mesh.BC_TypeNum, &Mesh.BCFaceNum, meshfile ) ;

		// Mesh.set_datasize_NodeNum( Mesh.NodeNum ) ;
		// Mesh.set_datasize_CellNum( Mesh.CellNum ) ;

		// boost::shared_array<int>			_BCFace_Type ;
		// boost::shared_array<string> 		_BCFace_Typename ;
		// boost::shared_array<vector<int>>	_BCFace_Node ;	

		// _BCFace_Type 		=	boost::shared_array<int> ( new int [ Mesh.BCFaceNum ] ) ;
		// _BCFace_Typename 	=	boost::shared_array<string> ( new string [ Mesh.BCFaceNum ] ) ;
		// _BCFace_Node 		=	boost::shared_array<vector<int>> ( new vector<int> [ Mesh.BCFaceNum ] ) ;

		// for ( int i = 0 ; i < Mesh.BCFaceNum ; i++ )
		// {
		// 	_BCFace_Type[ i ]		=	0 ;
		// 	_BCFace_Typename[ i ]	=	"Bulk" ;
		// }

		// // Read Girds Information, including Node, Cell and BCFace.
		// // Create Cell-Node, Noce-Cell, BCFace-Node & Node-BCFace relations
		// ReadGrid_FieldView ( dimension, mesh_scale, Mesh.NodeNum, Mesh.CellNum, Mesh.Node_Position[ 0 ].get(), Mesh.Node_Position[ 1 ].get(), Mesh.Node_Position[ 2 ].get(), Mesh.Cell_Position[ 0 ].get(), Mesh.Cell_Position[ 1 ].get(), Mesh.Cell_Position[ 2 ].get(), &Mesh.BC_TypeNum, &Mesh.BC_TypeNum_forFace, &Mesh.BC_Typename, &Mesh.BCFaceNum, _BCFace_Type.get(), _BCFace_Typename.get(), _BCFace_Node.get(), Mesh.Cell_Form.get(), Mesh.Cell_Type.get(), Mesh.Cell_Typename.get(), Mesh.Cell_FaceNum.get(), Mesh.Cell_Node.get(), Mesh.Node_Cell.get(), Mesh.Node_BCFace.get(), meshfile ) ;
		
		// Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;
		// for ( int i = 0 ; i < Mesh.BCFaceNum ; i++ )
		// {
		// 	Mesh.BCFace_Type[ i ]		=	_BCFace_Type[ i ] ;
		// 	Mesh.BCFace_Typename[ i ]	=	_BCFace_Typename[ i ] ;

		// 	for ( int j = 0 ; j < _BCFace_Node[ i ].size() ; j++ )
		// 		Mesh.BCFace_Node[ i ].push_back( _BCFace_Node[ i ][ j ] ) ;
		// }		
	} else if ( Get_filetype( meshfile ) == 1 )	// format: cgns
	{	
		cout<<"Domain::GetMeshInformation( ), format: cgns, block by Kuan-Lin"<<endl;exit(1) ;
		// boost::shared_array<int>	_MeshType_Num ;
		// _MeshType_Num 		=	boost::shared_array<int> ( new int [ 7 ] ) ;
		// for ( int i = 0 ; i < 7 ; i++ )
		// 	_MeshType_Num[ i ]	=	0 ;

		// GetNumbers_CGNS ( dimension, &Mesh.NodeNum, &Mesh.CellNum, &Mesh.BC_TypeNum, &Mesh.BCFaceNum, _MeshType_Num.get(), meshfile ) ;

		// Mesh.set_datasize_NodeNum( Mesh.NodeNum ) ;
		// Mesh.set_datasize_CellNum( Mesh.CellNum ) ;
		// Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;

		// ReadGrid_CGNS ( dimension, mesh_scale, Mesh.NodeNum, Mesh.CellNum, _MeshType_Num.get(), Mesh.Node_Position[ 0 ].get(), Mesh.Node_Position[ 1 ].get(), Mesh.Node_Position[ 2 ].get(), Mesh.Cell_Position[ 0 ].get(), Mesh.Cell_Position[ 1 ].get(), Mesh.Cell_Position[ 2 ].get(), &Mesh.BC_TypeNum, &Mesh.BC_TypeNum_forFace, &Mesh.BC_Typename, Mesh.BCFaceNum, Mesh.BCFace_Type.get(), Mesh.BCFace_Typename.get(), Mesh.BCFace_Node.get(), Mesh.Cell_Form.get(), Mesh.Cell_Type.get(), Mesh.Cell_Typename.get(), Mesh.Cell_FaceNum.get(), Mesh.Cell_Node.get(), Mesh.Node_Cell.get(), Mesh.Node_BCFace.get(), meshfile ) ;

	} else if ( Get_filetype( meshfile ) == 2 )	//	format: msh
	{
		// To get the array size from grids information.
		//cout<<"Read Gemsh NodeNum, CellNum, BCFaceNum..."<<endl;
		GetNumbers_GMSH ( dimension, &Mesh.NodeNum, &Mesh.CellNum, &Mesh.BC_TypeNum, &Mesh.BCFaceNum, meshfile ) ;
		//cout<<"&Read Gemsh NodeNum, CellNum, BCFaceNum..."<<endl<<endl;
		
		//cout<<"Mesh.set_datasize_NodeNum( Mesh.NodeNum )"<<endl;
		Mesh.set_datasize_NodeNum( Mesh.NodeNum ) ;
		//cout<<"&Mesh.set_datasize_NodeNum( Mesh.NodeNum )"<<endl<<endl;

		//cout<<"Mesh.set_datasize_CellNum( Mesh.CellNum ) ;"<<endl;
		Mesh.set_datasize_CellNum( Mesh.CellNum ) ;
		//cout<<"&Mesh.set_datasize_CellNum( Mesh.CellNum ) ;"<<endl<<endl;

		//cout<<"Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;"	<<endl;
		Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;
		//cout<<"&Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;"<<endl<<endl;

		// Read Girds Information, including Node, Cell and BCFace.
		// Create Cell-Node, Noce-Cell, BCFace-Node & Node-BCFace relations
		//cout<<"ReadGrid_GMSH"<<endl;
		ReadGrid_GMSH ( dimension, mesh_scale, Mesh.Node_Position[ 0 ].get(), Mesh.Node_Position[ 1 ].get(), Mesh.Node_Position[ 2 ].get(), Mesh.Cell_Position[ 0 ].get(), Mesh.Cell_Position[ 1 ].get(), Mesh.Cell_Position[ 2 ].get(), &Mesh.BC_TypeNum, &Mesh.BC_TypeNum_forFace, &Mesh.BC_Typename, Mesh.Cell_Form.get(), Mesh.Cell_FaceNum.get(), Mesh.Cell_Type.get(), Mesh.Cell_Typename.get(), Mesh.Cell_Node.get(), Mesh.Node_Cell.get(), Mesh.BCFace_Type.get(), Mesh.BCFace_Typename.get(), Mesh.BCFace_Node.get(), Mesh.Node_BCFace.get(), meshfile ) ;
		//cout<<"&ReadGrid_GMSH"<<endl<<endl;

	} else if ( Get_filetype( meshfile ) == 3 )	//	format: pti
	{
		cout<<"Domain::GetMeshInformation( ), format: pti, block by Kuan-Lin"<<endl;exit(1) ;
		// boost::shared_array<int>	_MeshType_Num ;
		// boost::shared_array<int>	_MeshType_Num ;
		// _MeshType_Num 		=	boost::shared_array<int> ( new int [ 7 ] ) ;
		// for ( int i = 0 ; i < 7 ; i++ )
		// 	_MeshType_Num[ i ]	=	0 ;

		// GetNumbers_PTI ( dimension, &Mesh.NodeNum, &Mesh.CellNum, &Mesh.BC_TypeNum, &Mesh.BCFaceNum, _MeshType_Num.get(), meshfile ) ;

		// Mesh.set_datasize_NodeNum( Mesh.NodeNum ) ;
		// Mesh.set_datasize_CellNum( Mesh.CellNum ) ;
		// Mesh.set_datasize_BCFaceNum( Mesh.BCFaceNum ) ;

		// ReadGrid_PTI ( dimension, mesh_scale, Mesh.NodeNum, _MeshType_Num.get(), Mesh.Node_Position[ 0 ].get(), Mesh.Node_Position[ 1 ].get(), Mesh.Node_Position[ 2 ].get(), Mesh.Cell_Position[ 0 ].get(), Mesh.Cell_Position[ 1 ].get(), Mesh.Cell_Position[ 2 ].get(), &Mesh.BC_TypeNum, &Mesh.BC_TypeNum_forFace, &Mesh.BC_Typename, Mesh.Cell_Form.get(), Mesh.Cell_FaceNum.get(), Mesh.Cell_Type.get(), Mesh.Cell_Typename.get(), Mesh.Cell_Node.get(), Mesh.Node_Cell.get(), Mesh.BCFace_Type.get(), Mesh.BCFace_Typename.get(), Mesh.BCFace_Node.get(), Mesh.Node_BCFace.get(), meshfile ) ;
	}

	/*---	Create Cell-Cell relation, and return inner FaceNum  ---*/
	//CreateCellCellRelation( dimension, Mesh.NodeNum, Mesh.CellNum, &Mesh.InnerFaceNum, &Mesh.InterFaceNum, Mesh.Node_Cell.get(), Mesh.Cell_Type.get(), Mesh.Cell_Node.get(), Mesh.Cell_Cell.get() ) ;

	//cout<<"CreateCellCellRelation"<<endl;
	CreateCellCellRelation( Mesh.NodeNum, Mesh.CellNum, &Mesh.InnerFaceNum, &Mesh.InterFaceNum, &Mesh.FaceNum, Mesh.Node_Cell.get(), Mesh.Cell_Form.get(), Mesh.Cell_Type.get(), Mesh.Cell_Node.get(), Mesh.Cell_Cell.get() ) ;
	//cout<<"&CreateCellCellRelation"<<endl<<endl;
	//cout<<"set_datasize_FaceNum"<<endl;
	Mesh.set_datasize_FaceNum( Mesh.FaceNum ) ;
	//cout<<"&set_datasize_FaceNum"<<endl<<endl;

	//cout<<"&Domain::GetMeshInformation( )"<<endl;
}

void Domain::DomainPartition( int *CellWgt )
{
	set_table_processor( comm_size ) ;

	if ( comm_size != 1 )
	{
		int		MetisState ;
		idx_t	*MetisCellProcessorNo, *MetisXadj, *MetisAdjncy, *MetisVwgt ;
		idx_t	MetisCellNum, MetisProcessorNum, MetisObjval, MetisNcon ;

		MetisState				=	0 ;
		MetisCellNum			=	Mesh.CellNum ;
		MetisProcessorNum		=	comm_size ;
		MetisObjval				=	0 ;
		MetisNcon				=	1 ;

		MetisCellProcessorNo	=	new idx_t[ Mesh.CellNum ] ;
		MetisXadj				=	new idx_t[ Mesh.CellNum + 1 ] ;
		MetisAdjncy				=	new idx_t[ Mesh.InnerFaceNum * 2 ] ;
		MetisVwgt				=	new idx_t[ Mesh.CellNum ] ;

		// Create CSR format for graph data.
		MetisXadj[ 0 ]	= 0 ;
		for ( int i = 0 ; i < Mesh.CellNum ; i++ )
		{
			MetisXadj[ i + 1 ]	=	MetisXadj[ i ] ;

			for ( int j = 0 ; j < Mesh.Cell_Cell[ i ].size() ; j++ )
			{
				MetisAdjncy[ MetisXadj[ i + 1 ] ]	=	Mesh.Cell_Cell[ i ][ j ] ;
				MetisXadj[ i + 1 ]++ ;
			}
		}
		
		// Meitis API for graph partition.
		MetisState	=	METIS_PartGraphRecursive( &MetisCellNum, &MetisNcon, MetisXadj, MetisAdjncy, CellWgt, NULL, NULL, &MetisProcessorNum, NULL, NULL, NULL, &MetisObjval, MetisCellProcessorNo ) ;

		for ( int i = 0 ; i < Mesh.CellNum ; i++ )
		{
			Mesh.Cell_ProcessorNo[ i ]	=	MetisCellProcessorNo[ i ] ;
			Processor_MeshCell[ MetisCellProcessorNo[ i ] ].push_back( i ) ;
		}

		delete [] MetisCellProcessorNo ;
		delete [] MetisXadj ;
		delete [] MetisAdjncy ;
		delete [] MetisVwgt ;
	} else
	{
		for ( int i = 0 ; i < Mesh.CellNum ; i++ )
		{
			Mesh.Cell_ProcessorNo[ i ]	=	0 ;
			Processor_MeshCell[ 0 ].push_back( i ) ;
		}
	}	
}

void Domain::CreateCellInformation( )
{
	int MeshNo, No[ 2 ], ProcNo, GlobalNo = 0, LocalNo = 0, NodeNo ;
	map<int, int> C, N ;
	CellMapping	pMapping ;

	set_table_globalcell( Mesh.CellNum ) ;

	for ( int i = 0 ; i < comm_size ; i++ )
	{
		for ( int j = 0 ; j < Processor_MeshCell[ i ].size() ; j++ )
		{
			MeshNo =	Processor_MeshCell[ i ][ j ] ;

			GlobalCell_MeshCellNo[ GlobalNo ]	=	MeshNo ;
			MeshCell_GlobalCellNo[ MeshNo ]		=	GlobalNo ;

			if ( i == comm_rank )
			{
				for ( int k = 0 ; k < Mesh.Cell_Cell[ MeshNo ].size() ; k++ )
				{
					No[ 0 ] =	Mesh.Cell_Cell[ MeshNo ][ k ] ;
					ProcNo 	=	Mesh.Cell_ProcessorNo[ No[ 0 ] ] ;

					if ( ( ProcNo != i ) && ( C[ No[ 0 ] ] == 0 ) )
					{
						Processor_GhostMeshCell_lv1[ i ].push_back( No[ 0 ] ) ;
					}
					C[ No[ 0 ] ]++ ;
				}
			}
			C[ MeshNo ]++ ;
			GlobalNo++ ;
		}

		if ( i == comm_rank )
		{
			for ( int j = 0 ; j < Processor_GhostMeshCell_lv1[ i ].size() ; j++ )
			{
				No[ 0 ] =	Processor_GhostMeshCell_lv1[ i ][ j ] ;

				for ( int k = 0 ; k < Mesh.Cell_Cell[ No[ 0 ] ].size() ; k++ )
				{
					No[ 1 ]	=	Mesh.Cell_Cell[ No[ 0 ] ][ k ] ;

					if ( C[ No[ 1 ] ] == 0 )
					{
						Processor_GhostMeshCell_lv2[ i ].push_back( No[ 1 ] ) ;
						C[ No[ 1 ] ]++ ;
					}
				}
			}

			for ( int j = 0 ; j < Processor_MeshCell[ i ].size() ; j++ )
			{
				MeshNo =	Processor_MeshCell[ i ][ j ] ;

				for ( int k = 0 ; k < Mesh.Cell_Node[ MeshNo ].size() ; k++ )
				{
					No[ 0 ]	=	Mesh.Cell_Node[ MeshNo ][ k ] ;

					if ( N[ No[ 0 ] ] == 0 )
					{
						for ( int m = 0 ; m < Mesh.Node_Cell[ No[ 0 ] ].size() ; m++ )
						{
							No[ 1 ]	=	Mesh.Node_Cell[ No[ 0 ] ][ m ] ;

							if ( C[ No[ 1 ] ] == 0 )
							{
								Processor_GhostMeshCell_lv3[ i ].push_back( No[ 1 ] ) ;
								C[ No[ 1 ] ]++ ;
							}
						}
						N[ No[ 0 ] ]++ ;
					}
				}
			}
		}
		C.clear() ;
		N.clear() ;
	}

	set_table_localcell( Processor_MeshCell[ comm_rank ].size() + Processor_GhostMeshCell_lv1[ comm_rank ].size() + Processor_GhostMeshCell_lv2[ comm_rank ].size() + Processor_GhostMeshCell_lv3[ comm_rank ].size() ) ;

	for ( int i = 0 ; i < Processor_MeshCell[ comm_rank ].size() ; i++ )
	{
		MeshNo 		=	Processor_MeshCell[ comm_rank ][ i ] ;
		GlobalNo 	=	MeshCell_GlobalCellNo[ MeshNo ] ;

		LocalCell_GlobalCellNo[ i ]			=	GlobalNo ;
		GlobalCell_LocalCellNo[ GlobalNo ]	=	i ;

		LocalCell_MeshCellNo[ i ]			=	MeshNo ;
		MeshCell_LocalCellNo[ MeshNo ]		=	i ;

		for ( int j = 0 ; j < Mesh.Cell_Node[ MeshNo ].size() ; j++ )
		{
			NodeNo =	Mesh.Cell_Node[ MeshNo ][ j ] ;

			if ( N[ NodeNo ] == 0 )
			{
				Processor_MeshNode[ comm_rank ].push_back( NodeNo ) ;
			}
			N[ NodeNo ]++ ;
		}
	}

	for ( int i = 0 ; i < Processor_GhostMeshCell_lv1[ comm_rank ].size() ; i++ )
	{
		No[ 0 ]	=	Processor_MeshCell[ comm_rank ].size() + i ;

		MeshNo 		=	Processor_GhostMeshCell_lv1[ comm_rank ][ i ] ;
		GlobalNo 	=	MeshCell_GlobalCellNo[ MeshNo ] ;

		LocalCell_GlobalCellNo[ No[ 0 ] ]	=	GlobalNo ;
		GlobalCell_LocalCellNo[ GlobalNo ]	=	No[ 0 ] ;

		LocalCell_MeshCellNo[ No[ 0 ] ]		=	MeshNo ;
		MeshCell_LocalCellNo[ MeshNo ]		=	No[ 0 ] ;

		for ( int j = 0 ; j < Mesh.Cell_Node[ MeshNo ].size() ; j++ )
		{
			NodeNo =	Mesh.Cell_Node[ MeshNo ][ j ] ;
			if ( N[ NodeNo ] == 0 )
			{
				Processor_GhostMeshNode_lv1[ comm_rank ].push_back( NodeNo ) ;
			}
			N[ NodeNo ]++ ;
		}
	}

	for ( int i = 0 ; i < Processor_GhostMeshCell_lv2[ comm_rank ].size() ; i++ )
	{
		No[ 0 ]	=	Processor_MeshCell[ comm_rank ].size() + Processor_GhostMeshCell_lv1[ comm_rank ].size() + i ;

		MeshNo 	=	Processor_GhostMeshCell_lv2[ comm_rank ][ i ] ;
		GlobalNo =	MeshCell_GlobalCellNo[ MeshNo ] ;

		LocalCell_GlobalCellNo[ No[ 0 ] ]	=	GlobalNo ;
		GlobalCell_LocalCellNo[ GlobalNo ]	=	No[ 0 ] ;

		LocalCell_MeshCellNo[ No[ 0 ] ]		=	MeshNo ;
		MeshCell_LocalCellNo[ MeshNo ]		=	No[ 0 ] ;
	}

	for ( int i = 0 ; i < Processor_GhostMeshCell_lv3[ comm_rank ].size() ; i++ )
	{
		No[ 0 ]	=	Processor_MeshCell[ comm_rank ].size() + Processor_GhostMeshCell_lv1[ comm_rank ].size() + Processor_GhostMeshCell_lv2[ comm_rank ].size() + i ;

		MeshNo 	=	Processor_GhostMeshCell_lv3[ comm_rank ][ i ] ;
		GlobalNo =	MeshCell_GlobalCellNo[ MeshNo ] ;

		LocalCell_GlobalCellNo[ No[ 0 ] ]	=	GlobalNo ;
		GlobalCell_LocalCellNo[ GlobalNo ]	=	No[ 0 ] ;

		LocalCell_MeshCellNo[ No[ 0 ] ]		=	MeshNo ;
		MeshCell_LocalCellNo[ MeshNo ]		=	No[ 0 ] ;
	}

	set_table_node( Mesh.NodeNum, Processor_MeshNode[ comm_rank ].size() + Processor_GhostMeshNode_lv1[ comm_rank ].size() ) ;

	for ( int i = 0 ; i < Processor_MeshNode[ comm_rank ].size() ; i++ )
	{
		NodeNo =	Processor_MeshNode[ comm_rank ][ i ] ;

		MeshNode_LocalNodeNo[ NodeNo ]	=	i ;
		LocalNode_MeshNodeNo[ i ]		=	NodeNo ;
	}

	for ( int i = 0 ; i < Processor_GhostMeshNode_lv1[ comm_rank ].size() ; i++ )
	{
		No[ 0 ] =	Processor_MeshNode[ comm_rank ].size() + i ;
		NodeNo 	=	Processor_GhostMeshNode_lv1[ comm_rank ][ i ] ;

		MeshNode_LocalNodeNo[ NodeNo ]	=	No[ 0 ] ;
		LocalNode_MeshNodeNo[ No[ 0 ] ]	=	NodeNo ;
	}

	CalculateCellVolume( geometry, Mesh.CellNum, Mesh.Cell_Form.get(), &min_cell_length, Mesh.Node_Position[ 0 ].get(), Mesh.Node_Position[ 1 ].get(), Mesh.Node_Position[ 2 ].get(), Mesh.Cell_Position[ 0 ].get(), Mesh.Cell_Position[ 1 ].get(), Mesh.Cell_Node.get(), Mesh.Cell_Volume.get(), &pMapping ) ;
}

void Domain::CreateFaceInformation( )
{
	int c1, c2, Type1, Type2, j2, NodeNo, BCFaceNo, CellNo, FaceNo, No ;

	/*--- Modifyed by Kuan-Lin ---*/
	//int Check[ Mesh.CellNum ][ 6 ] ;
	//int Check_BCFace[ Mesh.BCFaceNum ] ;
	int **Check = new int *[ Mesh.CellNum ] ;
	for (int i=0; i < Mesh.CellNum ; i++ ) Check[ i ] = new int [6] ;
	int *Check_BCFace = new int [Mesh.BCFaceNum] ;

	/*----------------------------*/
	map<int, int>	F ;
	CellMapping	pMapping ;

	for ( int i = 0 ; i < Mesh.CellNum ; i++ )
	{
		for ( int j = 0 ; j < 6 ; j++ )
		{
			Check[ i ][ j ] = -999 ;
		}
	}

	for ( int i = 0 ; i < Mesh.BCFaceNum ; i++ )
	{
		Check_BCFace[ i ]	=	-999 ;
	}

	// InnerFace
	FaceNo = 0 ;
	for ( int i = 0 ; i < Mesh.NodeNum ; i++ )
	{
		for ( int j1 = 0 ; j1 < Mesh.Node_Cell[ i ].size() ; j1++ )
		{
			c1	=	Mesh.Node_Cell[ i ][ j1 ] ;

			Type1	=	Mesh.Cell_Form[ c1 ] ;
			for ( int k1 = 0 ; k1 < Mesh.Cell_FaceNum[ c1 ] ; k1++ )
			{
				if ( Check[ c1 ][ k1 ] == -999 )
				{
					j2	= 0 ;
					while ( j2 < Mesh.Node_Cell[ i ].size() )
					{
						if ( j2 != j1  )
						{
							c2	=	Mesh.Node_Cell[ i ][ j2 ] ;

							Type2	=	Mesh.Cell_Form[ c2 ] ;
							for ( int k2 = 0 ; k2 < Mesh.Cell_FaceNum[ c2 ] ; k2++ )
							{							
								if ( CheckNeighbor( Type1, Type2, &Mesh.Cell_Node[ c1 ], &Mesh.Cell_Node[ c2 ], k1, k2, &Mesh.Face_Node[ FaceNo ], &pMapping ) )
								{								
									if ( Mesh.Cell_Type[ c1 ] != Mesh.Cell_Type[ c2 ] )
									{
										for ( int n = 0 ; n < Mesh.Node_BCFace[ i ].size() ; n++ )
										{
											BCFaceNo =	Mesh.Node_BCFace[ i ][ n ] ;

											if ( Check_BCFace[ BCFaceNo ] == -999 )
											{
												if ( CheckNeighbor( Type1, &Mesh.Cell_Node[ c1 ], k1, &Mesh.BCFace_Node[ BCFaceNo ], &pMapping ) )
												{
													Mesh.Face_Type[ FaceNo ]		=	Mesh.BCFace_Type[ BCFaceNo ] ;
													Mesh.Face_Typename[ FaceNo ]	=	Mesh.BCFace_Typename[ BCFaceNo ] ;

													Check_BCFace[ BCFaceNo ]		=	0 ;
													break ;
												}
											}
										}

										// without interface condition, create new one.
										if ( Mesh.Face_Type[ FaceNo ] == 0 )
										{
											Mesh.Face_Type[ FaceNo ]		=	Mesh.BC_TypeNum_forFace - 1 ;
											Mesh.Face_Typename[ FaceNo ]	=	"Interface" ;
										}
									} else
									{
										Mesh.Face_Type[ FaceNo ]		=	0 ;
										Mesh.Face_Typename[ FaceNo ]	=	"Bulk" ;
									}

									Mesh.Face_Cell[ FaceNo ].push_back( c1 ) ;
									Mesh.Face_Cell[ FaceNo ].push_back( c2 ) ;

									Mesh.Cell_Face[ c1 ].push_back( FaceNo ) ;
									Mesh.Cell_Face[ c2 ].push_back( FaceNo ) ;

									Mesh.Face_Position[ 0 ][ FaceNo ] = 0 ;
									Mesh.Face_Position[ 1 ][ FaceNo ] = 0 ;
									Mesh.Face_Position[ 2 ][ FaceNo ] = 0 ;																		
									for ( int n = 0 ; n < Mesh.Face_Node[ FaceNo ].size() ; n++ )
									{
										NodeNo = Mesh.Face_Node[ FaceNo ][ n ] ;
										Mesh.Node_Face[ NodeNo ].push_back( FaceNo ) ;

										Mesh.Face_Position[ 0 ][ FaceNo ]	+=	Mesh.Node_Position[ 0 ][ NodeNo ] ;
										Mesh.Face_Position[ 1 ][ FaceNo ]	+=	Mesh.Node_Position[ 1 ][ NodeNo ] ;
										Mesh.Face_Position[ 2 ][ FaceNo ]	+=	Mesh.Node_Position[ 2 ][ NodeNo ] ;										
									}
									Mesh.Face_Position[ 0 ][ FaceNo ] /= Mesh.Face_Node[ FaceNo ].size() ;
									Mesh.Face_Position[ 1 ][ FaceNo ] /= Mesh.Face_Node[ FaceNo ].size() ;
									Mesh.Face_Position[ 2 ][ FaceNo ] /= Mesh.Face_Node[ FaceNo ].size() ;

									Check[ c1 ][ k1 ]	=	0 ;
									Check[ c2 ][ k2 ]	=	0 ;

									FaceNo++ ;
									j2	=	Mesh.Node_Cell[ i ].size() ;
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

	// BCFace
	for ( int i = 0 ; i < Mesh.CellNum ; i++ )
	{		
		Type1	=	Mesh.Cell_Form[ i ] ;

		for ( int j = 0 ; j < Mesh.Cell_FaceNum[ i ] ; j++ )
		{			
			if ( Check[ i ][ j ] == -999 )
			{			
				NodeNo	=	pMapping.Node[ Type1 ][ j ][ 0 ] ;
				NodeNo	=	Mesh.Cell_Node[ i ][ NodeNo ] ;

				for ( int k = 0 ; k < Mesh.Node_BCFace[ NodeNo ].size() ; k++ )
				{
					BCFaceNo	=	Mesh.Node_BCFace[ NodeNo ][ k ] ;

					if ( Check_BCFace[ BCFaceNo ] == -999 )
					{
						if ( CheckNeighbor( Type1, &Mesh.Cell_Node[ i ], j, &Mesh.BCFace_Node[ BCFaceNo ], &Mesh.Face_Node[ FaceNo ], &pMapping ) )
						{
							Mesh.Face_Type[ FaceNo ]		=	Mesh.BCFace_Type[ BCFaceNo ] ;
							Mesh.Face_Typename[ FaceNo ]	=	Mesh.BCFace_Typename[ BCFaceNo ] ;

							Mesh.Face_Cell[ FaceNo ].push_back( i ) ;
							Mesh.Cell_Face[ i ].push_back( FaceNo ) ;

							Mesh.Face_Position[ 0 ][ FaceNo ] = 0 ;
							Mesh.Face_Position[ 1 ][ FaceNo ] = 0 ;
							Mesh.Face_Position[ 2 ][ FaceNo ] = 0 ;	
							for ( int n = 0 ; n < Mesh.Face_Node[ FaceNo ].size() ; n++ )
							{
								No =	Mesh.Face_Node[ FaceNo ][ n ] ;
								Mesh.Node_Face[ No ].push_back( FaceNo ) ;

								Mesh.Face_Position[ 0 ][ FaceNo ]	+=	Mesh.Node_Position[ 0 ][ No ] ;
								Mesh.Face_Position[ 1 ][ FaceNo ]	+=	Mesh.Node_Position[ 1 ][ No ] ;
								Mesh.Face_Position[ 2 ][ FaceNo ]	+=	Mesh.Node_Position[ 2 ][ No ] ;								
							}
							Mesh.Face_Position[ 0 ][ FaceNo ]	/=	Mesh.Face_Node[ FaceNo ].size() ;
							Mesh.Face_Position[ 1 ][ FaceNo ]	/=	Mesh.Face_Node[ FaceNo ].size() ;
							Mesh.Face_Position[ 2 ][ FaceNo ]	/=	Mesh.Face_Node[ FaceNo ].size() ;

							Check[ i ][ j ]				=	0 ;
							Check_BCFace[ BCFaceNo ]	=	0 ;

							FaceNo++ ;
							break ;
						}
					}
				}
			}
		}
	}

	for ( int i = 0 ; i < Processor_MeshCell[ comm_rank ].size() ; i++ )
	{
		CellNo =	Processor_MeshCell[ comm_rank ][ i ] ;

		for ( int j = 0 ; j < Mesh.Cell_FaceNum[ CellNo ] ; j++ )
		{
			FaceNo =	Mesh.Cell_Face[ CellNo ][ j ] ;

			if ( F[ FaceNo ] == 0 )
			{
				Processor_MeshFace[ comm_rank ].push_back( FaceNo ) ;
				Mesh.Face_ProcessorNo[ FaceNo ]	=	comm_rank ;
			}
			F[ FaceNo ]++ ;
		}
	}

	for ( int i = 0 ; i < Processor_GhostMeshCell_lv1[ comm_rank ].size() ; i++ )
	{
		CellNo =	Processor_GhostMeshCell_lv1[ comm_rank ][ i ] ;

		for ( int j = 0 ; j < Mesh.Cell_FaceNum[ CellNo ] ; j++ )
		{
			FaceNo =	Mesh.Cell_Face[ CellNo ][ j ] ;
			if ( F[ FaceNo ] == 0 )
			{
				Processor_GhostMeshFace_lv1[ comm_rank ].push_back( FaceNo ) ;
				Mesh.Face_ProcessorNo[ FaceNo ]	=	Mesh.Cell_ProcessorNo[ CellNo ] ;
			}
			F[ FaceNo ]++ ;
		}	
	}

	set_table_face( Mesh.FaceNum, Processor_MeshFace[ comm_rank ].size() + Processor_GhostMeshFace_lv1[ comm_rank ].size() ) ;

	for ( int i = 0 ; i < Processor_MeshFace[ comm_rank ].size() ; i++ )
	{
		FaceNo =	 Processor_MeshFace[ comm_rank ][ i ] ;

		MeshFace_LocalFaceNo[ FaceNo ]	=	i ;
		LocalFace_MeshFaceNo[ i ]		=	FaceNo ;
	}

	for ( int i = 0 ; i < Processor_GhostMeshFace_lv1[ comm_rank ].size() ; i++ )
	{
		No 		=	Processor_MeshFace[ comm_rank ].size() + i ;
		FaceNo 	=	Processor_GhostMeshFace_lv1[ comm_rank ][ i ] ;

		MeshFace_LocalFaceNo[ FaceNo ]	=	No ;
		LocalFace_MeshFaceNo[ No ]		=	FaceNo ;
	}
}

/*
void Domain::CreateFaceInformation( )
{
	int Type1, Type2, c1, c2, f1, f2, n1, m, CCNum ; 
	int NodeNo, BCFaceNo, CellNo, FaceNo, No ;
	int Check[ Mesh.CellNum ][ 6 ] ;
	map<int, int>	F ;
	CellMapping	pMapping ;

	for ( int i = 0 ; i < Mesh.CellNum ; i++ )
	{
		for ( int j = 0 ; j < 6 ; j++ )
		{
			Check[ i ][ j ] = -999 ;
		}
	}

	// InnerFace
	FaceNo =	0 ;
	for ( c1 = 0 ; c1 < Mesh.CellNum ; c1++ )
	{
		Type1 =	Mesh.Cell_Form[ c1 ] ;
		for ( f1 = 0 ; f1 < Mesh.Cell_FaceNum[ c1 ] ; f1++ )
		{
			if ( Check[ c1 ][ f1 ] == -999 )
			{
				CCNum = 0 ;
				while ( CCNum < Mesh.Cell_Cell[ c1 ].size() )
				{
					c2 		=	Mesh.Cell_Cell[ c1 ][ CCNum ] ;
					Type2 	=	Mesh.Cell_Form[ c2 ] ;

					for ( f2 = 0 ; f2 < Mesh.Cell_FaceNum[ c2 ] ; f2++ )
					{
						if ( CheckNeighbor( Type1, Type2, &Mesh.Cell_Node[ c1 ], &Mesh.Cell_Node[ c2 ], f1, f2, &Mesh.Face_Node[ FaceNo ], &pMapping ) )
						{
							n1 =	pMapping.Node[ Type1 ][ f1 ][ 0 ] ;
							n1 =	Mesh.Cell_Node[ c1 ][ n1 ] ;

							if ( Mesh.Cell_Type[ c1 ] != Mesh.Cell_Type[ c2 ] )
							{
								for ( int m = 0 ; m < Mesh.Node_BCFace[ n1 ].size() ; m++ )
								{
									BCFaceNo =	Mesh.Node_BCFace[ n1 ][ m ] ;

									if ( CheckNeighbor( Type1, &Mesh.Cell_Node[ c1 ], f1, &Mesh.BCFace_Node[ BCFaceNo ], &pMapping ) )
									{
										Mesh.Face_Type[ FaceNo ]		=	Mesh.BCFace_Type[ BCFaceNo ] ;
										Mesh.Face_Typename[ FaceNo ]	=	Mesh.BCFace_Typename[ BCFaceNo ] ;
									}
								}

								// without interface condition, create new one.
								if ( Mesh.Face_Type[ FaceNo ] == 0 )
								{
									Mesh.Face_Type[ FaceNo ]		=	Mesh.BC_TypeNum ;
									Mesh.Face_Typename[ FaceNo ]	=	"Interface" ;
								}
							} else
							{
								Mesh.Face_Type[ FaceNo ]		=	Mesh.Cell_Type[ c1 ] ;
								Mesh.Face_Typename[ FaceNo ]	=	Mesh.Cell_Typename[ c1 ] ;
							}

							Mesh.Face_Cell[ FaceNo ].push_back( c1 ) ;
							Mesh.Face_Cell[ FaceNo ].push_back( c2 ) ;

							Mesh.Cell_Face[ c1 ].push_back( FaceNo ) ;
							Mesh.Cell_Face[ c2 ].push_back( FaceNo ) ;

							Mesh.Face_Position[ 0 ][ FaceNo ] = 0 ;
							Mesh.Face_Position[ 1 ][ FaceNo ] = 0 ;
							Mesh.Face_Position[ 2 ][ FaceNo ] = 0 ;																		
							for ( int n = 0 ; n < Mesh.Face_Node[ FaceNo ].size() ; n++ )
							{
								NodeNo = Mesh.Face_Node[ FaceNo ][ n ] ;
								Mesh.Node_Face[ NodeNo ].push_back( FaceNo ) ;

								Mesh.Face_Position[ 0 ][ FaceNo ]	+=	Mesh.Node_Position[ 0 ][ NodeNo ] ;
								Mesh.Face_Position[ 1 ][ FaceNo ]	+=	Mesh.Node_Position[ 1 ][ NodeNo ] ;
								Mesh.Face_Position[ 2 ][ FaceNo ]	+=	Mesh.Node_Position[ 2 ][ NodeNo ] ;										
							}
							Mesh.Face_Position[ 0 ][ FaceNo ] /= Mesh.Face_Node[ FaceNo ].size() ;
							Mesh.Face_Position[ 1 ][ FaceNo ] /= Mesh.Face_Node[ FaceNo ].size() ;
							Mesh.Face_Position[ 2 ][ FaceNo ] /= Mesh.Face_Node[ FaceNo ].size() ;

							Check[ c1 ][ f1 ]	=	0 ;
							Check[ c2 ][ f2 ]	=	0 ;

							FaceNo++ ;
							CCNum =	Mesh.Cell_Cell[ c1 ].size() ;
							break ;
						}					
					}
					CCNum++ ;
				}
			}
		}
	}

	// BCFace
	for ( int i = 0 ; i < Mesh.CellNum ; i++ )
	{
		Type1	=	Mesh.Cell_Form[ i ] ;
		for ( int j = 0 ; j < Mesh.Cell_FaceNum[ i ] ; j++ )
		{
			if ( Check[ i ][ j ] == -999 )
			{
				NodeNo	=	pMapping.Node[ Type1 ][ j ][ 0 ] ;
				NodeNo	=	Mesh.Cell_Node[ i ][ NodeNo ] ;
				for ( int k = 0 ; k < Mesh.Node_BCFace[ NodeNo ].size() ; k++ )
				{
					BCFaceNo	=	Mesh.Node_BCFace[ NodeNo ][ k ] ;
					if ( CheckNeighbor( Type1, &Mesh.Cell_Node[ i ], j, &Mesh.BCFace_Node[ BCFaceNo ], &Mesh.Face_Node[ FaceNo ], &pMapping ) )
					{
						Mesh.Face_Type[ FaceNo ]		=	Mesh.BCFace_Type[ BCFaceNo ] ;
						Mesh.Face_Typename[ FaceNo ]	=	Mesh.BCFace_Typename[ BCFaceNo ] ;

						Mesh.Face_Cell[ FaceNo ].push_back( i ) ;
						Mesh.Cell_Face[ i ].push_back( FaceNo ) ;

						Mesh.Face_Position[ 0 ][ FaceNo ] = 0 ;
						Mesh.Face_Position[ 1 ][ FaceNo ] = 0 ;
						Mesh.Face_Position[ 2 ][ FaceNo ] = 0 ;	
						for ( int n = 0 ; n < Mesh.Face_Node[ FaceNo ].size() ; n++ )
						{
							No =	Mesh.Face_Node[ FaceNo ][ n ] ;
							Mesh.Node_Face[ No ].push_back( FaceNo ) ;

							Mesh.Face_Position[ 0 ][ FaceNo ]	+=	Mesh.Node_Position[ 0 ][ No ] ;
							Mesh.Face_Position[ 1 ][ FaceNo ]	+=	Mesh.Node_Position[ 1 ][ No ] ;
							Mesh.Face_Position[ 2 ][ FaceNo ]	+=	Mesh.Node_Position[ 2 ][ No ] ;								
						}
						Mesh.Face_Position[ 0 ][ FaceNo ]	/=	Mesh.Face_Node[ FaceNo ].size() ;
						Mesh.Face_Position[ 1 ][ FaceNo ]	/=	Mesh.Face_Node[ FaceNo ].size() ;
						Mesh.Face_Position[ 2 ][ FaceNo ]	/=	Mesh.Face_Node[ FaceNo ].size() ;

						Check[ i ][ j ]	=	0 ;
						FaceNo++ ;
						break ;
					}
				}
			}
		}
	}

	for ( int i = 0 ; i < Processor_MeshCell[ comm_rank ].size() ; i++ )
	{
		CellNo =	Processor_MeshCell[ comm_rank ][ i ] ;
		for ( int j = 0 ; j < Mesh.Cell_FaceNum[ i ] ;j++ )
		{
			FaceNo =	Mesh.Cell_Face[ CellNo ][ j ] ;
			if ( F[ FaceNo ] == 0 )
			{
				Processor_MeshFace[ comm_rank ].push_back( FaceNo ) ;
			}
			F[ FaceNo ]++ ;
		}
	}

	for ( int i = 0 ; i < Processor_GhostMeshCell_lv1[ comm_rank ].size() ; i++ )
	{
		CellNo =	Processor_GhostMeshCell_lv1[ comm_rank ][ i ] ;
		for ( int j = 0 ; j < Mesh.Cell_FaceNum[ CellNo ] ; j++ )
		{
			FaceNo =	Mesh.Cell_Face[ CellNo ][ j ] ;
			if ( F[ FaceNo ] == 0 )
			{
				Processor_GhostMeshFace_lv1[ comm_rank ].push_back( FaceNo ) ;
			}
			F[ FaceNo ]++ ;
		}	
	}

	set_table_face( Mesh.FaceNum, Processor_MeshFace[ comm_rank ].size() + Processor_GhostMeshFace_lv1[ comm_rank ].size() ) ;

	for ( int i = 0 ; i < Processor_MeshFace[ comm_rank ].size() ; i++ )
	{
		FaceNo =	 Processor_MeshFace[ comm_rank ][ i ] ;

		MeshFace_LocalFaceNo[ FaceNo ]	=	i ;
		LocalFace_MeshFaceNo[ i ]		=	FaceNo ;
	}

	for ( int i = 0 ; i < Processor_GhostMeshFace_lv1[ comm_rank ].size() ; i++ )
	{
		No 		=	Processor_MeshFace[ comm_rank ].size() + i ;
		FaceNo 	=	Processor_GhostMeshFace_lv1[ comm_rank ][ i ] ;

		MeshFace_LocalFaceNo[ FaceNo ]	=	No ;
		LocalFace_MeshFaceNo[ No ]		=	FaceNo ;
	}
}
*/

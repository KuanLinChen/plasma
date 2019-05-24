#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include "vtk_output.h"
#include "domain.h"
//#include "main.h"

/*!

@file tecplot_output.cpp Transfer memory data into Tecplot formate. The data should follow the domain structure.

*/

/*! \brief Tecplot format output.

@param zone_info The zone name in tecplot ASCII file
@param data_number Numbers of output data. Ex. ( X, Y, Z, V1, V1 ), data_number = 5
@param variable The pointer to the real data. Note that we only record the address of memory, which means the output data is the latest not the moment of when use assign address.
@param output_filename Filename
@param  base_domain The reference domain.
@is_new If ( is_new == true ), means we want to have a new output file no matter the file has content already, the old data will be clean. If ( is_new == false ), the data will append to the original one, with the new zone name  defined in 'zone_info'.

*/

void vtk_output( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
{
	int i, j, k, MeshNo ;
	int dimension = base_domain->dimension ;
	int mpi_rank, mpi_size ;
	int OutNum, BuffSize;
	int NodeNum = base_domain->global_node_number ;
	int CellNum = base_domain->global_cell_number ;
	MPI_Comm    comm ;
	MPI_Status  status ;
	ofstream	Output ;
	boost::shared_array <double> value ;

	mpi_size = base_domain->comm_size ;
	mpi_rank = base_domain->comm_rank ;
	comm = base_domain->comm ;

	if ( mpi_rank == 0 )
	{
		value = boost::shared_array <double> (  new double [ base_domain->global_cell_number] ) ;
		for ( i = 0 ; i < base_domain->global_cell_number ; i++ ) 
			value[i] = 0.;
		
		if ( is_new)
		{
			Output.open( output_filename.c_str(), ios::out | ios::trunc ) ;
		} else
		{
			Output.open( output_filename.c_str(), ios::out | ios::app ) ;
		}

		// Output paraview data
		// Hearder begin
		Output << "# vtk DataFile Version 1.0" << endl ;
		// Title
		Output << "Plasma properties" << endl ;
		// Data type, either ASCII or BINARY
		Output  << "ASCII" << endl ;
		// Dataset, UNSTRUCTURED_GRID only for now
		Output  << "DATASET UNSTRUCTURED_GRID" << endl << endl ;

		// POINTS n dataType
		Output  << "POINTS " << setw( 20 ) << NodeNum << setw( 20 ) << "float" << endl ;
        // coordinate: x y z
		OutNum = 0 ;
		for ( j = 0 ; j < NodeNum ; j++ )
		{
			OutNum++ ;
			Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 0 ][ j ] ;
			Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 1 ][ j ] ;
			Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 2 ][ j ] ;
			if ( OutNum % 5 == 0 ) Output << endl ;
		}
		Output << endl << endl ;

		// CELLS n size
		// n:       number of cell
		// size:    the total number of integer values required to represent the list
		// i.e., sum of numPoints and connectivity indices over each cell
		BuffSize = 0 ;
		for ( j = 0 ; j < CellNum ; j++ )
		{
			MeshNo = base_domain->GlobalCell_MeshCellNo[ j ] ;

			BuffSize = BuffSize + base_domain->Mesh.Cell_Node[ MeshNo ].size() + 1 ;
		}
		Output << "CELLS " << setw( 20 ) << CellNum << setw( 20 ) << BuffSize << '\n' ;
		// numPoints, point index i, point index j, point index k, ... of the cell
		OutNum = 0 ;
		for ( j = 0 ; j < CellNum ; j++ )
		{
			MeshNo = base_domain->GlobalCell_MeshCellNo[ j ] ;

			OutNum++ ;
			Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ].size() ;
			if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 6 ) 
			{
				// Pyramid Cell, the indices should be reordering
				Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] ;
				Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] ;
				Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] ;
				Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] ;
				Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] ;
				Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ 5 ] ;
			} else
			{
				// Tetra, Hexa, Prism; no need to reordering
				for ( k = 0 ; k < base_domain->Mesh.Cell_Node[ MeshNo ].size() ; k++ )
					Output << setw( 20 ) << base_domain->Mesh.Cell_Node[ MeshNo ][ k ] ;
			}
			if ( OutNum % 5 == 0 ) Output << endl ;
		}
		Output << endl << endl ;

		// CELL_TYPES n
		// n:       number of cell
		Output  << "CELL_TYPES " << setw( 20 ) << CellNum << '\n' ;
		// 5 for triangle, 8, 9 for quadrilatera
		OutNum = 0 ;
		for ( j = 0 ; j < CellNum ; j++ )
		{
			MeshNo = base_domain->GlobalCell_MeshCellNo[ j ] ;

			OutNum++ ;
			if ( dimension == 2  )
			{
				// Triangular Cell type: 5
				if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 3 ) Output << setw( 20 ) << "5"; 
				// Quadrilateral Cell type: 9
				else if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 4 ) Output << setw( 20 ) << "9"; 
			} else if ( dimension == 3 )
			{
				// Tetra Cell type: 10
				if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 4 ) Output << setw( 20 ) << "10"; 
				// Pyramid Cell type: 14
				else if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 5 ) Output << setw( 20 ) << "14"; 
				// Prism Cell type: 13
				else if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 6 ) Output << setw( 20 ) << "13"; 
				// Hexa Cell type: 12
				else if ( base_domain->Mesh.Cell_Node[ MeshNo ].size() == 8 ) Output << setw( 20 ) << "12"; 

			}
			if ( OutNum % 5 == 0 ) Output << endl ;
        }
		Output << endl << endl ;
	}

	// CELL_DATA n
	// CELL_DATA for element value; POINT_DATA for nodal value
	Output  << "CELL_DATA " << setw( 20 ) << CellNum << endl ;

	for ( k = 0 ; k < data_number ; k++  )
	{
		// Send and recieve

		if ( mpi_rank == 0 )
		{
			for ( i = 0 ; i < base_domain->local_cell_number ; i++ )
				value[ i ]	=	variable[k][i] ;

			j	=	base_domain->local_cell_number ;
			for ( i = 1 ; i < mpi_size ; i++ )
			{
				MPI_Recv( value.get() + j , base_domain->Processor_MeshCell[ i ].size(), MPI_DOUBLE, i, 1, comm, &status ) ;
				j	+=	base_domain->Processor_MeshCell[ i ].size() ;
			}
		} else
		{
			MPI_Send( variable[k], base_domain->local_cell_number, MPI_DOUBLE, 0, 1, comm ) ;
		}
		MPI_Barrier( comm ) ;

		// Dump data to file
		if ( mpi_rank == 0 )
		{
			// DataAttributes Name Type
			// DataAttributes: SCALARS, VECTORS, ...
			Output  << "SCALARS " << setw( 20 ) << field_name[ k ] << setw( 20 ) << "float" << endl ;
			// LOOKUP_TABLE default, necessary for SCALARS
			Output  << "LOOKUP_TABLE default" << endl ;
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 5 == 0 ) Output << endl ;
			}
			Output << endl << endl ;
		}
	}
	Output << endl << endl ;

}

//void tecplot_output_root ( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
//{
//	cout << "vtk_output.cpp: vtk_output_root" << endl;
//}


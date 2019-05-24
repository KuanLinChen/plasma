#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include "tecplot_output.h"
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
@param base_domain The reference domain.
@is_new If ( is_new == true ), means we want to have a new output file no matter the file has content already, the old data will be clean. If ( is_new == false ), the data will append to the original one, with the new zone name  defined in 'zone_info'.

*/

void tecplot_output ( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
{
	int i, j, k, MeshNo ;
	int dimension ;
	int mpi_rank, mpi_size ;
	int OutNum;
	int NodeNum = base_domain->global_node_number ;
	int CellNum = base_domain->global_cell_number ;
	MPI_Comm 	comm ;
	MPI_Status	status ;
	ofstream	Output ;
	boost::shared_array <double> value ;

	mpi_size = base_domain->comm_size ;
	mpi_rank = base_domain->comm_rank ;
	comm = base_domain->comm ;
	dimension = base_domain->dimension ;


	value =	boost::shared_array <double> (  new double [ base_domain->global_cell_number] ) ;

	if ( mpi_rank == 0 )
	{
		// Output tecplot data
		// Hearder begin
		if ( is_new)
		{
			Output.open( output_filename.c_str(), ios::out | ios::trunc ) ;
		} else
		{
			Output.open( output_filename.c_str(), ios::out | ios::app ) ;
		}

		Output << "TITLE = \"Finite volume dataset\"" << endl ;
		if ( dimension == 2 )
		{
			Output << "VARIABLES = \"X\", \"Y\"" ;
		} else if ( dimension == 3 )
		{
			Output << "VARIABLES = \"X\", \"Y\", \"Z\"" ;
		}

		for ( i = 0 ; i < data_number ; i++ )
			Output << ",\"" << field_name[i] << "\"" ;
		Output << endl ;

		if ( dimension == 2)
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-2]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 3 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		} else if ( dimension == 3 )
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-3]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 4 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2, 3]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		}

		if ( is_new )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 0 ][ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl  ;

			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 1 ][ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;

			if ( dimension == 3 )
			{
				OutNum = 0 ;
				for ( j = 0 ; j < NodeNum ; j++ )
				{
					OutNum++ ;
					Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 2 ][ j ] ;
					if ( OutNum % 6 == 0 ) Output << endl ;
				}
				Output << endl ;
			}
		}

		// Header end
	}
	Output << endl ;

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
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;
		}
	}
	Output << endl ;

	// Tail begin
	if ( mpi_rank == 0 && is_new )
	{
		if ( dimension == 2 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{
				MeshNo =	base_domain->GlobalCell_MeshCellNo[ j ] ;

				for ( k = 0 ; k < base_domain->Mesh.Cell_Node[ MeshNo ].size() ; k++ )
					Output << setw( 20 ) << ( base_domain->Mesh.Cell_Node[ MeshNo ][ k ] + 1 ) ;
				if (  base_domain->Mesh.Cell_Node[ MeshNo ].size() == 3 )
					Output << setw( 20 ) << ( base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ) ;
				Output << endl ;
			}
			Output << endl ;
			Output.clear() ;
			Output.close() ;
		} else if ( dimension == 3 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{
				MeshNo =	base_domain->GlobalCell_MeshCellNo[ j ] ;

				// Tetra Cell.
				if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 0 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
				// Hex. Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 1 )
				{
					for ( k = 0 ; k < 8 ; k++ )
						Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ k ] + 1 ;
				// Prism Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 2 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 5 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 5 ] + 1 ;					
				// Pyramid Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 3 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;										
				// 2D Triangular Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 4 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;					
				// 2D Quadrilateral Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 5 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;						
				}

				Output << endl ;
			}
		}
	} else
	{
		Output << endl ;
		Output.clear() ;
		Output.close() ;
	}
	// Tail end
	MPI_Barrier( comm ) ;
} ;


// Out put only data from root mpi
void tecplot_output_root ( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
{
	int i, j, k, MeshNo ;
	int dimension ;
	int mpi_rank, mpi_size ;
	int OutNum;
	int NodeNum = base_domain->global_node_number ;
	int CellNum = base_domain->global_cell_number ;
	MPI_Comm 	comm ;
	MPI_Status	status ;
	ofstream	Output ;
	boost::shared_array <double> value ;

	mpi_size = base_domain->comm_size ;
	mpi_rank = base_domain->comm_rank ;
	comm = base_domain->comm ;
	dimension = base_domain->dimension ;

	value =	boost::shared_array <double> (  new double [ base_domain->global_cell_number] ) ;

	if ( mpi_rank == 0 )
	{
		// Output tecplot data
		// Hearder begin
		if ( is_new)
		{
			Output.open( output_filename.c_str(), ios::out | ios::trunc ) ;
		} else
		{
			Output.open( output_filename.c_str(), ios::out | ios::app ) ;
		}

		Output << "TITLE = \"Finite volume dataset\"" << endl ;
		if ( dimension == 2 )
		{
			Output << "VARIABLES = \"X\", \"Y\"" ;
		} else if ( dimension == 3 )
		{
			Output << "VARIABLES = \"X\", \"Y\", \"Z\"" ;
		}

		for ( i = 0 ; i < data_number ; i++ )
			Output << ",\"" << field_name[i] << "\"" ;
		Output << endl ;

		if ( dimension == 2)
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-2]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 3 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		} else if ( dimension == 3 )
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-3]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 4 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2, 3]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		}

		if ( is_new )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 0 ][ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl  ;

			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 1 ][ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;

			if ( dimension == 3 )
			{
				OutNum = 0 ;
				for ( j = 0 ; j < NodeNum ; j++ )
				{
					OutNum++ ;
					Output << setw( 20 ) <<  base_domain->Mesh.Node_Position[ 2 ][ j ] ;
					if ( OutNum % 6 == 0 ) Output << endl ;
				}
				Output << endl ;
			}
		}

		// Header end
	}
	Output << endl ;

	for ( k = 0 ; k < data_number ; k++  )
	{
		// Send and recieve
		if ( mpi_rank == 0 )
		{
			for ( i = 0 ; i < base_domain->global_cell_number ; i++ )
				value[ i ]	=	variable[k][i] ;
		}

		// Dump data to file
		if ( mpi_rank == 0 )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;
		}
	}
	Output << endl ;

	// Tail begin
	if ( mpi_rank == 0 && is_new )
	{
		if ( dimension == 2 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{
				MeshNo =	base_domain->GlobalCell_MeshCellNo[ j ] ;

				for ( k = 0 ; k < base_domain->Mesh.Cell_Node[ MeshNo ].size() ; k++ )
					Output << setw( 20 ) << ( base_domain->Mesh.Cell_Node[ MeshNo ][ k ] + 1 ) ;
				if (   base_domain->Mesh.Cell_Node[ MeshNo ].size() == 3 )
					Output << setw( 20 ) << ( base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ) ; // base_domain->cell_node[ j ][ 2 ] + 1  ;
				Output << endl ;
			}
			Output << endl ;
			Output.clear() ;
			Output.close() ;
		} else if ( dimension == 3 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{
				MeshNo =	base_domain->GlobalCell_MeshCellNo[ j ] ;

				// Tetra Cell.
				if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 0 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
				// Hex. Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 1 )
				{
					for ( k = 0 ; k < 8 ; k++ )
						Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ k ] + 1 ;
				// Prism Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 2 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 5 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 5 ] + 1 ;					
				// Pyramid Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 3 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 4 ] + 1 ;										
				// 2D Triangular Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 4 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;					
				// 2D Quadrilateral Cell
				} else if ( base_domain->Mesh.Cell_Form[ MeshNo ] == 5 )
				{
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 0 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 1 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 2 ] + 1 ;
					Output << setw(12) << base_domain->Mesh.Cell_Node[ MeshNo ][ 3 ] + 1 ;						
				}

				Output << endl ;
			}
		}
	} else
	{
		Output << endl ;
		Output.clear() ;
		Output.close() ;
	}
	// Tail end
	MPI_Barrier( comm ) ;
};


/*
void tecplot_output ( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
{
	int i, j, k;
	int dimension ;
	int mpi_rank, mpi_size ;
	int OutNum;
	int NodeNum = base_domain->global_node_number ;
	int CellNum = base_domain->global_cell_number ;
	MPI_Comm 	comm ;
	MPI_Status	status ;
	ofstream	Output ;
	boost::shared_array <double> value ;

	mpi_size = base_domain->comm_size ;
	mpi_rank = base_domain->comm_rank ;
	comm = base_domain->comm ;
	dimension = base_domain->dimension ;


	value =	boost::shared_array <double> (  new double [ base_domain->global_cell_number] ) ;

	if ( mpi_rank == 0 )
	{
		// Output tecplot data
		// Hearder begin
		if ( is_new)
		{
			Output.open( output_filename.c_str(), ios::out | ios::trunc ) ;
		}
		else
		{
			Output.open( output_filename.c_str(), ios::out | ios::app ) ;
		}
		Output << "TITLE = \"Finite volume dataset\"" << endl ;
		if ( dimension == 2 )
		{
			Output << "VARIABLES = \"X\", \"Y\"" ;
		}
		else if ( dimension == 3 )
		{
			Output << "VARIABLES = \"X\", \"Y\", \"Z\"" ;
		}

		for ( i = 0 ; i < data_number ; i++ )
			Output << ",\"" << field_name[i] << "\"" ;
		Output << endl ;

		if ( dimension == 2)
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-2]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 3 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		}
		else if ( dimension == 3 )
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-3]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 4 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2, 3]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		}

		if ( is_new )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->pre_Node[ j ].x ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl  ;

			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->pre_Node[ j ].y ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;

			if ( dimension == 3 )
			{
				OutNum = 0 ;
				for ( j = 0 ; j < NodeNum ; j++ )
				{
					OutNum++ ;
					Output << setw( 20 ) <<  base_domain->pre_Node[ j ].z ;
					if ( OutNum % 6 == 0 ) Output << endl ;
				}
				Output << endl ;
			}
		}

		// Header end
	}
	Output << endl ;

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
				MPI_Recv( value.get() + j , base_domain->Processor[i].Cell.size(), MPI_DOUBLE, i, 1, comm, &status ) ;
				j	+=	base_domain->Processor[i].Cell.size() ;
			}
		}
		else
		{
			MPI_Send( variable[k], base_domain->local_cell_number, MPI_DOUBLE, 0, 1, comm ) ;
		}
		MPI_Barrier( comm ) ;

		// Dump data to file
		if ( mpi_rank == 0 )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;
		}
	}
	Output << endl ;

	// Tail begin
	if ( mpi_rank == 0 && is_new )
	{
		if ( dimension == 2 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{
				for ( k = 0 ; k < base_domain->pre_Cell[ j ].NodeNum ; k++ )
					Output << setw( 20 ) << ( base_domain->pre_Cell[ j ].pNode[k]->id + 1 ) ;
				if (   base_domain->pre_Cell[ j ].NodeNum == 3 )
					Output << setw( 20 ) << ( base_domain->pre_Cell[ j ].pNode[2]->id + 1 ) ; // base_domain->cell_node[ j ][ 2 ] + 1  ;
				Output << endl ;
			}
			Output << endl ;
			Output.clear() ;
			Output.close() ;
		}
		else if ( dimension == 3 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{

				// Tetra Cell.
				if ( base_domain->pre_Cell[j].Type == 0 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					// Hex. Cell
				}
				else if ( base_domain->pre_Cell[j].Type == 1 )
				{
					for ( k = 0 ; k < 8 ; k++ )
						Output << setw(12) << base_domain->pre_Cell[j].pNode[k]->id + 1 ;
					// Pyramid Cell
				}
				else if ( base_domain->pre_Cell[j].Type == 2 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					// Prism Cell
				}
				else if (  base_domain->pre_Cell[j].Type == 3 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 5 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 5 ]->id + 1 ;
					// 2D Triangular Cell
				}
				else if (  base_domain->pre_Cell[j].Type == 4 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					// 2D Quadrilateral Cell
				}
				else if (  base_domain->pre_Cell[j].Type == 5 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
				}

				Output << endl ;
			}
		}
	}
	else
	{
		Output << endl ;
		Output.clear() ;
		Output.close() ;
	}
	// Tail end
	MPI_Barrier( comm ) ;
};



// Out put only data from root mpi
void tecplot_output_root ( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
{
	int i, j, k;
	int dimension ;
	int mpi_rank, mpi_size ;
	int OutNum;
	int NodeNum = base_domain->global_node_number ;
	int CellNum = base_domain->global_cell_number ;
	MPI_Comm 	comm ;
	MPI_Status	status ;
	ofstream	Output ;
	boost::shared_array <double> value ;

	mpi_size = base_domain->comm_size ;
	mpi_rank = base_domain->comm_rank ;
	comm = base_domain->comm ;
	dimension = base_domain->dimension ;


	value =	boost::shared_array <double> (  new double [ base_domain->global_cell_number] ) ;

	if ( mpi_rank == 0 )
	{
		// Output tecplot data
		// Hearder begin
		if ( is_new)
		{
			Output.open( output_filename.c_str(), ios::out | ios::trunc ) ;
		}
		else
		{
			Output.open( output_filename.c_str(), ios::out | ios::app ) ;
		}
		Output << "TITLE = \"Finite volume dataset\"" << endl ;
		if ( dimension == 2 )
		{
			Output << "VARIABLES = \"X\", \"Y\"" ;
		}
		else if ( dimension == 3 )
		{
			Output << "VARIABLES = \"X\", \"Y\", \"Z\"" ;
		}

		for ( i = 0 ; i < data_number ; i++ )
			Output << ",\"" << field_name[i] << "\"" ;
		Output << endl ;

		if ( dimension == 2)
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-2]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 3 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		}
		else if ( dimension == 3 )
		{
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << endl ;
			Output << "T = \"" << zone_info << "\"" << endl;
			Output << "VARLOCATION=([1-3]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 4 << "]=CELLCENTERED" ;
			Output << ")" << endl ;
			if ( is_new == false ) Output << "VARSHARELIST = ([1, 2, 3]=1), CONNECTIVITYSHAREZONE = 1" << endl ;
		}

		if ( is_new )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->pre_Node[ j ].x ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl  ;

			OutNum = 0 ;
			for ( j = 0 ; j < NodeNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->pre_Node[ j ].y ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;

			if ( dimension == 3 )
			{
				OutNum = 0 ;
				for ( j = 0 ; j < NodeNum ; j++ )
				{
					OutNum++ ;
					Output << setw( 20 ) <<  base_domain->pre_Node[ j ].z ;
					if ( OutNum % 6 == 0 ) Output << endl ;
				}
				Output << endl ;
			}
		}

		// Header end
	}
	Output << endl ;

	for ( k = 0 ; k < data_number ; k++  )
	{
		// Send and recieve

		if ( mpi_rank == 0 )
		{
			for ( i = 0 ; i < base_domain->global_cell_number ; i++ )
				value[ i ]	=	variable[k][i] ;
		}

		// Dump data to file
		if ( mpi_rank == 0 )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;
		}
	}
	Output << endl ;

	// Tail begin
	if ( mpi_rank == 0 && is_new )
	{
		if ( dimension == 2 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{
				for ( k = 0 ; k < base_domain->pre_Cell[ j ].NodeNum ; k++ )
					Output << setw( 20 ) << ( base_domain->pre_Cell[ j ].pNode[k]->id + 1 ) ;
				if (   base_domain->pre_Cell[ j ].NodeNum == 3 )
					Output << setw( 20 ) << ( base_domain->pre_Cell[ j ].pNode[2]->id + 1 ) ; // base_domain->cell_node[ j ][ 2 ] + 1  ;
				Output << endl ;
			}
			Output << endl ;
			Output.clear() ;
			Output.close() ;
		}
		else if ( dimension == 3 )
		{
			for ( j = 0 ; j < CellNum ; j++ )
			{

				// Tetra Cell.
				if ( base_domain->pre_Cell[j].Type == 0 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					// Hex. Cell
				}
				else if ( base_domain->pre_Cell[j].Type == 1 )
				{
					for ( k = 0 ; k < 8 ; k++ )
						Output << setw(12) << base_domain->pre_Cell[j].pNode[k]->id + 1 ;
					// Pyramid Cell
				}
				else if ( base_domain->pre_Cell[j].Type == 2 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					// Prism Cell
				}
				else if (  base_domain->pre_Cell[j].Type == 3 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 4 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 5 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 5 ]->id + 1 ;
					// 2D Triangular Cell
				}
				else if (  base_domain->pre_Cell[j].Type == 4 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					// 2D Quadrilateral Cell
				}
				else if (  base_domain->pre_Cell[j].Type == 5 )
				{
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 0 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 1 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 2 ]->id + 1 ;
					Output << setw(12) << base_domain->pre_Cell[j].pNode[ 3 ]->id + 1 ;
				}

				Output << endl ;
			}
		}
	}
	else
	{
		Output << endl ;
		Output.clear() ;
		Output.close() ;
	}
	// Tail end
	MPI_Barrier( comm ) ;
};
*/
void tecplot_output1D ( string zone_info, int data_number, string *field_name, double **variable, string output_filename, Domain *base_domain , bool is_new )
{
	int i, j, k, MeshNo ;
	int dimension ;
	int mpi_rank, mpi_size ;
	int OutNum;
	int NodeNum = base_domain->global_node_number ;
	int CellNum = base_domain->global_cell_number ;
	MPI_Comm 	comm ;
	MPI_Status	status ;
	ofstream	Output ;
	boost::shared_array <double> value ;

	mpi_size = base_domain->comm_size ;
	mpi_rank = base_domain->comm_rank ;
	comm = base_domain->comm ;
	dimension = base_domain->dimension ;


	value =	boost::shared_array <double> (  new double [ base_domain->global_cell_number] ) ;

	if ( mpi_rank == 0 )
	{
		// Output tecplot data
		// Hearder begin
		if ( is_new)
		{
			Output.open( output_filename.c_str(), ios::out | ios::trunc ) ;

		} else
		{
			Output.open( output_filename.c_str(), ios::out | ios::app ) ;
		}

		Output << "TITLE = \"Finite volume dataset\"" << endl ;

		Output << "VARIABLES = \"Position\"" ;

		for ( i = 0 ; i < data_number ; i++ )
			Output << ",\"" << field_name[i] << "\"" ;

			Output << endl ;

			Output << "ZONE I=" << CellNum <<", J=1" << ", DATAPACKING=BLOCK" << endl ;
			Output << "VARLOCATION=([1]=NODAL" ;
			for ( i = 0 ; i < data_number ; i++ )
				Output << " ,[" << i + 2 << "]=NODAL" ;
			Output << ")" << endl ;

		if ( is_new )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) <<  base_domain->Mesh.Cell_Position[ 0 ][ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl  ;
		}
		// Header end
	}
	Output << endl ;

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
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				MeshNo =	base_domain->GlobalCell_MeshCellNo[ j ] + 1 ;
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;
		}
	}
	Output << endl ;


	for ( k = 0 ; k < data_number ; k++  )
	{
		// Send and recieve
		if ( mpi_rank == 0 )
		{
			for ( i = 0 ; i < base_domain->global_cell_number ; i++ )
				value[ i ]	=	variable[k][i] ;
		}

		// Dump data to file
		if ( mpi_rank == 0 )
		{
			OutNum = 0 ;
			for ( j = 0 ; j < CellNum ; j++ )
			{
				OutNum++ ;
				Output << setw( 20 ) << value[ j ] ;
				if ( OutNum % 6 == 0 ) Output << endl ;
			}
			Output << endl ;
		}
	}
	Output << endl ;

	// Tail begin
	if ( mpi_rank == 0 && is_new )
	{
		Output << endl ;
		Output.clear() ;
		Output.close() ;
	} else
	{
		Output << endl ;
		Output.clear() ;
		Output.close() ;
	}
	// Tail end
	MPI_Barrier( comm ) ;
} ;
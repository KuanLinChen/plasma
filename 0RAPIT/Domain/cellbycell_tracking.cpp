#include <iostream>
#include <mpi.h>
#include <cmath>
#include "cellbycell_tracking.h"
#include "domain.h"
#include "sys_log.h"

using namespace std;


CellByCellRayTracking::CellByCellRayTracking( Domain *BseaG )
{
	int i ;
	BaseGrid        =	BseaG ;
	NumOfRandTry    =	1 ;
	nTry            =	0 ;

	MPI_Comm_size( BaseGrid->comm, &mpi_size ) ;
	MPI_Comm_rank( BaseGrid->comm, &mpi_id ) ;

	MainID 			=	mpi_id ;
	iNumDataSize 	=	9 ;
	iDataNum 		=	4 ;
	dDataNum 		=	6 ;

	SendOffSet 		=	boost::shared_array< vector<int> > ( new vector<int> [ iNumDataSize ] ) ;
	RecvOffSet 		=	boost::shared_array< vector<int> > ( new vector<int> [ iNumDataSize ] ) ;

	for( i = 0 ; i < mpi_size + 1 ; i++ ) 
	{
		SendOffSet[ 1 ].push_back( 0 ) ;
		RecvOffSet[ 1 ].push_back( 0 ) ;
	}

	//for ( int i = 0 ; i < iNumDataSize ; i++ )
	//{
	//	SendOffSet[ i ] = boost::shared_array<int> > ( new int [ mpi_size + 1 ] ) ;
	//	RecvOffSet[ i ] = boost::shared_array<int> > ( new int [ mpi_size + 1 ] ) ;
	//}
		
	iNumSend 	=	boost::shared_array<int> ( new int[ mpi_size ] ) ;
	iNumRecv 	=	boost::shared_array<int> ( new int[ mpi_size ] ) ;

	CellShift 	=	boost::shared_array<int> ( new int[ mpi_size + 1 ] ) ;

	MPI_Allgather( &(BaseGrid->local_cell_number), 1, MPI_INT, CellShift.get(), 1, MPI_INT, BaseGrid->comm ) ;

	for ( int i = 1 ; i < mpi_size ; i++ )		CellShift[ i ] += CellShift[ i - 1 ] ;
	for ( int i = mpi_size - 1 ; i > 0 ; i-- )	CellShift[ i ]  = CellShift[ i - 1 ] ;

	CellShift[ mpi_size ]	=	BaseGrid->global_cell_number ;
	CellShift[ 0 ]			=	0 ;
}

CellByCellRayTracking::~CellByCellRayTracking()
{
	/*
	    for( int i=0; i <mpi_size+1; i++) delete [] SendOffSet[i];
	    for( int i=0; i <mpi_size+1; i++) delete [] RecvOffSet[i];

	    delete [] SendOffSet;
	    delete [] RecvOffSet;

	    delete [] iNumSend;
	    delete [] iNumRecv;

	    delete [] CellShift;

	    if(PntID        != NULL) delete [] PntID;
	    if(ToCellID     != NULL) delete [] ToCellID;
	    if(FromCellID   != NULL) delete [] FromCellID;
	    if(CpuID        != NULL) delete [] CpuID;
	    if(SendTag      != NULL) delete [] SendTag;
	    for(int i=0; i< 3; i++){
	        if(PntPosDest[i] != NULL) delete [] PntPosDest[i];
	        if(PntPosFrom[i] != NULL) delete [] PntPosFrom[i];
	    }
	*/

}


CellByCellRayTracking & CellByCellRayTracking::operator << ( int iTry )
{
	NumOfRandTry = iTry ;
	return *this ;
} ;

CellByCellRayTracking & CellByCellRayTracking::operator << ( Domain *Grid2 )
{
	InGrid 		=	Grid2 ;
	iNumPnts 	=	InGrid->local_cell_number ;

	//if (PntID != NULL)
	//{
		//delete [] PntID;
		//delete [] ToCellID;
		//delete [] FromCellID;
		//delete [] CpuID;
		//delete [] PntPosFrom[0];
		//delete [] PntPosFrom[1];
		//delete [] PntPosFrom[2];
		//delete [] SendTag;
		//delete [] PntPosDest[0];
		//delete [] PntPosDest[1];
		//delete [] PntPosDest[2];
	//}

	//PntID 		= new int [ iNumPnts ] ;
	//ToCellID 		= new int [ iNumPnts ] ;
	//FromCellID 	= new int [ iNumPnts ] ;
	//CpuID 		= new int [ iNumPnts ] ;
	PntPosFrom[ 0 ]	= boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosFrom[ 1 ]	= boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosFrom[ 2 ]	= boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	//SendTag		= new int [ iNumPnts ] ;
	PntID           = boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	ToCellID        = boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	FromCellID      = boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	CpuID           = boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	SendTag         = boost::shared_array<int>  ( new int [ iNumPnts ] ) ;

	PntPosDest[ 0 ]	= boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosDest[ 1 ]	= boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosDest[ 2 ]	= boost::shared_array<double> ( new double [ iNumPnts ] ) ;

	if ( InGrid != NULL )
	{
		for ( int i = 0 ; i < iNumPnts ; i++ )
		{
			PntPosDest[ 0 ][ i ] = InGrid->cell[ i ].x ;
			PntPosDest[ 1 ][ i ] = InGrid->cell[ i ].y ;
			PntPosDest[ 2 ][ i ] = InGrid->cell[ i ].z ;
		}
	}

	for ( int i = 0; i < iNumPnts ; i++)
	{
		if ( InGrid != NULL )
		{
			PntID[ i ]	=	InGrid->cell[ i ].id ;
		} else
		{
			PntID[ i ]	=	i ;
		}

		CpuID[ i ]				=	mpi_id ;
		ToCellID[ i ]			=	-1 ;
		FromCellID[ i ]			=	BaseGrid->cell[ 0 ].id ;
		PntPosFrom[ 0 ][ i ]	=	BaseGrid->cell[ 0 ].x ;
		PntPosFrom[ 1 ][ i ]	=	BaseGrid->cell[ 0 ].y ;
		PntPosFrom[ 2 ][ i ]	=	BaseGrid->cell[ 0 ].z ;
	}

	return *this ;
} ;

void CellByCellRayTracking::PointsMapping( int iNumPnt, double *Xoc, double *Yoc, double *Zoc, int *CellIDMap )
{
	
	std::string filename;
	std::ostringstream Sconvert;

	Sconvert.clear();
	Sconvert.str("");
	Sconvert << mpi_id ;
	filename = "debug"+ Sconvert.str() + ".txt";

	ofstream degug ;
	debug.open(filename.c_str()) ;
	

	int 	i ;
	InGrid 			=	NULL ;
	iNumPnts 		=	iNumPnt ;

	PntID 			=	boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	ToCellID 		=	boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	FromCellID 		=	boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	CpuID 			=	boost::shared_array<int>  ( new int [ iNumPnts ] ) ;
	PntPosFrom[ 0 ]	=	boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosFrom[ 1 ]	=	boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosFrom[ 2 ]	=	boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	SendTag 		=	boost::shared_array<int>  ( new int [ iNumPnts ] ) ;

	PntPosDest[ 0 ]	=	boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosDest[ 1 ]	=	boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	PntPosDest[ 2 ]	=	boost::shared_array<double> ( new double [ iNumPnts ] ) ;
	
	for ( i = 0 ; i < iNumPnts ; i++ )
	{
		PntPosDest[ 0 ][ i ]	=	Xoc[ i ] ;
		PntPosDest[ 1 ][ i ]	=	Yoc[ i ] ;
		PntPosDest[ 2 ][ i ]	=	Zoc[ i ] ;
	}

	for ( i = 0 ; i < iNumPnts ; i++ )
	{
		PntID[ i ]				=	i ;
		CpuID[ i ]				=	mpi_id ;
		ToCellID[ i ]			=	-1 ;
		FromCellID[ i ]			=	BaseGrid->cell[ 0 ].id ;
		PntPosFrom[ 0 ][ i ]	=	BaseGrid->cell[ 0 ].r[ 0 ] ;
		PntPosFrom[ 1 ][ i ]	=	BaseGrid->cell[ 0 ].r[ 1 ] ;
		PntPosFrom[ 2 ][ i ]	=	BaseGrid->cell[ 0 ].r[ 2 ] ;
	}

	int NumOfNotFound = 0, NotFound = 0 ;
	int NumOfHitBndy = 0,  HitBndy = 0 ;
	int NotFoundBKP, HitBKP, TotUmMaped, TotUmMapedBKP, TotUmMapedBKP2 ;
	int iNumTotalPnt ;

	//nTry = 0 ;
	//for( int i = 0 ; i< iNumPnts ; i++ ) if( ToCellID[ i ] < 0 ) NotFound++ ;
	//MPI_Allreduce( &NotFound, &NumOfNotFound, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD ) ;
	//MPI_Allreduce( &iNumPnts, &iNumTotalPnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD ) ;
	//if( mpi_id == 0 ) cout << "To be mapped:" << NumOfNotFound << "/" << iNumTotalPnt << endl ;
	//if( NumOfNotFound == 0 ) return ;

	int itag ;
	for ( int jj = 0 ; jj < 40 ; jj++ )
	{
		itag =	mpi_size * jj + 1 ;
		*this << itag + 1 ;
		nTry = 0 ;
		NotFound = 0 ;

		for ( i = 0 ; i < iNumPnts ; i++ ) if ( ToCellID[ i ] < 0 ) NotFound++ ;
		MPI_Allreduce( &NotFound, &NumOfNotFound, 1, MPI_INT, MPI_SUM, BaseGrid->comm ) ;
		MPI_Allreduce( &iNumPnts, &iNumTotalPnt, 1, MPI_INT, MPI_SUM, BaseGrid->comm ) ;

		//if ( mpi_id == 0 ) cout << "To be mapped:" << NumOfNotFound << "/" << iNumTotalPnt << endl ;
		Log().TagDump( logLEVEL4 ) << "To be mapped:" << NumOfNotFound << "/" << iNumTotalPnt ;
		if ( NumOfNotFound == 0 )
		{
			for ( i = 0 ; i < iNumPnts ; i++ )
			{
				CellIDMap[ i ]	=	ToCellID[ i ] ;
			}
			return ;
		}

		for (  i = 0 ; i < itag ; i++ )
		{
            NotFoundBKP =	NumOfNotFound ;
            HitBKP      =	NumOfHitBndy ;

			Tracking( iNumPnts, PntPosDest[ 0 ].get(), PntPosDest[ 1 ].get(), PntPosDest[ 2 ].get(), PntPosFrom[ 0 ].get(), PntPosFrom[ 1 ].get(), PntPosFrom[ 2 ].get(), PntID.get(), ToCellID.get(), FromCellID.get(), CpuID.get() ) ;
			CountingBuffer() ;
			SendToBuffer() ;
			Tracking( iNumTotRecv, PntPosDest_r[ 0 ], PntPosDest_r[ 1 ], PntPosDest_r[ 2 ], PntPosFrom_r[ 0 ], PntPosFrom_r[ 1 ], PntPosFrom_r[ 2 ], PntID_r, ToCellID_r, FromCellID_r, CpuID_r ) ;
			ReturnBuffer( &NumOfNotFound, &NumOfHitBndy ) ;
			TotUmMaped = NumOfNotFound + NumOfHitBndy ;

			if ( ( NotFoundBKP == NumOfNotFound && HitBKP == NumOfHitBndy && nTry >= NumOfRandTry ) || ( NumOfNotFound == 0 && NumOfHitBndy == 0 ) ) break ;

			if ( NotFoundBKP == NumOfNotFound && HitBKP == NumOfHitBndy )
			{
				//if( TotUmMapedBKP2 == TotUmMaped ) break ;
				nTry++ ;
				MainID = nTry + mpi_id ;
				while ( MainID > NumOfRandTry )
				{
					MainID -= NumOfRandTry ;
				}

				int cid = MainID * double( BaseGrid->global_cell_number ) / double( NumOfRandTry ) - 1 ;
				//int gcid = BaseGrid->pCell_Cell[ cid ] ;
				int gcid	=	BaseGrid->GlobalCell_LocalCellNo[ cid ] ;
				int pid 	=	0 ;
				int gpid 	=	0 ;

				for ( int j = 0 ; j < mpi_size ; j++ )
				{
					if ( cid >= CellShift[ j ] )
					{
						pid = j ;
					} else
					{
						//cout <<"myid = " << mpi_id << "cid = " << cid << ", shift = " <<  CellShift[ j ] << endl ;
						break ;
					}
				}

				//gpid = BaseGrid->pre_Cell[ cid ].ProcessorNo ;
				gpid =	BaseGrid->Mesh.Cell_ProcessorNo[ BaseGrid->GlobalCell_MeshCellNo[ cid ] ]  ;
				if ( NumOfRandTry > nTry )
				{
					for ( i = 0 ; i < iNumPnts ; i++ )
					{
						if ( ToCellID[ i ] < 0 )
						{
							ToCellID[ i ]	=	-1 ;
							CpuID[ i ]		=	pid ;
							FromCellID[ i ]	=	cid ;
						}
					}
				}
			}
		}

		//if( mpi_id == 0 ) cout << "nTry: " << nTry << " Total un-mapped: " << NumOfNotFound + NumOfHitBndy; // << endl;
		//if( mpi_id == 0 ) cout << ", To be mapped:" << NumOfNotFound << "/" << iNumTotalPnt ;//<< ", after " << i+1 <<"th mapping";// << endl;
		//if( mpi_id == 0 ) cout << ", Hit Boundary:" << NumOfHitBndy << "/" << iNumTotalPnt << endl; //<< ", after " << i+1 <<"th mapping" << endl;

		//if( mpi_id == 0 ) (*BaseGrid->debug_out) << "nTry: " << nTry << " Total un-mapped: " << NumOfNotFound + NumOfHitBndy; // << endl;
		//if( mpi_id == 0 ) (*BaseGrid->debug_out) << ", To be mapped:" << NumOfNotFound << "/" << iNumTotalPnt ;//<< ", after " << i+1 <<"th mapping";// << endl;
		//if( mpi_id == 0 ) (*BaseGrid->debug_out) << ", Hit Boundary:" << NumOfHitBndy << "/" << iNumTotalPnt << endl; //<< ", after " << i+1 <<"th mapping" << endl;
	}

	for ( i = 0 ; i < iNumPnts ; i++ )
	{
		CellIDMap[ i ]	=	ToCellID[ i ] ;
	}
} ;


void CellByCellRayTracking::TryMap( int itag )
{
	int NumOfNotFound = 0, NotFound = 0 ;
	int NumOfHitBndy = 0,  HitBndy = 0 ;
	int NotFoundBKP, HitBKP, TotUmMaped, TotUmMapedBKP, TotUmMapedBKP2 ;

	nTry =	0 ;

	for ( int i = 0 ; i < InGrid->local_cell_number ; i++ ) if ( ToCellID[ i ] < 0 ) NotFound++ ;
	MPI_Allreduce( &NotFound, &NumOfNotFound, 1, MPI_INT, MPI_SUM, BaseGrid->comm ) ;
	// if ( mpi_id == 0 ) cout << "To be mapped:" << NumOfNotFound << "/" << InGrid->global_cell_number << endl ;
	Log().TagDump( logLEVEL4 )  << "To be mapped:" << NumOfNotFound << "/" << InGrid->global_cell_number ;
	if ( NumOfNotFound == 0 ) return ;

	for ( int i = 0 ; i < itag ; i++ )
	{
		NotFoundBKP 	=	NumOfNotFound ;
		HitBKP 			=	NumOfHitBndy ;
		//TotUmMapedBKP =	NotFoundBKP + HitBKP ;
		//if( itag == 0 )	TotUmMapedBKP2 = TotUmMapedBKP ;

		Tracking( InGrid->local_cell_number, PntPosDest[ 0 ].get(), PntPosDest[ 1 ].get(), PntPosDest[ 2 ].get(), PntPosFrom[ 0 ].get(), PntPosFrom[ 1 ].get(), PntPosFrom[ 2 ].get(), PntID.get(), ToCellID.get(), FromCellID.get(), CpuID.get() ) ;
		CountingBuffer() ;
		SendToBuffer() ;
		Tracking( iNumTotRecv, PntPosDest_r[ 0 ], PntPosDest_r[ 1 ], PntPosDest_r[ 2 ],
							   PntPosFrom_r[ 0 ], PntPosFrom_r[ 1 ], PntPosFrom_r[ 2 ], PntID_r, ToCellID_r, FromCellID_r, CpuID_r ) ;
		ReturnBuffer( &NumOfNotFound, &NumOfHitBndy ) ;
		TotUmMaped 	=	NumOfNotFound + NumOfHitBndy ;

		if ( ( NotFoundBKP == NumOfNotFound && HitBKP == NumOfHitBndy && nTry >= NumOfRandTry ) || ( NumOfNotFound == 0 && NumOfHitBndy == 0 ) ) break ;

		if ( NotFoundBKP == NumOfNotFound && HitBKP == NumOfHitBndy )
		{
			//if(TotUmMapedBKP2 == TotUmMaped ) break;
			nTry++ ;
			MainID =	nTry + mpi_id ;
			while ( MainID > NumOfRandTry )
			{
				MainID -=	NumOfRandTry ;
			}
			//if( nTry%3 == 0 ) TotUmMapedBKP2 = TotUmMaped ;

			int cid = MainID * double( BaseGrid->global_cell_number ) / double( NumOfRandTry ) - 1 ;
			//int gcid = BaseGrid->pCell_Cell[ cid ] ;
			int gcid 	=	BaseGrid->GlobalCell_LocalCellNo[ cid ] ;
			int pid 	=	0 ;
			int gpid 	=	0 ;

			for ( int j = 0 ; j < mpi_size ; j++ )
			{
				if ( cid >= CellShift[ j ] )
				{
					pid = j ;
				} else
				{
					//cout <<"myid = " << mpi_id <<"cid = " << cid << ", shift=" <<  CellShift[j] <<endl;
					break ;
				}
			}

			//gpid 	=	BaseGrid->pre_Cell[ cid ].ProcessorNo ;
			gpid	=	BaseGrid->Mesh.Cell_ProcessorNo[ BaseGrid->GlobalCell_MeshCellNo[ cid ] ] ;

			//cout <<"myid = " << mpi_id << ", pid = " << pid << ", gpid = " << gpid << ", cid = " << cid << ", gcid = "<< gcid << endl;

			if ( NumOfRandTry > nTry )
			{
				for ( int i = 0 ; i < InGrid->local_cell_number ; i++ )
				{
					if ( ToCellID[ i ] < 0 )
					{
						ToCellID[ i ]	=	-1 ;
						CpuID[ i ]		=	pid ;
						FromCellID[ i ]	=	cid ;
					}
				}
			}
		}
	}

	if ( mpi_id == 0 ) cout << "nTry: " << nTry << " Total un-mapped: " << NumOfNotFound + NumOfHitBndy ; // << endl;
	if ( mpi_id == 0 ) cout << ", To be mapped:" << NumOfNotFound << "/" << InGrid->global_cell_number ; //<< ", after " << i+1 <<"th mapping";// << endl;
	if ( mpi_id == 0 ) cout << ", Hit Boundary:" << NumOfHitBndy << "/" << InGrid->global_cell_number << endl ; //<< ", after " << i+1 <<"th mapping" << endl;

	//if(mpi_id == 0)(*BaseGrid->debug_out) << "nTry: " << nTry << " Total un-mapped: " << NumOfNotFound + NumOfHitBndy; // << endl;
	//if(mpi_id == 0)(*BaseGrid->debug_out) << ", To be mapped:" << NumOfNotFound << "/" << InGrid->global_cell_number ;//<< ", after " << i+1 <<"th mapping";// << endl;
	//if(mpi_id == 0)(*BaseGrid->debug_out) << ", Hit Boundary:" << NumOfHitBndy << "/" << InGrid->global_cell_number << endl; //<< ", after " << i+1 <<"th mapping" << endl;
}

void CellByCellRayTracking::Tracking(int NumPoint, double *ToPntPosX, double *ToPntPosY, double *ToPntPosZ, double *FrPntPosX, double *FrPntPosY, double *FrPntPosZ, int *ID, int *ToID, int *FrID, int *MPIID)
{
	int 	nTry, ShiftCell, LocalCellID, FromCpu ;
	Cell 	*FromCell, *ToCell ;

	for ( int i = 0 ; i < NumPoint ; i++ )
	{		
		if ( ToID[ i ] == -1 && MPIID[ i ] == mpi_id )
		{
			if ( FrID[ i ] > -1 )
			{
				FromCell 		=	&(BaseGrid->cell[ FrID[ i ] - CellShift[ mpi_id ] ] ) ;
				FrPntPosX[ i ]	=	FromCell->x ;
				FrPntPosY[ i ]	=	FromCell->y ;
				FrPntPosZ[ i ]	=	FromCell->z ;
			} else
			{
				continue ;
			}

			ToCell 	=	NULL ;
			nTry 	=	0 ;

			//    *BaseGrid->debug_out << endl;
			//    *BaseGrid->debug_out << " ToPoint:" << setw(5) << ID[i] << " : " << setw(12) << ToPntPosX[i]<< setw(12) << ToPntPosY[i]<< setw(12) << ToPntPosZ[i] << endl;

			while ( ToCell == NULL )
			{
				ToCell 	=	ToNextCell( FromCell, ToPntPosX[ i ], ToPntPosY[ i ], ToPntPosZ[ i ], &(FrPntPosX[ i ]), &(FrPntPosY[ i ]), &(FrPntPosZ[ i ]) ) ;

				//if ( ToCell != NULL) cout << "TC " <<  ToCell->id << " FC " << FromCell->id << " tx " << ToPntPosX[i] << " ty "<< ToPntPosY[i] << " fx " << FrPntPosX[i] << " fy " << FrPntPosY[i] << endl ;

				if ( ToCell == NULL )
				{
					//-- hit the boundary, next try
					//-- break the while loop
					FrID[ i ]	=	-10000 ;
					ToID[ i ]	=	-10000 ;
					break ;
				} else if ( ToCell != FromCell )
				{
					if ( ToCell->mpi_id == mpi_id )
					{
						FromCell 	=	ToCell ;
						ToCell 		=	NULL ;
						nTry++ ;
					} else
					{
						// put in buffer
						MPIID[ i ]	=	ToCell->mpi_id ;
						FrID[ i ]	=	ToCell->id ;
						break ;
					}
				}

				if ( nTry > 500 )
				{
					/*
					debug << "FromCell:" << setw(5) << FromCell->id << " : ";
					debug << setw(12) << FrPntPosX[i] << setw(12) << FrPntPosY[i] << setw(12) << FrPntPosZ[i] << " C ";
					debug << setw(12) << FromCell->x  << setw(12) << FromCell->y  << setw(12) << FromCell->z << " D ";
					debug << setw(12) << FromCell->x - FrPntPosX[i] << setw(12) << FromCell->y - FrPntPosY[i] << setw(12) << FromCell->z - FrPntPosZ[i];
					debug << endl;
					*/
					
					cout << "Hi ";
				} ;
			} ;

			if ( ToCell == NULL )
			{
				//-- hit the boundary, next try
				//MPIID[i] = MainID;
			} else if ( FromCell == ToCell )
			{
				// cell found
				ToID[ i ]	=	ToCell->id ;
				MPIID[ i ]	=	mpi_id ;
				//cout << "Found " << i << " " << ToID[i] << endl ;
			}
		} else if ( ToID[ i ] == -1 && MPIID[ i ] != mpi_id )
		{

		}
	}

}

void CellByCellRayTracking::CountingBuffer()
{
	int i, j ;
	ReSetOffset() ;
	for ( int i = 0 ; i < iNumPnts ; i++ )
	{
		if ( ToCellID[ i ] == -1 && CpuID[ i ] != mpi_id )
		{
			iNumSend[ CpuID[ i ] ]++ ;
		}
	}

	iNumTotRecv = 0 ;
	iNumTotSend = 0 ;

	for ( i = 0 ; i < mpi_size ; i++)
	{
		j =	iNumSend[ i ] ;
		MPI_Gather( &j, 1, MPI_INT, iNumRecv.get(), 1, MPI_INT, i, BaseGrid->comm ) ;
		iNumTotSend 	+=	iNumSend[ i ] ;
	}
	
	for ( int i = 0 ; i < mpi_size ; i++ )
	{
		SendOffSet[ 1 ][ i + 1 ]	=	SendOffSet[ 1 ][ i ] + iNumSend[ i ] ;
		RecvOffSet[ 1 ][ i + 1 ]	=	RecvOffSet[ 1 ][ i ] + iNumRecv[ i ] ;

		//SendOffSet[ 1 ].push_back ( SendOffSet[ 1 ][ i ] + iNumSend[ i ] ) ;
		//RecvOffSet[ 1 ].push_back ( RecvOffSet[ 1 ][ i ] + iNumRecv[ i ] ) ;
		iNumTotRecv += iNumRecv[ i ] ;
	}

	for ( int j = 2 ; j < iNumDataSize ; j++ )
	{
		for ( int i = 0 ; i < mpi_size + 1 ; i++ )
		{
			//SendOffSet[ j ][ i ]	=	SendOffSet[ 1 ][ i ] * j ;
			//RecvOffSet[ j ][ i ]	=	RecvOffSet[ 1 ][ i ] * j ;
			SendOffSet[ j ].push_back( SendOffSet[ 1 ][ i ] * j ) ;
			RecvOffSet[ j ].push_back( RecvOffSet[ 1 ][ i ] * j ) ;
		}
	}

	SendDataI = new int[ iNumDataSize * iNumTotSend ] ;
	RecvDataI = new int[ iNumDataSize * iNumTotRecv ] ;
	SendDataD = new double[ iNumDataSize * iNumTotSend ] ;
	RecvDataD = new double[ iNumDataSize * iNumTotRecv ] ;
}


void CellByCellRayTracking::SendToBuffer()
{
	MPI_Status status ;
	int iShift, dShift ;

	for ( int i = 0 ; i < mpi_size ; i++ ) iNumSend[ i ] = 0 ;

	for ( int i = 0 ; i < iNumPnts ; i++ )
	{
		if ( ToCellID[ i ] == -1 && CpuID[ i ] != mpi_id )
		{
			SendTag[ i ]	=	1 ;
			iShift =	SendOffSet[ iDataNum ][ CpuID[ i ] ] + iNumSend[ CpuID[ i ] ] * iDataNum ;
			dShift =	SendOffSet[ dDataNum ][ CpuID[ i ] ] + iNumSend[ CpuID[ i ] ] * dDataNum ;

			for ( int j = 0 ; j < 3 ; j++) SendDataD[ dShift + j     ]	=	PntPosDest[ j ][ i ] ;
			for ( int j = 0 ; j < 3 ; j++) SendDataD[ dShift + j + 3 ]	=	PntPosFrom[ j ][ i ] ;

			SendDataI[ iShift + 0 ]	=	PntID[ i ] ;
			SendDataI[ iShift + 1 ]	=	FromCellID[ i ] ;
			SendDataI[ iShift + 2 ]	=	CpuID[ i ] ;
			SendDataI[ iShift + 3 ]	=	ToCellID[ i ] ;

			iNumSend[ CpuID[ i ] ]++ ;
		}
	}

	for ( int i = 0 ; i < mpi_size ; i++ )
	{
		int idest 	= ( mpi_id + i + 1 ) ;
		int isrc 	= ( mpi_id - i - 1 ) ;
		if ( idest >= mpi_size ) idest 	= idest - mpi_size ;
		if ( isrc  < 0         ) isrc 	= isrc  + mpi_size ;
		int itag_in = ( mpi_id + 1 ) * 100000 * ( i + 1 ) ;
		int itag_out = ( idest + 1 ) * 100000 * ( i + 1 ) ;
		MPI_Sendrecv( &SendDataD[ SendOffSet[ dDataNum ][ idest ] ], iNumSend[ idest ] * dDataNum, MPI_DOUBLE, idest, itag_out,
					  &RecvDataD[ RecvOffSet[ dDataNum ][ isrc  ] ], iNumRecv[ isrc ]  * dDataNum, MPI_DOUBLE, isrc,  itag_in,  BaseGrid->comm, &status ) ;
	}

	for ( int i = 0 ; i < mpi_size ; i++ )
	{
		int idest 	= ( mpi_id + i + 1 ) ;
		int isrc 	= ( mpi_id - i - 1 ) ;
		if ( idest >= mpi_size ) idest 	= idest - mpi_size ;
		if ( isrc  < 0         ) isrc 	= isrc  + mpi_size ;
		int itag_in = ( mpi_id + 1 ) * 200000 * ( i + 1 ) ;
		int itag_out = ( idest + 1 ) * 200000 * ( i + 1 ) ;
		MPI_Sendrecv( &SendDataI[ SendOffSet[ iDataNum ][ idest ] ], iNumSend[ idest ] * iDataNum, MPI_INT, idest, itag_out,
					  &RecvDataI[ RecvOffSet[ iDataNum ][ isrc  ] ], iNumRecv[ isrc ]  * iDataNum, MPI_INT, isrc,  itag_in,  BaseGrid->comm, &status ) ;
	}


	for ( int i = 0 ; i < 3 ; i++) PntPosDest_r[ i ] = new double [ iNumTotRecv ] ;
	for ( int i = 0 ; i < 3 ; i++) PntPosFrom_r[ i ] = new double [ iNumTotRecv ] ;

	PntID_r 		=	new int[ iNumTotRecv ] ;
	ToCellID_r 		=	new int[ iNumTotRecv ] ;
	FromCellID_r 	=	new int[ iNumTotRecv ] ;
	CpuID_r 		=	new int[ iNumTotRecv ] ;

	iShift 	=	0 ;
	dShift 	=	0 ;
	for ( int i = 0 ; i < iNumTotRecv ; i++ )
	{
		for ( int j = 0 ; j < 3 ; j++ ) PntPosDest_r[ j ][ i ] = RecvDataD[ dShift + j     ] ;
		for ( int j = 0 ; j < 3 ; j++ ) PntPosFrom_r[ j ][ i ] = RecvDataD[ dShift + j + 3 ] ;
		PntID_r[ i ]		=	RecvDataI[ iShift + 0 ] ;
		FromCellID_r[ i ]	=	RecvDataI[ iShift + 1 ] ;
		CpuID_r[ i ]		=	RecvDataI[ iShift + 2 ] ;
		ToCellID_r[ i ]		=	RecvDataI[ iShift + 3 ] ;
		ToCellID_r[ i ]		=	-1 ;
		iShift +=	iDataNum ;
		dShift +=	dDataNum ;
	}
} ;


void CellByCellRayTracking::ReturnBuffer( int *NumOfNotFound, int *NumOfHitBndy )
{
	MPI_Status status ;
	int NotFound = 0, HitBndy = 0 ;
	int iShift, dShift ;

	/*
	    (*BaseGrid->debug_out) << "Send1 :" << endl ;
	    for ( int i = 0 ; i < mpi_size ; i++ ) (*BaseGrid->debug_out) << setw( 3 ) << i << setw( 6 ) << iNumSend[ i ] << ", " ;
	    (*BaseGrid->debug_out) << endl;
	    for( int i = 0 ; i < iNumTotSend ; i++ )
	    (*BaseGrid->debug_out)
	         << setw(8) << i
	         << setw(8) << SendDataI[ i * iDataNum + 0 ]
	         << setw(8) << SendDataI[ i * iDataNum + 2 ]
	         << setw(8) << SendDataI[ i * iDataNum + 3 ]
	         << endl ;


	    (*BaseGrid->debug_out) << "Get1 :" << iNumTotRecv << endl ;
	    for ( int i = 0 ; i < mpi_size ; i++ ) (*BaseGrid->debug_out) << setw(3) << i << setw(6) << iNumRecv[ i ] << ", " ;
	    (*BaseGrid->debug_out) << endl ;

	    for ( int i = 0 ; i < iNumTotRecv ; i++ )
	    {
	        (*BaseGrid->debug_out)
	             << " ith:" << setw(8) << i
	             << " OID:" << setw(8) << PntID_r[ i ]    << " - " << setw(8) << RecvDataI[ i * iDataNum + 0 ]
	             << " CPU:" << setw(8) << CpuID_r[ i ]    << " - " << setw(8) << RecvDataI[ i * iDataNum + 2 ]
	             << " NID:" << setw(8) << ToCellID_r[ i ] << " - " << setw(8) << RecvDataI[ i * iDataNum + 3 ]
	             << endl;
	    }
	*/

	iShift = 0 ;
	dShift = 0 ;
	for ( int i = 0 ; i < iNumTotRecv ; i++ )
	{
		for ( int j = 0 ; j < 3 ; j++ ) RecvDataD[ dShift + j     ] = PntPosDest_r[ j ][ i ] ;
		for ( int j = 0 ; j < 3 ; j++ ) RecvDataD[ dShift + j + 3 ] = PntPosFrom_r[ j ][ i ] ;
		RecvDataI[ iShift + 0 ]	=	PntID_r[ i ] ;
		RecvDataI[ iShift + 1 ]	=	FromCellID_r[ i ] ;
		RecvDataI[ iShift + 2 ]	=	CpuID_r[ i ] ;
		RecvDataI[ iShift + 3 ]	=	ToCellID_r[ i ] ;
		iShift +=	iDataNum ;
		dShift +=	dDataNum ;
	}

	for ( int i = 0 ; i < mpi_size ; i++ )
	{
		int idest	=	( mpi_id + i + 1 ) ;
		int isrc 	=	( mpi_id - i - 1 ) ;
		if ( idest >= mpi_size ) idest = idest - mpi_size ;
		if ( isrc  < 0         ) isrc = isrc + mpi_size ;
		int itag_in = ( mpi_id + 1 ) * 300000 * ( i + 1 ) ;
		int itag_out = ( idest + 1 ) * 300000 * ( i + 1 ) ;
		MPI_Sendrecv( &RecvDataI[ RecvOffSet[ iDataNum ][ idest ] ], iNumRecv[ idest ] * iDataNum, MPI_INT, idest, itag_out,
					  &SendDataI[ SendOffSet[ iDataNum ][ isrc ] ],  iNumSend[ isrc ]  * iDataNum, MPI_INT, isrc,  itag_in,  BaseGrid->comm, &status ) ;
	}

	for ( int i = 0 ; i < mpi_size ; i++ )
	{
		int idest = ( mpi_id + i + 1 ) ;
		int isrc = ( mpi_id - i - 1 ) ;
		if ( idest >= mpi_size ) idest 	= idest - mpi_size ;
		if ( isrc  < 0         ) isrc 	= isrc  + mpi_size ;
		int itag_in = ( mpi_id + 1 ) * 400000 * ( i + 1 ) ;
		int itag_out = ( idest + 1 ) * 400000 * ( i + 1 ) ;
		MPI_Sendrecv( &RecvDataD[ RecvOffSet[ dDataNum ][ idest ] ], iNumRecv[ idest ] * dDataNum, MPI_DOUBLE, idest, itag_out,
					  &SendDataD[ SendOffSet[ dDataNum ][ isrc ] ],  iNumSend[ isrc ]  * dDataNum, MPI_DOUBLE, isrc,  itag_in,  BaseGrid->comm, &status ) ; 
	}


	for ( int i = 0 ; i < mpi_size ; i++ ) iNumSend[ i ] = 0 ;

	//(*BaseGrid->debug_out) <<"send_back: " << iNumTotSend <<endl;
	int ith = 0 ;
	for ( int i = 0 ; i < iNumPnts ; i++ )
	{
		if ( SendTag[ i ] == 1 )
		{
			iShift 	=	SendOffSet[ iDataNum ][ CpuID[ i ] ] + iNumSend[ CpuID[ i ] ] * iDataNum ;
			dShift 	=	SendOffSet[ dDataNum ][ CpuID[ i ] ] + iNumSend[ CpuID[ i ] ] * dDataNum ;
			iNumSend[ CpuID[ i ] ]++ ;
			/* 
			(*BaseGrid->debug_out)
			 << " ith:" << setw(8)		<< ith << " iShift: "		<< setw(8)	<< iShift
			 << " OID:" << setw(8)		<< PntID[ i ]				<< " - "	<< setw(8) << SendDataI[ iShift + 0 ]
			 //<< "  FC:" << setw(8)	<< FromCellID[ i ]			<< " - "	<< setw(8) << SendDataI[ iShift + 1 ]
			 << " CPU:" << setw(8)		<< CpuID[ i ]				<< " - "	<< setw(8) << SendDataI[ iShift + 2 ]
			 << " NID:" << setw(8)		<< ToCellID[ i ]			<< " - "	<< setw(8) << SendDataI[ iShift + 3 ]
			 << " Num:" << setw(8)		<< iNumSend[ CpuID[ i ] ]	<< endl ;
			 */
			ith++ ;

			if ( PntID[ i ] != SendDataI[ iShift + 0 ] )
			{
				cout << "wrong mapping " << i << endl ;
				//    (*BaseGrid->debug_out) << "w ";
			}

			FromCellID[ i ]			=	SendDataI[ iShift + 1 ] ;
			CpuID[ i ]				=	SendDataI[ iShift + 2 ] ;
			ToCellID[ i ]			=	SendDataI[ iShift + 3 ] ;

			PntPosFrom[ 0 ][ i ]	=	SendDataD[ dShift + 3 ] ;
			PntPosFrom[ 1 ][ i ]	=	SendDataD[ dShift + 4 ] ;
			PntPosFrom[ 2 ][ i ]	=	SendDataD[ dShift + 5 ] ;
		}
	}

	/*
	    for(int i=0; i<iNumTotSend; i++ )
	    {
	            (*BaseGrid->debug_out)
	            << setw(12) << ReturnDataI[0][i]
	            << setw(12) << ReturnDataI[1][i]
	            << setw(12) << ReturnDataI[2][i]
	            << setw(12) << ReturnDataI[3][i] << endl;

	    }

	*/

	for ( int i = 0 ; i < iNumPnts ; i++ )
	{
		if ( ToCellID[ i ] == -1     ) NotFound++ ;
		if ( ToCellID[ i ] == -10000 ) HitBndy++ ;
	}

	//cout << mpi_id << " NotFound " << NotFound  <<endl ;
	MPI_Allreduce( &NotFound, NumOfNotFound, 1, MPI_INT, MPI_SUM, BaseGrid->comm ) ;
	MPI_Allreduce( &HitBndy, NumOfHitBndy, 1, MPI_INT, MPI_SUM, BaseGrid->comm ) ;

	delete [] RecvDataD ;
	delete [] SendDataD ;
	delete [] RecvDataI ;
	delete [] SendDataI ;
	delete [] PntPosDest_r[ 0 ] ;
	delete [] PntPosDest_r[ 1 ] ;
	delete [] PntPosDest_r[ 2 ] ;
	delete [] PntPosFrom_r[ 0 ] ;
	delete [] PntPosFrom_r[ 1 ] ;
	delete [] PntPosFrom_r[ 2 ] ;
	delete [] PntID_r ;
	delete [] CpuID_r ;
	delete [] ToCellID_r ;
	delete [] FromCellID_r ;
}



void CellByCellRayTracking::ReSetOffset()
{
	int i;

	//SendOffSet[ 1 ][ 0 ]	=	0 ;
	//RecvOffSet[ 1 ][ 0 ]	=	0 ;
	
	SendOffSet.reset() ;
	RecvOffSet.reset() ;

	SendOffSet  = boost::shared_array< vector<int> > ( new vector<int> [ iNumDataSize ] ) ;
	RecvOffSet  = boost::shared_array< vector<int> > ( new vector<int> [ iNumDataSize ] ) ;

	SendOffSet[ 1 ].push_back( 0 ) ;
	RecvOffSet[ 1 ].push_back( 0 ) ;

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		//SendOffSet[ 1 ][ i + 1 ]	=	0 ;
		//RecvOffSet[ 1 ][ i + 1 ]	=	0 ;

		SendOffSet[ 1 ].push_back( 0 ) ;
		RecvOffSet[ 1 ].push_back( 0 ) ;

		iNumSend[ i ]	=	0 ;
		iNumRecv[ i ]	=	0 ;
	}

	for ( i = 0 ; i < iNumPnts ; i++ ) SendTag[ i ] = 0 ;
}

Cell* CellByCellRayTracking::ToNextCell( Cell *InCell, double ToPntPosX, double ToPntPosY, double ToPntPosZ, double *FrPntPosX, double *FrPntPosY, double *FrPntPosZ )
{
	int 	i, itag, outtag, outtag2, outtag3 ;
	double 	P1_vec[ 3 ], P1_len, iSign, Normal_vec[ 3 ] ;
	double 	vel_normal[ 6 ], dis_normal[ 6 ], dis_normal_old[ 6 ], time_normal[ 6 ], out_time ;
	Face 	*iFace ;

	itag 		=	0 ;
	outtag 		=	-1 ;
	outtag2 	=	-1 ;
	outtag3 	=	-1 ;
	out_time    =	2.0 ;

	P1_len 		=	0.0 ;

	//P1_vec[ 0 ]	=	ToPntPosX - InCell->x ;
	//P1_vec[ 1 ]	=	ToPntPosY - InCell->y ;
	//P1_vec[ 2 ]	=	ToPntPosZ - InCell->z ;
	//*FrPntPosX	=	InCell->x ;
	//*FrPntPosY	=	InCell->y ;
	//*FrPntPosZ	=	InCell->z ;

	P1_vec[ 0 ]	=	ToPntPosX - *FrPntPosX ;
	P1_vec[ 1 ]	=	ToPntPosY - *FrPntPosY ;
	P1_vec[ 2 ]	=	ToPntPosZ - *FrPntPosZ ;

	P1_len 	+=	P1_vec[ 0 ] * P1_vec[ 0 ] ;
	P1_len 	+=	P1_vec[ 1 ] * P1_vec[ 1 ] ;
	P1_len 	+=	P1_vec[ 2 ] * P1_vec[ 2 ] ;

	P1_len 	=	sqrt( P1_len ) ;

	//debug << "X = " << 		ToPntPosX << ", Y = " << ToPntPosY << endl ;
	//debug << "X = " << (*FrPntPosX) << ", Y = " << (*FrPntPosY) << endl ;
	for ( i = 0 ; i < InCell->face_number ; i++ )
	{
		iFace 	=	InCell->face[ i ] ;
		iSign 	=	InCell->face_sign[ i ] ;

		//dis_normal[ i ]		=	iSign * ( iFace->nx *   ToPntPosX  + iFace->ny *   ToPntPosY  + iFace->nz *   ToPntPosZ  - iFace->nd ) ;
		//dis_normal_old[ i ]	=	iSign * ( iFace->nx * (*FrPntPosX) + iFace->ny * (*FrPntPosY) + iFace->nz * (*FrPntPosZ) - iFace->nd ) ;
		dis_normal[ i ]		=	iSign * ( iFace->nA[ 0 ] *   ToPntPosX  + iFace->nA[ 1 ] *   ToPntPosY  + iFace->nA[ 2 ] *   ToPntPosZ  - iFace->nAd ) ;
		dis_normal_old[ i ]	=	iSign * ( iFace->nA[ 0 ] * (*FrPntPosX) + iFace->nA[ 1 ] * (*FrPntPosY) + iFace->nA[ 2 ] * (*FrPntPosZ) - iFace->nAd ) ;
		vel_normal[ i ]		=	iSign * ( iFace->nA[ 0 ] *  P1_vec[ 0 ] + iFace->nA[ 1 ] * P1_vec[ 1 ]  + iFace->nA[ 2 ] * P1_vec[ 2 ] ) ;
		time_normal[ i ]	=	-dis_normal_old[ i ] / vel_normal[ i ] ;

		if ( dis_normal[ i ] > 0.0 )
		{
			itag++ ;
			outtag3 = i ;
			if ( vel_normal[ i ] > 0.0 )
			{
				if ( out_time > time_normal[ i ] )
				{
					out_time 	=	time_normal[ i ] ;
					outtag 		=	i ;
				}
				outtag2 	=	i ;
			}
		}

		//debug << "\tF = " << i << ", id = " << iFace->id << ", d = " << dis_normal[ i ] << ", v = " <<vel_normal[ i ]  << ", t = " << 	time_normal[ i ] << endl ;
	}
	//debug << " tag = " << outtag << endl ;

	if ( itag == 0 )
	{
		return InCell ;
	} else if ( outtag > -1 )
	{
		if ( InCell->face[ outtag ]->type == 0 )
		{
			//cout << "B " << *FrPntPosX << " p1v "<< P1_vec[0]  << " ot " << out_time << endl ;
			(*FrPntPosX) +=	P1_vec[ 0 ] * out_time * 1.00000001 ;
			(*FrPntPosY) +=	P1_vec[ 1 ] * out_time * 1.00000001 ;
			(*FrPntPosZ) +=	P1_vec[ 2 ] * out_time * 1.00000001 ;
			//cout << "A " << *FrPntPosX << " p1v "<< P1_vec[0]  << " ot " << out_time << endl ;
			return InCell->cell[ outtag ] ;
		} else
		{
			//cout << " hit the boundary";
			return NULL ;
		}
	} else if ( itag == 1 && outtag > -1  )
	{
		return InCell->cell[ outtag ] ;
	}
	cout << " wrong tag in cell by cell tracking "  << itag << " " << outtag << " " << outtag2 << " " << outtag3 << endl ;
	return NULL ;
} ;

void CellByCellRayTracking::InGridToBaseGrid( double *BaseData, double *InData )
{
	MPI_Status status ;

	double 	*SBuffer, *RBuffer ;
	int 	*SCellID, *RCellID ;
	int 	i, j ;

	ReSetOffset() ;
	for ( i = 0 ; i < InGrid->local_cell_number; i++ )
		if ( CpuID[ i ] > -1 ) iNumSend[ CpuID[ i ] ]++ ;

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		j = iNumSend[ i ] ;
		MPI_Gather( &j, 1, MPI_INT, iNumRecv.get(), 1, MPI_INT, i, BaseGrid->comm ) ;
	}

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		SendOffSet[ 1 ][ i + 1 ]	=	SendOffSet[ 1 ][ i ] + iNumSend[ i ] ;
		RecvOffSet[ 1 ][ i + 1 ]	=	RecvOffSet[ 1 ][ i ] + iNumRecv[ i ] ;
	}

	SBuffer 	=	new double[ SendOffSet[ 1 ][ mpi_size ] ] ;
	RBuffer 	=	new double[ RecvOffSet[ 1 ][ mpi_size ] ] ;

	SCellID 	=	new int[ SendOffSet[ 1 ][ mpi_size ] ] ;
	RCellID 	=	new int[ RecvOffSet[ 1 ][ mpi_size ] ] ;

	for ( i = 0 ; i < mpi_size ; i++ ) iNumSend[ i ] = 0 ;

	for ( i = 0 ; i < InGrid->local_cell_number ; i++ )
	{
		if ( CpuID[ i ] > -1 )
		{
			int iShift = SendOffSet[ 1 ][ CpuID[ i ] ] + iNumSend[ CpuID[ i ] ] ;
			SBuffer[ iShift ]	=	InData[ i ] ;
			SCellID[ iShift ]	=	ToCellID[ i ] ;
			iNumSend[ CpuID[ i ] ]++ ;
		}
	}

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		int idest 	=	( mpi_id + i + 1 ) ;
		int isrc 	=	( mpi_id - i - 1 ) ;
		if ( idest >= mpi_size ) idest = idest - mpi_size ;
		if ( isrc  < 0         ) isrc = isrc  + mpi_size ;
		int itag_in = ( mpi_id + 1 ) * 2000000 * ( i + 1 ) ;
		int itag_out = ( idest  + 1 ) * 2000000 * ( i + 1 ) ;
		MPI_Sendrecv( &SBuffer[ SendOffSet[ 1 ][ idest ] ], iNumSend[ idest ], MPI_DOUBLE, idest, itag_out,
					  &RBuffer[ RecvOffSet[ 1 ][ isrc  ] ], iNumRecv[ isrc ],  MPI_DOUBLE, isrc,  itag_in,  BaseGrid->comm, &status ) ;
	}

	for ( i = 0 ; i < mpi_size ; i++ )
	{
		int idest 	=	( mpi_id + i + 1 ) ;
		int isrc 	=	( mpi_id - i - 1 ) ;
		if ( idest >= mpi_size ) idest = idest - mpi_size ;
		if ( isrc  < 0         ) isrc = isrc  + mpi_size ;
		int itag_in 	=	( mpi_id + 1 ) * 3000000 * ( i + 1 ) ;
		int itag_out 	=	( idest  + 1 ) * 3000000 * ( i + 1 ) ;
		MPI_Sendrecv( &SCellID[ SendOffSet[ 1 ][ idest ] ], iNumSend[ idest ], MPI_INT, idest, itag_out,
					  &RCellID[ RecvOffSet[ 1 ][ isrc  ] ], iNumRecv[ isrc ],  MPI_INT, isrc,  itag_in,  BaseGrid->comm, &status ) ;
	}

	for ( i = 0 ; i < RecvOffSet[ 1 ][ mpi_size ] ; i++ )
		BaseData[ RCellID[ i ] - CellShift[ mpi_id ] ]	=	RBuffer[ i ] ;

	delete [] SBuffer ;
	delete [] RBuffer ;
	delete [] SCellID ;
	delete [] RCellID ;
} ;
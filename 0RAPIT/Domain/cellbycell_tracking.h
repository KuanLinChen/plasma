#include "domain.h"


#if !defined(__CellByCell_Tracking_H)
#define __CellByCell_Tracking_H

class CellByCellRayTracking
{

public:

    CellByCellRayTracking( Domain *BaseGrid ) ;
    ~CellByCellRayTracking() ;

    ofstream debug ;

    Domain *BaseGrid, *InGrid ;

    int mpi_size, mpi_id, MainID ;

    boost::shared_array< vector<int> > SendOffSet, RecvOffSet ;

    int *ReturnOffset, *GetOffset ;
    boost::shared_array<int> iNumSend, iNumRecv, CellShift ;
    int iNumTotSend, iNumTotRecv ;

    int iDataNum, dDataNum ;
    int iNumDataSize ;
    int iNumPnts ;

    int    *SendDataI, *RecvDataI ;
    double *SendDataD, *RecvDataD ;

    int NumOfRandTry, nTry;
    boost::shared_array<double> PntPosDest[ 3 ], PntPosFrom[ 3 ] ;
    boost::shared_array<int>    PntID, ToCellID, FromCellID, CpuID, SendTag ;
    double *PntPosDest_r[ 3 ], *PntPosFrom_r[ 3 ] ;
    int    *PntID_r, *ToCellID_r, *FromCellID_r, *CpuID_r;

    void PointsMapping( int iNumPnt, double *Xoc, double *Yoc, double *Zoc, int *CellIDMap ) ;

    void ReSetOffset() ;
    void Tracking( int NumPoint, double *PntPosXt, double *PntPosYt, double *PntPosZt, double *PntPosXf, double *PntPosYf, double *PntPosZf, int *ID, int *CellID, int *OCellID, int *MPIID ) ;
    void CountingBuffer() ;
    void SendToBuffer() ;
    void ReturnBuffer( int *NotFound, int *NumOfHitBndy ) ;
    void TryMap( int itag ) ;
    void InGridToBaseGrid( double *BaseData, double *InData ) ;

    CellByCellRayTracking & operator << ( Domain *Grid2 ) ;
    CellByCellRayTracking & operator << ( int iTry ) ;

    Cell* ToNextCell( Cell *InCell, double Xto, double Yto, double Zto, double *Xfr, double *Yfr, double *Zfr ) ;
} ;


#endif

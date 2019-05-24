#include <string>
#include "preprocessing.h"

#if !defined(__PREPROCESSINGTOOL_H)
#define __PREPROCESSINGTOOL_H

#define FormatFIELDVIEW 0
#define FormatCGNS 1
#define FormatMSH 2
#define FormatPTI 3

int Get_filetype ( string _n ) ;
bool Compare_string ( string s1, string s2 ) ;
bool Find_string ( string s1, string s2 ) ;

void CellNodeCCW( int NodeNum, double *Node_x, double *Node_y, int *NodeNo ) ;

double CalculateAngleDeg( double a, double b, double c ) ;

template <class TYPE>
void BubbleSort( TYPE *a, int *Index, int Num ) ;

template <class TYPE>
void Swap( TYPE *a, TYPE *b ) ;

// For Cell-Cell relation
bool CheckNeighbor( int TypeA, int TypeB, vector<int> *_NodeA, vector<int> *_NodeB, int SurfaceA, int SurfaceB, CellMapping *pMapping ) ;

// For FaceInformation: InnerFace
bool CheckNeighbor( int TypeA, int TypeB, vector<int> *_NodeA, vector<int> *_NodeB, int SurfaceA, int SurfaceB, vector<int> *_Face_Node, CellMapping *pMapping ) ;
bool CheckNeighbor( int TypeA, vector<int> *_NodeA, int SurfaceA, vector<int> *_NodeB, CellMapping *pMapping ) ;

// For FaceInformation: BCFace
bool CheckNeighbor( int TypeA, vector<int> *_NodeA, int SurfaceA, vector<int> *_NodeB, vector<int> *_Face_Node, CellMapping *pMapping ) ;

void CreateCellCellRelation( int NodeNum, int CellNum, int *InnerFaceNum, int *InterFaceNum, int *FaceNum, vector<int> *Node_Cell, int *Cell_Form, int *Cell_Type, vector<int> *Cell_Node, vector<int> *Cell_Cell ) ;

void CalculateCellVolume( string geometry, int CellNum, int *Cell_Form, double *min_cell_length, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, vector<int> *Cell_Node, double *Cell_Volume, CellMapping *pMapping ) ;



//void CreateCellCellRelation( int dimension, int NodeNum, int CellNum, int *InnerFaceNum, int *InterFaceNum, vector<int> *Node_Cell, int *Cell_Type, vector<int> *Cell_Node, vector<int> *Cell_Cell ) ;

#endif

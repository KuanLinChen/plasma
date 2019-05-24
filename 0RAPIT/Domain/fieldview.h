#include <string>
#include <vector>

using namespace std ; 

#if !defined(__FIELDVIEW_H)
#define __FIELDVIEW_H

void GetNumbers_FieldView ( int dimension, int *NodeNum, int *CellNum, int *BCNum, int *BCFaceNum, string Filename ) ;

void ReadGrid_FieldView ( int dimension, double scale, int NodeNum, int CellNum, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int *BCFaceNum, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, int *Cell_Form, int *Cell_Type, string *Cell_Typename, int *Cell_FaceNum, vector<int> *Cell_Node, vector<int> *Node_Cell, vector<int> *Node_BCFace, string Filename ) ;

#endif
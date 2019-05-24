#include <string>
#include <vector>


using namespace std ; 

#if !defined(__PTI_H)
#define __PTI_H

void GetNumbers_PTI( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, int *MeshType_Num, string Filename ) ;
void ReadGrid_PTI ( int dimension, double scale, int NodeNum, int *MeshType_Num, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int *Cell_Form, int *Cell_FaceNum, int *Cell_Type, string *Cell_Typename, vector<int> *Cell_Node, vector<int> *Node_Cell, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, vector<int> *Node_BCFace, string Filename ) ;

#endif
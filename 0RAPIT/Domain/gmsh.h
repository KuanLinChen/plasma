#include <string>
#include <vector>


using namespace std ; 

#if !defined(__GMSH_H)
#define __GMSH_H

void GetNumbers_GMSH( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, string Filename ) ;
void ReadGrid_GMSH ( int dimension, double scale, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int *Cell_Form, int *Cell_FaceNum, int *Cell_Type, string *Cell_Typename, vector<int> *Cell_Node, vector<int> *Node_Cell, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, vector<int> *Node_BCFace, string Filename ) ;

#endif
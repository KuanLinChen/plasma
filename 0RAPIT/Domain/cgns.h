#include <string>
#include <vector>

using namespace std; 

#if !defined(__CGNS_H)
#define __CGNS_H

void GetNumbers_CGNS( int dimension, int *NodeNum, int *CellNum, int *BC_TypeNum, int *BCFaceNum, int *MeshType_Num, string Filename ) ;
//void ReadGrid_CGNS ( int dimension, double scale, int NodeNum, int CellNum, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BCFaceNum, string Filename ) ;
void ReadGrid_CGNS ( int dimension, double scale, int NodeNum, int CellNum, int *MeshType_Num, double *Node_x, double *Node_y, double *Node_z, double *Cell_x, double *Cell_y, double *Cell_z, int *BC_TypeNum, int *BC_TypeNum_forFace, vector<string> *BC_Typename, int BCFaceNum, int *BCFace_Type, string *BCFace_Typename, vector<int> *BCFace_Node, int *Cell_Form, int *Cell_Type, string *Cell_Typename, int *Cell_FaceNum, vector<int> *Cell_Node, vector<int> *Node_Cell, vector<int> *Node_BCFace, string Filename ) ;

/*
void read_cgns_uns( ) ;
void read_coordinates(char *name, double *x, double *y, double *z, int *con, int **isize) ;
void read_bc_conditions(char *name, int **bc) ;
int get_bc_num(char *name) ;
int get_bc_ind_num(char *name, int ib) ;
*/

#endif

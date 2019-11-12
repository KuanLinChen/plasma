#pragma once
#include "ultraMPP.h"
#include "petscsys.h" 
using namespace std;

extern int gargc2 ;
extern char **gargv2 ;

extern int mpi_rank, mpi_size ;
extern int nDim ;
extern ultraMPP plasma ;
extern map<string,int> var_name ; 

extern json &json_bc_setting    ;//= *(Surface_charge_test.get_json_input_parameter("boundary_setting") );
extern json &json_cell_setting  ;//= *(Surface_charge_test.get_json_input_parameter("volume_setting") );
extern map<string,double> cell_parameter ;
extern map<string,double> face_parameter ;
extern map<string,int>    MPP_face_tag ;
extern map<string,int>    MPP_cell_tag ;

extern map<int,int> face_type ;
extern map<int,int> cell_type ;

extern map<int, string>	type_typename ;
extern map<string, int>	typename_type ;

const double		vacuum_permittivity 	=  8.8541878176e-12 ; /* Unit in F/m */

#define MASTER_NODE		0

#define BULK			  0
#define PLASMA			9958
//#define SOLID			9958

#define Debug_MatCoeff_DD_Zero 1

#define ZERO			1.E-9


#define POWER					100
#define GROUND 				200
#define DIELECTRIC 		300

#define NEUMANN				400
#define SYMMETRIC		  400

#define SOLID_POWER		500
#define SOLID_GROUND 	600


/*--- Equation Type ---*/
#define ELECTRON	 0
#define ION				 1
#define NEUTRAL		 2
#define BACKGROUND 3
#define POISSON 	 4


#define INLET		1110
#define OUTLET		2220
#define EULER_WALL	3330
#define FAR_FIELD	4440
#define EMITTER 	5550
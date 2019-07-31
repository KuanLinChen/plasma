#pragma once
#include "ultraMPP.h"
#include "petscsys.h" 
using namespace std;

extern int mpi_rank, mpi_size ;
extern int nDim ;
extern ultraMPP plasma ;
extern map<string,int> var_name ; 

const double		vacuum_permittivity 	=  8.8541878176e-12 ; /* Unit in F/m */

#define MASTER_NODE		0

#define BULK			0
#define PLASMA			0
#define SOLID			9958

#define Debug_MatCoeff_DD_Zero 1

#define ZERO			1.E-9

#define ELECTRODE_PLASMA_INTERFACE		10
#define ELECTRODE_DIELECTRIC_INTERFACE	20
#define DIELEC_PLASMA_INTERFACE			30

#define POWER_DIELEC_INTERFACE	13
#define GROUND_DIELEC_INTERFACE	23


#define POWER			100
#define POWER_1			101
#define POWER_2			102
#define POWER_3			103
#define POWER_4			104
#define POWER_5			105
#define POWER_6			106
#define POWER_7			107
#define POWER_8			108
#define POWER_9			109

#define GROUND 			200
#define GROUND_1		201
#define GROUND_2		202
#define GROUND_3		203
#define GROUND_4		204
#define GROUND_5		205
#define GROUND_6		206
#define GROUND_7		207
#define GROUND_8		208
#define GROUND_9		209

#define DIELECTRIC 		300
#define DIELECTRIC_1	301
#define DIELECTRIC_2	302
#define DIELECTRIC_3	303
#define DIELECTRIC_4	304
#define DIELECTRIC_5	305
#define DIELECTRIC_6	306
#define DIELECTRIC_7	307
#define DIELECTRIC_8	308
#define DIELECTRIC_9	309

/*--- Equation Type ---*/
#define ELECTRON	0
#define ION			1
#define NEUTRAL		2
#define BACKGROUND	3
#define POISSON 	4

#define NEUMANN		400
#define SYMMETRIC	400

#define INLET		1110
#define OUTLET		2220
#define EULER_WALL	3330
#define FAR_FIELD	4440
#define EMITTER 	5550
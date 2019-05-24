#include <iostream>
#include "domain.h"


#if !defined(__TECPLOT_OUTPUT_H)
#define __TECPLOT_OUTPUT_H

void tecplot_output ( string , int, string *, double **, string , Domain *, bool) ;
void tecplot_output1D ( string , int, string *, double **, string , Domain *, bool) ;
void tecplot_output_root ( string , int, string *, double **, string , Domain *, bool) ;

//void tecplot_output_3D ( string , int, string *, double **, string , Domain *, bool) ;

#endif

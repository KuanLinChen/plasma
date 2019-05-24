#include <iostream>
#include "domain.h"


#if !defined(__VTK_OUTPUT_H)
#define __VTK_OUTPUT_H

//void tecplot_output ( string , int, string *, double **, string , Domain *, bool) ;
//void tecplot_output_root ( string , int, string *, double **, string , Domain *, bool) ;
void vtk_output ( string , int, string *, double **, string , Domain *, bool) ;
void vtk_output_root ( string , int, string *, double **, string , Domain *, bool) ;

//void tecplot_output_3D ( string , int, string *, double **, string , Domain *, bool) ;

#endif

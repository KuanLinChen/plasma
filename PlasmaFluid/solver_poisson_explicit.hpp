#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <metis.h>
#include <algorithm>
#include <fstream>
#include <set>
#include "main.h"
#include "PFM.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "PETSc_solver.h"

#include "domain.h"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"
#include "io.h"
#include "scalar.h"

using namespace std;

class CEPoisson 
{


	public:
		CEPoisson();
		Scalar Res[3], Source[3], u_exact, p_exact ;
		void exact( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  ) ;
	/*--- Solver Control Parameter ---*/	
		int Correction, its ;  /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 
		void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> & ) ;
		void CalculateSourceTerm( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
		void ComputeResidue( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable ) ;
		void Solve( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &  ) ;
		void Plot( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  ) ;

};

#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
//#include <metis.h>
#include <algorithm>
#include <fstream>
//#include <set>
//#include "main.h"
#include "PFM.hpp"
#include "json.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

//#include "PETSc_solver.h"

//#include "domain.h"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"
#include "Table.hpp"

//#include "scalar.h"

using namespace std;

class CFDMaxwell
{


	public:
		
		CFDMaxwell();
		
//		CTable CollTable ;
		/*!
		 * \brief module initialization.
		 */	
		void Init( boost::shared_ptr<CConfig> &config ,boost::shared_ptr<CVariable> &var ) ;
		int its ;  /*!< \brief ksp iteration number */ 

		/*!
		 * \brief Compute the effective permittivity (electron only).
		 */	
		void SOLVE                                  ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;

		void UltraMPPComputeCurrentDenAndSourceTerm( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
		
		void UltraMPPComputePowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
		
		void UltraMPPComputeInstantPowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
};

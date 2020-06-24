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
		
		CTable FVFD_CollTable ;
		/*!
		 * \brief module initialization.
		 */	
		void Init( boost::shared_ptr<CConfig> &config ,boost::shared_ptr<CVariable> &var ) ;
		int its ;  /*!< \brief ksp iteration number */ 


		//Solve compled equation.	 
		void SOLVE                                  ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;

		//Compute right hand side, omega * mu0 * J.
		void UltraMPPComputeCurrentDenAndSourceTerm( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
		
		// Compute EM_power absorptiom [Wm^-3] in time average form.	
		void UltraMPPComputePowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
		
		// Compute EM_power absorptiom [Wm^-3] which P = Re{J}*Re{E}.		
		void UltraMPPComputeInstantPowerAbsorptionFromMaxwell( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
		
		// Compute EM_power & ES_power, which integral EM_power / ES_power absorption over all computation domain.
		void UltraMPPComputeTotalPower( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> &var ) ;
};

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

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

//#include "PETSc_solver.h"

//#include "domain.h"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"

//#include "scalar.h"

using namespace std;

class CPoisson 
{


	public:
		
		CPoisson();
		/*!
		 * \brief module initialization.
		 */	
		void Init( boost::shared_ptr<CConfig> &config ) ;
		int its ;  /*!< \brief ksp iteration number */ 

		/*!
		 * \brief Compute the net charged density sum_j( sgn(q)*e*n_j ). Note: e: elementary charge.
		 */	
		void UltraMPPComputeNetCharge               ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

		/*!
		 * \brief Compute the  permittivity =  effective permittivity.
		 */	
		void UltraMPPComputePermitt                 ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

		/*!
		 * \brief Compute the effective permittivity. Semi-Implicit accroding K. M. Lin et al. CPC 183 (2012) 1225–1236.
		 */	
		void UltraMPPComputeEffectivePermitt        ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

		/*!
		 * \brief Compute the effective permittivity (electron only).
		 */	
		void UltraMPPComputeEffectivePermittEleOnly ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

		/*!
		 * \brief Compute the surface charge for plasma-dielectric interface.
		 */	
		void UltraMPPComputeSurfaceCharge           ( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

		/*!
		 * \brief Compute the effective permittivity (electron only).
		 */	
		void SOLVE                                  ( boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

		/*!
		 * \brief Compute the displacement current.
		 */	
		void UltraMPPComputeDispCurrentDensity( boost::shared_ptr<CVariable> &var ) ;

		/*!
		 * \brief Update the electrical component map. 
		 * \param[in] var    - Variable.
		 */	
		void UpdateElectricalMap( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		double SineVoltage( string FaceType, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

	/* TEST FUNCTION */
	// void MatA_SourceB( boost::shared_ptr<CConfig> &config,boost::shared_ptr<CVariable> &var ) ;
	// void ComputeGradient( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	// void SOLVE_TEST( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  ) ;
};

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
#include "json.hpp"
#include "main.h"
#include "PFM.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "PETSc_solver.h"

#include "domain.h"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"
#include "scalar.h"

using namespace std;
class CEnergyDensity 
{


	public:
	/*! 
	 * \brief Constructor of the class.
	 */		
	CEnergyDensity();

	/*!
	 * \brief Allocate the solution variables with the configuration setup.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] index  - Index of this module, it should be match the species index.
	 */ 
	void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &config, int index ) ;


	int iSpecies ; /*!< \brief Species index for this module */ 
	/*--- PETSc Solver ---*/	
		PETScSolver s ;
		int *d_nnz, *o_nnz ;
		int row, col[ 5 ], ncol ;
		double C[ 5 ] ;
		double *NormalizeCoeff ;
	/*--- Solver Control Parameter ---*/	
	double GVarN[2], GVarP[2] ;
	double C53, C43, Reflec ;
	int Correction, its ; /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 
	bool fixTe, TG, Insert ;
	int eLOSS, WallType ;
	void Solver( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Calculate the matrix A and source term B w/o the non-orthogonal correction. 
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Bulid_A_B_1st_default( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Calculate the matrix A and source term B w/ the non-orthogonal correction. 
	 *		  Note: The S-G flux treat as convection flux, maybe don't need to correct the non-orthogonal term.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Bulid_A_B_2nd( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Reset the number density gradient to be zero.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate number density gradient using LSQ.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate temperature.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void CalculateTemperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	double DotProduct(double *A, double *B ){
		return A[0]*B[0] + A[1]*B[1] ;
	};
};

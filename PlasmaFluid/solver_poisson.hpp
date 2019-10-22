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

	/*--- PETSc Solver ---*/	
		map<string, int> BCs ;
		CScalarFace face_data ;
		void ultraMPP( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		void Calculate_NetCharge_minus( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		//PETScSolver s ;
		//int *d_nnz, *o_nnz ;
		//int row, col[ 5 ], ncol ;
		//double C[ 5 ] ;
	/*--- Solver Control Parameter ---*/	
		double GVarN[2], GVarP[2] ;
		int Correction, its ;  /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 
		void Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	/*!
	 * \brief Solver procedure.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */		
	void Solve( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	//void Solve0( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  ) ;
	/*!
	 * \brief Calculate the matrix A w/ the non-orthogonal correction and source term B.  
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Bulid_A_B_0th( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Bulid_A_B_Orthogonal( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Bulid_A_B_Orthogonal2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	/*!
	 * \brief Calculate the matrix A and source term B w/o the non-orthogonal correction. 
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Bulid_A_B_1st( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Calculate the matrix A and source term B w/ the non-orthogonal correction. 
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */	
	void Bulid_A_B_2nd( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	void ScalePermitt( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void ScalePermittBack( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Scale_NetCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Scale_NetChargeBack( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	/*!
	 * \brief Update the electrical component map. 
	 * \param[in] var    - Variable.
	 */	
	void UpdateElectricalMap( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;



	void Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate the gradient(Φ). Note: E = -∇Φ, so you need minus the value to obtain electric fields. 
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */		
	void Calculate_Gradient_Neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Calculate_Gradient_Neumann2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate the net density sum_j( sgn(q)*e*n_j ). Note: e: elementary charge.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */	
	void Calculate_NetCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Calculate_NetCharge_2N_minus_N( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

		void CalculateEffectivePermitt( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		void CalculateEffectivePermittEleOnly( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		void CalculatePermitt( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

		void CalculateDispCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		//void CalculateTotalCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		double SineVoltage( string FaceType, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		double DotProduct(double *A, double *B ){
			return A[0]*B[0] + A[1]*B[1] ;
		};
};

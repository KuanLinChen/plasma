#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <algorithm>
//#include <fstream>
//#include <set>
#include "json.hpp"
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
class CDriftDiffusion 
{


	public:
	/*! 
	 * \brief Constructor of the class.
	 */		
	CDriftDiffusion();

	/*!
	 * \brief Allocate the solution variables with the configuration setup.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] index  - Index of this module, it should be match the species index.
	 */ 
	void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &config, int index ) ;

	int iMatrix ;
	int iSpecies, SpeciesType ; /*!< \brief Species index for this module */ 
	/*--- PETSc Solver ---*/	
	/*--- Solver Control Parameter ---*/	
		double GVarN[2], GVarP[2] ;
		int Correction, WallType, its ;  /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 
		
	/*!
	 * \brief Solver procedure.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Solve( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Solver procedure for diffusion term.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Solve_Diffusion( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Calculate the matrix A and source term B w/o the non-orthogonal correction. Note: This is for reduce the code. 
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void Bulid_A_B_1st_default( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Bulid_A_B_1st_neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Bulid_A_B_1st_0D( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	//void Bulid_A_B_1st( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void CalculateSurfaceCharge_test( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Semi_Empirical_Temperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

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
	 * \brief Calculate number density gradient using LSQ. Note: Using simple shapefunction in the boundary.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void Calculate_Gradient_Neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate reconstruction flux @ cell centered according to Raja JCP 228 (2009) 4435–4443.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void CalculateAvgDDFlux_default( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void CalculateAvgDDFlux_neumann( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void CalculateAvgDDFlux_zero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

	void CalculateDDConvection( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void CalculateCondCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		//void CalculateTotalCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate the surface charge on dielectric interface for charge species.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void CalculateSurfaceCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;

		double DotProduct(double *A, double *B ){
			return A[0]*B[0] + A[1]*B[1] ;
		};

		double Bernoulli_Coeff( double sgn_X, double sgn_q, double mu, double D, double dx, double E_x)
		{
	        double X;
	        double B;

	        if ( mu <= 0.0 || sgn_q == 0 ) X = 0.0;
	        else X = (double) (-1.0) * sgn_X * sgn_q * mu * E_x * dx / D ;

	        if ( X != 0.0 && fabs(E_x) > 1e-4) B = Bernoulli_func(X);
	        //if ( X != .0 ) B = Bernoulli_func(X);
	        else B = 1.0;

	        return B;
		};
		 double Bernoulli_func( double x )
		{
		    if ( fabs(x) < 1e-9 ){
		        return 1.0;
		    }else{
		        return x/(exp(x) - 1.0 );
		    }
		}
	bool fixTe, TG, Insert ;
	double Reflec ;


};

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
class CFluidModel 
{


	public:
	/*! 
	 * \brief Constructor of the class.
	 */		
	CFluidModel();

	/*!
	 * \brief Allocate the solution variables with the configuration setup.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] index  - Index of this module, it should be match the species index.
	 */ 
	void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &config, int index ) ;


	boost::shared_array <double> Flux ;/*!< \brief Convection Flux[iVar] @ cell interface. */ 

	double **Res ;/*!< \brief Residue for each equation Res[iEqn][iCell] */ 
	int iSpecies, SpeciesType ; /*!< \brief Species index for this module */ 

	Scalar Pressure ; /*!< \brief Pressure of i ion species */ 
	Scalar Thermal2 ; /*!< \brief Thermal velocity square. */
	Scalar *CollisionFreq ;
	//Scalar Kappa ;
	
	double CrossSection ;
	bool Grave ;
	/*--- PETSc Solver ---*/	
		PETScSolver s ;
		int *d_nnz, *o_nnz ;
		int row, col[ 5 ], ncol ;
		double C[ 5 ] ;
		double Omaga ;
	/*--- Solver Control Parameter ---*/	
		double GVarN[2], GVarP[2] ;
		double *Mx, *My, *Mz, *CollisionIntegral ;
		int Correction, its ;  /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 

	/*!
	 * \brief Continuity solver procedure.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void Solve_Continuity( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Calculate the convection flux and adding to resudue.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void ComputeFlux_HLL( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Adding source/sink to residue and marching to next time level.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void ContinuityIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  ) ;

	/*!
	 * \brief Momentum solver procedure.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void Solve_Momentum( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Adding force, pressure and collision term to residue and marching to next time level.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void MomentumIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &variable  ) ;

	void EnergyDensityIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  ) ;

	/*!
	 * \brief Calculate the ion temperature from total energy.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void CalculateTemperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	/*!
	 * \brief Calculate the ion total energy according to the initial consonant ion temperature.
	 *		  Note: This function only call when ion temperature is assume be constant.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */
	void CalculateTotalEnergy( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void Calculate_Gradient_T( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate the surface charge on dielectric interface for charge species.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */		
	void CalculateSurfaceCharge( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;


		// void Bulid_A_B_2nd( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/*!
	 * \brief Calculate number density gradient using LSQ.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;
		
	void CalculateCondCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		//void CalculateTotalCurrentDensity( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void CalculateKappa( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  ) ;
	void CalculateIonNeutralCollisionFrequency( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void CalculateThermal2( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
	void CalculateCollisionIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  ) ;
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
	bool fixTe ;


};

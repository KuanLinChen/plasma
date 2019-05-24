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
//#include "variable_structure.hpp"
#include "variable_structure_NS.hpp"
#include "post_structure.hpp"
#include "config_structure.hpp"
#include "scalar.h"

using namespace std;
class CNavierStokes 
{


	public:
	/*! 
	 * \brief Constructor of the class.
	 */		
	CNavierStokes();

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

	double Mass ;
	/*--- PETSc Solver ---*/	
		//PETScSolver s ;
		//int *d_nnz, *o_nnz ;
		//int row, col[ 5 ], ncol ;
		//double C[ 5 ] ;
		double Omaga ;
	/*--- Solver Control Parameter ---*/	
		double GVarN[2], GVarP[2] ;
		int Correction, its ;  /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 
	double MaxERR[ 5 ], **buffer ;
	/*!
	 * \brief Continuity solver procedure.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void Solve( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariableNS> &, boost::shared_ptr<CPost> & ) ;

	/*!
	 * \brief Calculate the convection flux and adding to resudue.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void ComputeFlux_HLL( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariableNS> & ) ;

	/*!
	 * \brief Adding source/sink to residue and marching to next time level.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void ContinuityIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &variable  ) ;
	void CalculateLocalTimeStep( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var  );

	/*!
	 * \brief Adding force, pressure and collision term to residue and marching to next time level.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */ 
	void MomentumIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &variable  ) ;

	void EnergyIntegral( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariableNS> &var  ) ;

	/*!
	 * \brief Calculate the ion temperature from total energy.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] var    - Variable.
	 */


		// void Bulid_A_B_2nd( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariableNS> & ) ;
	void Zero_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariableNS> &var ) ;

	/*!
	 * \brief Calculate number density gradient using LSQ.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void Calculate_Gradient( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariableNS> &var ) ;
		
		double DotProduct(double *A, double *B ){
			return A[0]*B[0] + A[1]*B[1] ;
		};
};

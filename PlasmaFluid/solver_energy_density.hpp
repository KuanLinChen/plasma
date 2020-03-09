#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <fstream>
#include <set>
#include "PFM.hpp"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"

using namespace std;
class CEnergyDensity 
{


	public:
	/*! 
	 * \brief Constructor of the class.
	 */		
	CEnergyDensity();
	ultraMPP energy_density ;
	/*!
	 * \brief Allocate the solution variables with the configuration setup.
	 * \param[in] domain - Grid information.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] index  - Index of this module, it should be match the species index.
	 */ 
	void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &config, int index ) ;
	
	/*--- For Debug ---*/	
	int ncol ;
	double C[ 5 ] ;

	int iSpecies ; /*!< \brief Species index for this module */ 

	/*--- Solver Control Parameter ---*/	
	double GVarN[2], GVarP[2] ;
	double C53, C43, Reflec ;
	int its ; /*!< \brief expilicit correction of the non-orthogonal effect & ksp iteration number */ 
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
	void Bulid_A_B_1st_GEC( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Bulid_A_B_1st_BBC( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
					//BBC = Brezmes & Breitkopf, 2015 ; COMSOL, 2013
	void Bulid_A_B_1st_Hagelaar( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Bulid_A_B_1st_Hagelaar_Txy( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;
	void Bulid_A_B_1st_zero( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> &, boost::shared_ptr<CVariable> & ) ;

	/*!
	 * \brief Calculate temperature.
	 * \param[in] domain - Grid information.
	 * \param[in] var    - Variable.
	 */
	void CalculateTemperature( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	/**
	 * @brief     Calculate the poewe absorption. (JdotE times volume)
	 *
	 * @param     m     The mesh module.
	 * @param     var   The variable module.
	 * Note: return valule to variable module->PowerAbs.
	 */
	void CalculatePowerAbs( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CVariable> &var ) ;

	double DotProduct(double *A, double *B ){
		return A[0]*B[0] + A[1]*B[1] ;
	};
};

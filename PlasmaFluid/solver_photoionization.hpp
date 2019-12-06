#pragma once 
#include "iostream"
#include "PFM.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"

using namespace std;

class CHelmholtz
{
	public:
		CHelmholtz() ;


		ultraMPP helmholtz ; /*!< \brief ultraMPP PDE solver object. */
		variable_set photoionization ;/*!< \brief potential */

		double *RHS ;  /*!< \brief The RHS of the helmholtz equation. */
		int Tag_RHS, Tag_SphSum ;

		void Init(boost::shared_ptr<CConfig> &config) ;

		map<string,double> helmholtz_cell_parameter ;
		map<string,double> helmholtz_face_parameter ;
		map<string,int>    helmholtz_face_tag ;
		map<string,int>    helmholtz_cell_tag ;

		int numFitTerm ;     /*!< \brief How many terms to use to fit the curve */
		double *A, *lambda ; /*!< \brief Fitting coefficient. */
		double *SphSum ;     /*!< \brief Summation of photo ionization source term. */
		double O2PartialPressure ; /*!< \brief The partial pressure of oxygen gas. */
		void SOLVE( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var  ) ;
};

#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <metis.h>
#include "config_structure.hpp"
#include "domain_structure.hpp"
#include "chemistry.h"
#include "scalar.h"
#include "Table.hpp"
#include "PFM.hpp"

using namespace std;
/*! 
 * \class CVariable
 * \brief Contain the solution variables and chemistry module.
 */
class CVariableNS
{
	public:
	/*! 
	 * \brief Constructor of the class.
	 */	
	CVariableNS();

	/*!
	 * \brief Allocate the solution variables with the configuration setup.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> & ) ;

	void Calculate_LSQ_Coeff_Scalar( boost::shared_ptr<CDomain> &m ) ;
	void InitialCondition( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, int iSpecies ) ;
	void UpdateSolution( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	double DotProduct(double *A, double *B ){
		//cout<<"A"<<endl;
		return A[0]*B[0] + A[1]*B[1] ;
	};
	double PN[ 3 ], Pf[ 3 ], Nf[ 3 ], PPf[ 3 ], NPf[ 3 ], fPf[ 3 ], mf[ 3 ] ;
	Scalar LSQ_Cx[6], LSQ_Cy[6], LSQ_Cz[6] ;/*!< \brief Lease-Square Coefficient using scale. LSq_C_dir[iFace][cell]. */
	Scalar U0, U1, U2, U3, U4, Dt ;
	Scalar PreU0, PreU1, PreU2, PreU3, PreU4, MachNumber ;

};
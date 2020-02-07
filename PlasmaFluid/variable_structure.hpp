#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "config_structure.hpp"
#include "domain_structure.hpp"
#include "chemistry.h"
#include "Table.hpp"
#include "PFM.hpp"


using namespace std;
class CScalarFace
{
	public:
	double *data ;
	int data_id ;
	string data_name ;
	CScalarFace(){};
	void initial(string name)
	{
	  data_name = name ;
	  data_id = plasma.set_face_data(  &data, data_name ) ;
	}
	double & operator[] ( int i ) { return data[i] ;}
	const double & operator[] ( int i ) const { return data[i] ;}

	/*---*/
	int get_data_id()
	{
	        return data_id ;
	}
	void zero()
	{
		for ( int i = 0 ; i < plasma.Mesh.face_number ; i++ )
			data[ i ] = 0.0  ;
	}
	/*---*/
	string get_data_name()
	{
		return data_name ;
	}
	/*---*/
	inline CScalarFace & operator= ( const double value ) {
		for ( int i = 0 ; i < plasma.Mesh.face_number ; i++ ) {
			data[ i ] = value  ;
		}
	}
};
class CScalar
{
	public:
	double *data ;
	int data_id ;
	string data_name ;
	CScalar(){};
	void initial(string name)
	{
	        data_name = name ;
	        data_id = plasma.set_parallel_cell_data(  &data, data_name ) ;
	}
	inline double & operator[] ( int i ) { return data[i] ;}
	inline const double & operator[] ( int i ) const { return data[i] ;}

	/*---*/
	inline int get_data_id()
	{
	        return data_id ;
	}
	void zero()
	{
		for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ )
			data[ i ] = 0.0  ;
	}
	/*---*/
	inline string get_data_name()
	{
		return data_name ;
	}
	/*---*/
	inline CScalar & operator= ( const double value ) {
		for ( int i = 0 ; i < plasma.Mesh.cell_number ; i++ ) {
			data[ i ] = value  ;
		}
	}
	/*---*/
	inline CScalar & operator= ( const CScalar & c )
	{
	        plasma.syn_parallel_cell_data( data_id ) ;
	}
};
/*! 
 * \class CVariable
 * \brief Contain the solution variables and chemistry module.
 */
class CVariable
{
	public:
	/*! 
	 * \brief Constructor of the class.
	 */	
	CVariable();

	/*!
	 * \brief Allocate the solution variables with the configuration setup.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void Init( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> & ) ;

	/*! 
	 * \brief Bulid the cell perperties.
	 * \param[in] domain - Name of the file with the grid information.
	 */
	void CellProperties	( boost::shared_ptr<CDomain> & ) ;

	/*! 
	 * \brief Initial the solution variable array.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void InitialConditions( boost::shared_ptr<CDomain> &, boost::shared_ptr<CConfig> & ) ;

	/*! 
	 * \brief Calculate total background gas pressure.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CalTotalPressure( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;

	/*! 
	 * \brief Save the solution variable for next time step.
	 * \param[in] domain - Name of the file with the grid information.
	 */
	void UpdateSolution ( boost::shared_ptr<CDomain> & ) ;

	/*! 
	 * \brief Calculate reduced electric field: E/N (unit: Td = 10^-21 V*m^2)
	 * \param[in] domain - Name of the file with the grid information.
	 */
	void UltraMPPComputeReducedElectricField() ;

	/*! 
	 * \brief Use reduced electric field to lookup the electron energy density tabe.
	 * \param[in] domain - Name of the file with the grid information.
	 */
	void CalMeanEnergy( boost::shared_ptr<CDomain> & ) ;
	

	int SolutionFieldNum ;

	int NetQ_Tag ;

	/*--- UltraMPP variables ---*/
	void UltraMPPVarInit() ;

	map<string,int> VarTag ;

	variable_set Potential ;/*!< \brief potential */
	double *ChargeDen,      /*!< \brief net charge density */ //-rho/eps//
		   *RealChargeDen,      /*!< \brief net charge density */ //rho [C/m^3] //
				 *eps,            /*!< \brief Material permittivity */ 
				 *eps_eff ;       /*!< \brief Effective permittivity for semi-implicit poissiony */ 
	double *Ex, *PreEx,     /*!< \brief current & previous Electric fields in X-dir. */ 
	       *Ey, *PreEy,     /*!< \brief current & previous Electric fields in Y-dir. */ 
	       *Ez, *PreEz,     /*!< \brief current & previous Electric fields in Z-dir. */ 
				 *Etd,            /*!< \brief reduce Electric fields in unit: Td */ 
				 *Emag ;          /*!< \brief Electric fields magnitude */ 
	double *Kappa ;				  /*!< \brief the variable for modified sound speed. See My IEEE paper. */
	double *plot_var, *avg_plot_var ;

	void UltraMPPAvgVarInit() ;
	double *AvgPotential,   /*!< \brief cycle-averaged potential */
			*AvgEx,          /*!< \brief cycle-averaged electric fields in X-dir. */ 
			*AvgEy,          /*!< \brief cycle-averaged electric fields in Y-dir. */ 
			*AvgEz ;         /*!< \brief cycle-averaged electric fields in Z-dir. */ 

	void UltraMPPInitialCellParameter() ;

	//---------Variable definition of ICP simulation ---------------------------
	double Coil_frequency  ; // ICP coil freqency
	double omega  ; // ICP coil angular freqency
	double Coil_Current  ;
	variable_set E_phi_Re, E_phi_Im ;
	double *eps_FVFD ;
	double collision_frequency ; // mom_transfer_rate
	double *sigma_p_Re_plasma 	; // Real 		part of plasma conductivity in plasma mesh 
	double *sigma_p_Im_plasma 	; // Imaginary 	part of plasma conductivity in plasma mesh
	double *sigma_p_Re_FVFD 	; // Real 		part of plasma conductivity in FVFD   mesh 
	double *sigma_p_Im_FVFD 	; // Imaginary 	part of plasma conductivity in FVFD   mesh 
	double *k_square_Re	; // Real 		part of square of helmholtz constant
	double *k_square_Im	; // Imaginary 	part of square of helmholtz constant
	double *CurrentDen	;
	double *Re_eq_source; // Source term in Real 		part FD maxwell equation
	double *Im_eq_source; // Source term in Imaginary 	part FD maxwell equation
	double *Power_Absorption_plasma ;
	double *Power_Absorption_FVFD ;
	//--------------------------------------------------------------------------








	CScalar *DD_Convection ;/*!< \brief Drift-Diffusion Approximation convection term (for semi-implicit poisson's eqnuation). */ 

	CScalar TotalNumberDensity ;
	CScalar MPI_ID ;
	CScalar *Debug ;
	CScalar Beta ;/*!< \brief For Correct Ion sound Speed */ 
	double Ramp_factor ;
	CScalar  *T,  *AvgT,  *PreT ;/*!< \brief Present, cycle-averaged and Previoud time temperature. */ 
	CScalar *U0, *AvgU0, *PreU0 ;/*!< \brief Present, cycle-averaged and Previoud time number density. */ 
	CScalar *U1, *AvgU1, *PreU1 ;/*!< \brief Present, cycle-averaged and Previoud time X-mass flux */ 
	CScalar *U2, *AvgU2, *PreU2 ;/*!< \brief Present, cycle-averaged and Previoud time Y-mass flux */ 
	CScalar *U3, *AvgU3, *PreU3 ;/*!< \brief Present, cycle-averaged and Previoud time Z-mass flux */ 
	CScalar *U4, *AvgU4, *PreU4 ;/*!< \brief Present, cycle-averaged and Previoud time energy flux */ 
	CScalar *JouleHeating, *AvgJouleHeating ; /*!< \brief Present, cycle-averaged Joule heation term */


	CScalar *LFASourceSink ; /*!< \brief  */
	//CScalar Kappa ;
	//CScalar Force_x, Force_y ;/*!< \brief  Force in x & y direction. */
	CScalar ** GradT ;
	CScalar **GradU0 ;/*!< \brief Number density gradient [x,y,z]. */ 
	CScalar **GradU4 ;/*!< \brief Energy density gradient [x,y,z].*/ 

	CScalar *ProductionRate ;
	CScalar *Mobi, *Diff ;/*!< \brief Mobility and diffusivity */ 
		
	CScalar TotalGasPressure ;/*!< \brief Only condider the background gas. */ 

	CScalar DebyeLength, AvgDebyeLength ;/*!< \brief DebyeLength. */ 
	CScalar CFL, AvgCFL ;/*!< \brief DebyeLength. */ 


	CScalar Eps,/*!< \brief Effective permittivity for semi-implicit poissiony */ 
		   Eps0 ;/*!< \brief Material permittivity */ 




	CScalar NetQ ; /*!< \brief Net charged density for poisson source term */ 

	/*--- Chemistry Module ---*/
	chemistry Chemistry ;
	double PowerAbs;
	/*! 
	 * \brief Update the source/sink, electron energy source and transport coefficients.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void ChemistryUpdate( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> & ) ;

	double 	*InputNumberDensity, /*!< \brief Input number density point for YM chemistry module */ 
			*InputTemperature, /*!< \brief Input temperature point for YM chemistry module */ 
			*InputEField, /*!< \brief shoule be delete */
			*EnergySourcePoint, 
			*coll_freq ;
	
	//double*           EnergySourcePoint ;
    vector<double*>   DiffusivityPoint;
    vector<double*>   MobilityPoint;
    vector<double*>   ReactionRatePoint;

    /*! 
	 * \brief Update the electron transport coefficient (mobility and diffusivity).
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
    void UpdateElectronTransport( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config  ) ;

    /*! 
	 * \brief Update the ion and neutral transport coefficient (mobility and diffusivity).
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
    void UpdateIonNeutralTransport( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
    CTable CollTable ;/*!< \brief Collision table for electron transport coefficient calculation */ 
    CTable MeanEnergyTable ;/*!< \brief reduce electric field vs mean energy */ 
    double K2eV ;
    map< int, CTable> MobilityMap ;/*!< \brief Mobility Table */
    map< int, CTable> DiffusivityMap ;/*!< \brief Diffusion Table */
    CTable AlphaTable, EtaTable ;/*!< \brief Diffusion Table */
    double Qe, PI, Me, Kb ;

    CTable eEnergyLossTable ;
    CScalar eEnergyLoss, eAvgEnergyLoss ;
    double PhysicalTime ;/*!< \brief Physical time. */
    double Dt ;/*!< \brief Time step size. */
    /*! 
	 * \brief Reset the average variable to be zero.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
    void ResetAvgZero( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;

    /*! 
	 * \brief Adding each instantaneous variabe at each time to average variable.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
    void AddAverage  ( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;



	double AvgPowerAbs ;/*!< \brief Cycle averaged Power absorption. */
	void ResetAvgZero_PowerAbs( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	void AddAverage_PowerAbs ( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;






  void CalculateElectrodeCurrent( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	void ResetAvgZero_Electrode( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	void AddAverage_Electrode( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	CScalar   **CondJD ;/*!< \brief Conduction current density. */ 
	CScalar  DispJD[3] ;/*!< \brief Displacement current density[x,y,z]. */ 
	CScalar TotalJD[3] ;/*!< \brief Total Current density[x,y,z]. */ 

    double I_PowerElectrode_local, I_GroundElectrode_local ;

    double *CondI_PowerElectrode_local, *CondI_GroundElectrode_local ;

    double DispI_PowerElectrode_local, DispI_GroundElectrode_local ;

    double I_PowerElectrode_global_sum, I_GroundElectrode_global_sum ;

    double I_AvgPowerElectrode, I_AvgGroundElectrode ;

    double *CondI_PowerElectrode_global_sum, *CondI_GroundElectrode_global_sum ;

    double Disp_PowerElectrode_global_sum, Disp_GroundElectrode_global_sum ;

    double Volt ;
    
    void Alpha_Beta( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
    //void ComputeDebyeLengthRatio_CFL( boost::shared_ptr<CDomain> &m ) ;
		//void ComputeDebyeLengthRatio_CFL( boost::shared_ptr<CDomain> &m,  boost::shared_ptr<CConfig> &config  ) ;
	
    /*! 
	 * \brief Calculate Least-square coeff store in scalar.
	 * \param[in] domain - Name of the file with the grid information.
	 */ 
	void Calculate_LSQ_Coeff_Scalar( boost::shared_ptr<CDomain> &m ) ;
	double PN[ 3 ], Pf[ 3 ], Nf[ 3 ], PPf[ 3 ], NPf[ 3 ], fPf[ 3 ], mf[ 3 ] ;
	CScalar Global_Id, Local_Id ;
    CScalar LSQ_Cx[6], LSQ_Cy[6], LSQ_Cz[6] ;/*!< \brief Lease-Square Coefficient using scale. LSq_C_dir[iFace][cell]. */
    double DotProduct(double *A, double *B ){
		//cout<<"A"<<endl;
		return A[0]*B[0] + A[1]*B[1] ;
	};
	CScalar Energy_Term[6], Momentum_Term[3] ;

    /*! 
	 * \brief Compute source/sink for two-fluid (electron and positive ion )system.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SourceSink_2Fluid( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;

    /*! 
	 * \brief Compute source/sink for Three-fluid (electron, positive and negertive ion ) system.
	 * \param[in] domain - Name of the file with the grid information.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SourceSink_3Fluid( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;

	void SourceSink_PSST_2018( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	

	//	void CalculateEnergyLossFromTable( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
	void UltraMPPComputeEnergyLossFromTable( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;
    /*! 
	 * \brief Initialize reference valut.
	 */
	double Ref_L, Ref_N, Ref_SS, Ref_ES, Ref_Qe, Ref_Mass, Ref_Mu, Ref_Diff, Ref_Eps, Ref_Phi, Ref_V, Ref_Te, Ref_t, Ref_Flux, Ref_EN, Ref_Rho, Ref_JD, Ref_SQ, Ref_Kb, Ref_EField ;
	void ReferenceValue_Init( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config ) ;



	void UltraMPPComputeArgonIonTemperaturePHELPS( int iSpecies ) ;
	// CScalar dPhi, dEx, dEy ;
	// CScalar dPrePhi, dPreEx, dPreEy ;

};

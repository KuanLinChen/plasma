#pragma once
#include <cmath>
#include <vector>
#include <fstream>
#include <map>
#include "json.hpp"
#include "PFM.hpp"

using nlohmann::json;
using namespace std;
/*! 
 * \class CSpecies
 * \brief Species information class.
 */
class CSpecies 
{
	public:
	/*! 
	 * \brief Constructor of the class.
	 */	
	CSpecies(){} ;

	string Name ; /*!< \brief Species name for output and Scalar name */

	int Index,/*!< \brief Species index */
		Type,/*!< \brief Species type: 0: for electron, 1: for ion, 2: for neutral, 3: for background gas */
		MobilityType, /*!< \brief mobility type: 0: for constant, 1: for collision table */
		DiffusivityType ; /*!< \brief diffusivity type: 0: for constant, 1: for collision table */

	double Charge,/*!< \brief Species charge, electron: -X.0, positive ion: x.0, negative ion: -x.0, neutral: 0.0  */
		   sgn,/*!< \brief Species charge sign, electron: -1.0, positive ion: 1.0, negative ion: 1.0, neutral: 0.0  */ 
		   Mass_Kg,/*!< \brief Species mass unit kg */
		   InitialDensity,/*!< \brief Iniitial number density unit m^-3 */ 
		   InitialTemperature,/*!< \brief Iniitial temperature unit eV */
		   InitialPressure ;/*!< \brief Initial pressure unit torr */

	double MobilityValue,/*!< \brief Constant mobility value unit XXXXX */
		   DiffusivityValue ;/*!< \brief Constant mobility value unit XXXXX */

	string MobilityFile, DiffusivityFile, CollisionFreqFile ;

	bool ConstantMobilityUpdate,/*!< \brief For speedup, initially will be false */
		 ConstantDiffusivityUpdate ;/*!< \brief For speedup, initially will be false */

	double  LJ_Potential, /*!< \brief Lennard-Jones potential */
			LJ_Parameter ;/*!< \brief Lennard-Jones Parameter unit Angstroms */
	double 	Gamma ;/*!< \brief Cp/Cv, for monatomic is 5/3, for diatomic is 1.4 */
	bool 	Activate ;
};
/*!
 * \class CConfig
 * \brief Main class for defining the problem; basically this class reads the Species.inp and Solver.inp file, and
 *        stores all the information.
 */
class CElectrical
{
	public:
	/*! 
	 * \brief Constructor of the class.
	 */	
	CElectrical(){} ;
	string    Name ;
	int    Type ;
	double Frequency   ;/*!< \brief Apply frequency unit: Hz */
	double Period ;/*!< \brief AC period unit: sec. */
	double Voltage_p2p ;/*!< \brief Apply Voltage peak to peak unit: volt */
	double BiasVoltage ;/*!< \brief Bias voltage unit: volt */
	string VoltageFileName ;/*!< \brief apply voltage table: time v.s voltage */
	double BC_Voltage ;
	double RampFactor,T_min, T_max, RampTime, InitialRampFactot ;

	void UpdateRampFactor( double PhysicalTime ) {
		if( PhysicalTime > RampTime or RampTime < 0.0 ) RampFactor = 1.0 ;
		else RampFactor = (T_max - T_min)/RampTime*PhysicalTime+InitialRampFactot ;
	}
	void UpdateVoltage( double PhysicalTime ) {
		double PI = 4.0*atan(1.0);
		double time=0.0 ;
		double AC_Amplitude = 0.5*Voltage_p2p ;

		if ( PhysicalTime > Period ) {

			int n_period = (int)( PhysicalTime/Period ) ;
			time = PhysicalTime - (double)( n_period*Period ) ;

		} else {

			time = PhysicalTime ;
	
		}

		double duration = 2.0 * PI * time / Period ;

		BC_Voltage = AC_Amplitude*sin(duration) ;
		//BC_Voltage *= RampFactor ;

		if( fabs(BC_Voltage) < ZERO ) BC_Voltage = 0.0 ;
		BC_Voltage += BiasVoltage ;
	}
}; 
class CEquation
{
	public:
	/*! 
	 * \brief Constructor of the class.
	 */	
	CEquation(){} ;
	int Equation 	 ;
	int WallBoundaryType ;/*!< \brief Boundary condition type, 0:default, 1: neumann, 2: zero */
	/*! 
	 * \brief Initialization module.
	 */	
	void Init(){
		Equation   = -999 ;
		WallBoundaryType = 0 ;
	};
};
class CConfig 
{
	public:
	/*! 
	 * \brief Constructor of the class.
	 */	
	CConfig() ;

	/*! 
	 * \brief Initialization module. (read "Species.inp" and "solver.inp" file)
	 */	
	void Init( string ) ;

	/*--- Solver Information ---*/
	/*! 
	 * \brief Read the solver control information input file.
	 * \param[in] filename - The defult solver input file name is "Solver.inp".
	 */	
	void ReadSolverFile ( string ) ;

	string PFM_Assumption,/*!< \brief Plasma fluid model assumption: "LFA" for local field approximation, "LMEA" for local mean energy approximation*/
				 PFM_SubModel,/*!< \brief Submodel: "Beouf", "Nishida", "Cathode" ...etc */
				 eMeanEnergyFile,
				 eEnergyLossFile,/*!< \brief elastic/inelastic power loss file name */
				 CasePath ;/*!< \brief Case input file path. */
			
	bool Normalize ; 
			
	double BGMass_Kg, BGTemperature_K ;/*!< \brief  Background  mass (most havey one) and background temperature */

	CEquation Equation[5] ;/*!< \brief Equations informaiton vector, 0-Electron, 1-Ion, 2-Neutral, 3-Background, 4-Poisson's */
	double SecondaryElectronEmissionCoeff ;
	double SecondaryElectronEmissionEnergy ;
	double ReflectionCoefficient ;

	/*--- ElectricalComponent ---*/
	int ElectricalComponent ;/*!< \brief Number of Electrical Component (eg, nunber of Power or ground electrode ) */
	vector<CElectrical> Electrical ;/*!< \brief Electrical informaiton vector */
	map< int, CElectrical> ElectricalMap ;/*!< \brief Electrical informaiton mpa */
	double MaxFrequency ;/*!< \brief Maximum frequency for calculate the time step size. */
	double DC_TimeStep ;/*!< \brief Time step size for DC simulation. */
	int StepPerCycle ;/*!< \brief How many time step per one rf-cycle using max frequency. */
	//double Dt ;/*!< \brief Time step size. */
	int SimulationCycles ;/*!< \brief How many cycles in this simuation . */
	int Cycle ;
	int ExitCycle ; /*!< \brief Exit Cycle for this simulation . */

	int WRT_Insta_Freq ; /*!< \brief Instantaneous data output frequency per each rf-cycle. */
	int WRT_Cycle_Freq ; /*!< \brief Cycle average data output frequency. */

	int MON_Insta_Freq ; /*!< \brief ksp iteration monitors frequency per each rf-cycle. */
	int MON_Cycle_Freq ; /*!< \brief ksp iteration monitors frequency. */

	bool Average_Switch ;

	/*--- Species Information ---*/
	vector<CSpecies> Species ;/*!< \brief Specise informaiton vector */

	/*! 
	 * \brief Read the species information input file.
	 * \param[in] filename - The defult species input file name is "Species.inp".
	 */	
	void ReadSpeciesFile( string ) ;

	int ElelectronNum,/*!< \brief Number of electron species. (Only 1) */
		IonNum,/*!< \brief Number of ion species. */
		NeutralNum,/*!< \brief Number of neutral species. */ 
		BackGroundNum,/*!< \brief Number of Background species. */ 
		TotalSpeciesNum, /*!< \brief Number of total species. */ 
		ChargeSpeciesNum ;/*!< \brief Number of charge species. */ 
	int OutputFormat ;/*!< \brief IO_TECPLOT: 0, IO_VTK: 1, IO_TECPLOT: 2 */ 
	
	double P_back, T_back ;


};

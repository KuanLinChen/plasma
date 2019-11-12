/*************************************************************************************************
**************************************************************************************************
**     Title : Chemistry module                                                                 **
**     Work  : 1. To estimate the properties and parameter of reactions like source-sink,       **
**                electron energy loss due to collision, collisional frequence.                 **
**             2. Plus calculation of diffusion and mobility of each species in Feb. mid term.  **
**     Date  : 2007-12-19 to 2008-05-12                                                         **
**************************************************************************************************
*************************************************************************************************/
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <mpi.h>

using namespace std;

#ifndef chemistry_H
#define chemistry_H

#define FORMULA   1
#define TABLE 	  2
#define CONSTANT  3
#define GAS_TABLE 4

//inline double interpolation(double*,double*,double*,double*,double*);

class  GLOBAL
{   
	public: 
		vector<double>   global ;  

};


class  table_opt
{
	public:
	int      size_n;
	double   dT;
	double   *T;
	double   *value;
	double   max_T;
	double   min_T;
	table_opt(){};
	~table_opt(){};
};



/**
 * @brief Class for species identity. 
 *        The species_identity contains a lot of information about species(e, ions, neutrals).
 *        For example: names, mass, diffusion, mobility, source-sink term, etc.
 */
class species_identity
{
	public:
	vector<string>    name ;       /*!< \brief Name of Species.  */
	int               e_id ;       /*!< \brief electron species id. default is 0.  */
	vector<int>	      neutral_id ; /*!< \brief neutral species id . */
	vector<int>       ion_id ;     /*!< \brief ion species id.  */

	vector<GLOBAL>    sourcesink ; /*!< \brief source/sink of Species.  */
	vector<GLOBAL>    diffusion ;  /*!< \brief diffusion coefficient of Species.  */
	vector<GLOBAL>    mobility ;   /*!< \brief mobility coefficient of Species.  */

	/**
	 * @brief      Reads species infomations.
	 *
	 * @param[in]  <unnamed>  Sepcies file name
	 */
	void              ReadSpeciesInfomations( string ) ;

	vector<int>      	diff_type;
	vector<int>       mob_type;
	vector<double>   	diff_constant;
	vector<double>    mob_constant;
	vector<double>    mass;
	vector<double>    amu;
	vector<double>    self_diff_constatn;
	vector<double>    self_mob_constatn;
	vector<double>    binary_diameter;
	vector<double>    LJ;
	vector<double>    polarizability;
	vector<double>    viscosity;
	vector<double>    reduce_amu;
	vector<string>    mob_file;
	vector<string>    diff_file;
	vector<int>       interpolation;
	vector<table_opt> thermal_table;
	vector<string>    thermal_file;
	
	
//  Plus information of light radiation  2008-12-29	
	vector<string>                light_name;
	vector<double>                wavelength;
	vector<GLOBAL>                light_power;
		
	species_identity(){ }
	~species_identity(){ }
};


/*********************************************************************************************************************
**********************************************************************************************************************
        |--------------------------------------------------------------------------------------|
        | The single_reaction contains each reactions message and parameter                    |
	| For example: reaction rate,                                                          |
	| It is used by Object chemistry.  		                                       |
        |--------------------------------------------------------------------------------------|
**********************************************************************************************************************
**********************************************************************************************************************/
class  single_reaction
{
	public:
	string          filename;
	string          output_name;
	int             rate_type;                //1: formula 2: table 3: constant
	int             is_data_output;           //0: no 1: yes
	int             is_interpolation;         //0: no 1: yes
	double          threshold;
	double          rate_constant;
	vector<int>     source;
	vector<int>     sink;
	vector<double>  formula_coefficient;
	vector<double>  origin_eV;
	vector<double>  origin_rate;
	vector<double>  Channel_SourceSink;
//  Plus information of light radiation  2008-12-29	
	int             emit_id;
//  Additions at 2009-01-07 optimize source and sink term	
	vector<int>     source_ID;
	vector<int>     sink_ID;
//  Additions at 2009-02-03 record momentum reaction
	int             momentum_reactions;
	single_reaction(){}
	~single_reaction(){}
};
 

class  single_sourcesink
{
	public:
	int   sink_n;
	int   source_n;
	int   power_is_one;
	int   IS_momentum_Xsfer;
	int   *reactant_ID;
	int   *product_ID;
	int   *coefficient_reactant;
	int   *coefficient_product;
	single_sourcesink(){}
	~single_sourcesink(){}

};

class   table_new
{
	public:
	//static   int     n;
	static   int     channel_n;
	static   double  dT;
	double   T;
	double   *R;
	double   *dR;
	double   diff;
	double   ddiff;
	double   mob;
	double   dmob;
	
	table_new(){};
	~table_new(){};
};
//int table_new::n = 0;
//double table_new::dT = 0.0;

class  KNN
{
	public:
	int     sign;
	int     k;
	int     *nn;
	int     nn_size;
	int     partial_c;
	
	KNN(){};
	~KNN(){};
};
class partial_s
{
	public:
	int  n;
	KNN  *knn;
	partial_s(){};
	~partial_s(){};
};

class  apply_channel
{
	public:
	int   *sign;
	int   n;
	int   *join;
};
/*********************************************************************************************************************
**********************************************************************************************************************
    |--------------------------------------------------------------------------------------|
    | The species_identity contains Object single_reaction,Object species_identity.        |
	| It basic connect to applied program code.                                            |
	| Autorun function: rate_formula(), lookup_table(), read_reaction_table(),             |
	|                   loading_each_channel_file(), newtable_establish().                              |
	| Called by user: sourcesink() -                                                       |
	|                 diffusion_mobility_estimate() -                                      |
	|                 ptr_*() - return variable pointer to user                            |
    |--------------------------------------------------------------------------------------|
**********************************************************************************************************************
**********************************************************************************************************************/
class  chemistry
{
	public:
	species_identity          species;
	vector<single_reaction>   channel;
	GLOBAL                    energy_loss;
	GLOBAL	                  coll_frequency ;
	GLOBAL                    CollisionFrequency ;        
	int                       domain_size;
	string                    SpeciesFileName;
	string                    ChannelFileName;
	int                       Te_index;
	double                    T_Gas;
	
	double                    rate_formula(double*,int);
	void                      read_input( string ) ;
	void                      energy_loss_coll_fre_diff_mob_init();
	void                      read_reaction_table(string);
	void                      loading_each_channel_file();
	void                      newtable_establish();
	void                      CalSourceSinkRate_for_Te(const double *);
	void                      CalSourceSinkRate_for_gas(const double *);
	void                      CalTotalPressure(const double *, const double *);
	void                      SourceSink(const double *,const double *);
	vector<double*>           ptr_source_sink();
	double*                   ptr_coll_frequency();
	double*                   ptr_energy_loss();
	void                      UpdateElectronDiffMob(const double* ,const double* );
	void                      UpdateNeutralDiffMob(const double* ,const double* );
	void                      UpdateIonDiffMob(const double* ,const double* ,const double * );
	vector<double*>           ptr_diffusion();
	vector<double*>           ptr_mobility();
	vector<double*>           ptr_total_ChannelSourceSink();
	double*                   ptr_single_ChannelSourceSink(int);
	double*                   ptr_thermal(const double*);
	double*                   series_thermal;
	double*                   total_gas_density;
	double*                   total_gas_pressure;
//***********************************************************
// Plus at 2009-01-02                                       *
/**/	vector<double*>           ptr_light_power();      //*
/**/	int                       LightSize();            //*
//***********************************************************
//***************************************************************
// Plus at 2009-01-02  for optimization of source-sink          *
// Pre-calculate all not necessary parameter                    *
/**/    void                      OptimizationSourceSink();   //*
//***************************************************************
	
//***************************************************************
// Plus at 2009-04-01  for optimization of source-sink          *
// We figure out call vector::size() cost a lot of time         *	
/**/    int                       species_size;               //*
/**/    int                       channel_size;               //*
/**/    int                       ion_size;                   //*
/**/    int                       neutral_size;               //*
/**/    int                       sourcesink_size;            //*
/**/    int                       light_size;                 //*
/**/    single_sourcesink         *OPT_channel_sourcesink;     //*
//***************************************************************

//*****************************************************************
// Plus at 2009-04-01  for optimization of channel                *
// Plus at 2009-04-21  for partial D,mu / partial Te              *
// We figure out call vector::size() cost a lot of time           *
// channel_species: Which species act in channel                  *
// chem_rate: rates for sourcesink in sreies memory               *
/**/    int                       *channel_species;             //*
/**/    table_new                 *channel_new_table;           //*
/**/    int                       *which_is_table;              //*
/**/    int                       *channel_to_table;            //*
/**/    double                    *chem_rate_buf;                //
/**/    double                    *diff_buf;                      //
/**/    double                    *mob_buf;                      //
/**/    double                    *ddiff_buf;                      //
/**/    double                    *dmob_buf;                      //
/**/    int                       table_e_size;                    //
/**/    int                       table_ion_size;                 //
/**/    int                       table_neutral_size;             //
/**/    double                    table_Te_min;                   //
/**/    double                    table_Te_max;		         //
/**/    double                    table_gas_min;		 //
/**/    double                    table_gas_max;		 //
/**/    inline double*            search_table(const double *);  //
/**/    inline double*            search_perturb(const double *);//
/**/    inline double*            search_diffusion(const double *);//
/**/    inline double*            search_mobility(const double *);//
/**/    inline  double*            search_ddiff(const double *);
/**/    inline  double*            search_dmob (const double *);
/**/    apply_channel             *net_sourcesink;               //
//******************************************************************

//*******************************************************************
//    Plus at 2009-04-02  for Formfunction                          *
//    partial sourcesink / partial n_i                              *
        partial_s                 *for_Formfunction;              //*
/**/    double                    PS_PN(int,int,int,double*);     //*
/**/    apply_channel              *record_total_ss;               //
//*******************************************************************

//*******************************************************************
//    Plus at 2009-04-15  for Formfunction                          *
//    partial electron_energy / partial n_i                         *
/**/    double                    PE_PN(int,int,double*,double*); //*
/**/	double*                   e_loss_energy;                  //*
/**/    int*                      e_channel_ID;                   //*
/**/    partial_s                 *e_loss_for_Formfunction;        //*
//*******************************************************************



//*******************************************************************
//    Plus at 2009-04-14  for Formfunction                          *
//    partial sourcesink / partial n_i                              *
/**/    int                       ShowChemRate(const double*);     //
/**/    int                       ShowChemTable(const double*);    //
//*******************************************************************

//*******************************************************************
//    Plus at 2009-04-16  for Formfunction                          *
//    Te perturbation                                               *
/**/    double                    perturb_Te;                      //

/**/    double                    *chem_rate_perturb_buf;
/**/    int                       set_Te_perturbation(double);     //
/**/    apply_channel              e_loss_dTe_ID;                  //
/**/    apply_channel              *S_dTe_ID;                      //
/**/    double                    PS_PTe(int,int,double*);         //
/**/    double                    PE_PTe(int,double*,double*);     //
//*******************************************************************

	double            PD_PTe( int );
	double            PMu_PTe(int );

	//void                      init(string  ,string ,int );
	void                      init(string  ,int );
	void                      init(string  ,string);
	int                       SpeciesID(string);
	
	//int                       Record_n_gas_rate_table;
	
	chemistry(string  ,string ,int );
	chemistry(string  ,string );
	chemistry(){}
	~chemistry() { }
};

/*********************************************************************************************************************
**********************************************************************************************************************
        |--------------------------------------------------------------------------------------|
        | Xe diffusion and mobility object.                                                    |
        |--------------------------------------------------------------------------------------|
**********************************************************************************************************************
**********************************************************************************************************************/

class noble_drift
{
	public:
	double     mu_plus;
	double     C;
	double     k_plus;
	double     D;
	double     mobility_Ward(double,double);
	double     v_fun_1(double);
	double     v_fun_2(double);
	double     v_plus;
	//double     diffusion_McDaniel(double,double,double,double,double);
	
	noble_drift(){}
	~noble_drift(){}
	
};






#endif

#include "post_structure.hpp"
using namespace std ;
CPost::CPost()
{
}
void CPost::Init()
{
}
void CPost::OutputFlow( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int Cycle, int Step )
{

	plasma.set_output( "0Flow-"+to_string(Cycle)+"-"+to_string(Step) ) ;

	plasma.set_output( var->Potential.tag_current ) ;
	//plasma.set_output( var->VarTag["permittivity"] ) ;
	//plasma.set_output( var->VarTag["effective_permittivity"] ) ;
	plasma.set_output( var->VarTag["ChargeDen"             ] ) ;

	plasma.set_output( var->VarTag["Ex"] ) ;
	plasma.set_output( var->VarTag["Ey"] ) ;
	if ( nDim == 3 )
	plasma.set_output( var->VarTag["Ez"] ) ;

	plasma.set_output( var->VarTag["Emag"] ) ;


	/* Temperature */
	if ( config->PFM_SubModel != "Streamer" ) {
		for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
		plasma.set_output( var->T[ iSpecies ].data_id ) ;
	}//End streamer

	/* number density */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->U0[ iSpecies ].data_id ) ;
	


	if ( config->PFM_SubModel != "Streamer" ) {
		/* x- flux */
		for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
		plasma.set_output( var->U1[ iSpecies ].data_id ) ;
		/* y- flux */
		for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
		plasma.set_output( var->U2[ iSpecies ].data_id ) ;

		for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
		plasma.set_output( var->U4[ iSpecies ].data_id ) ;
	}//End streamer

	/* Production Rate */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->ProductionRate[ iSpecies ].data_id ) ;
	
	/* Joule Heating */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->JouleHeating[ iSpecies ].data_id ) ;

	/* Mobility */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->Mobi[ iSpecies ].data_id ) ;
	

	/* Diffusion coefficients */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->Diff[ iSpecies ].data_id ) ;
	plasma.write_output(to_string(Step)) ;

}
void CPost::OutputAverageFlow( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int Cycle )
{
	plasma.set_output( "1AvgFlow-"+to_string(Cycle) ) ;

	plasma.set_output( var->VarTag["AvgPotential"] ) ;

	plasma.set_output( var->VarTag["AvgEx"] ) ;
	plasma.set_output( var->VarTag["AvgEy"] ) ;
	if ( nDim == 3 )
	plasma.set_output( var->VarTag["AvgEz"] ) ;


	plasma.set_output( var->eAvgEnergyLoss.data_id ) ;
	
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var-> AvgT[ iSpecies ].data_id ) ;
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var->AvgU0[ iSpecies ].data_id ) ;
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var->AvgU1[ iSpecies ].data_id ) ;
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var->AvgU2[ iSpecies ].data_id ) ;
	if( nDim == 3 ) 
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var->AvgU3[ iSpecies ].data_id ) ;
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var->AvgU4[ iSpecies ].data_id ) ;
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) plasma.set_output( var->AvgJouleHeating[ iSpecies ].data_id ) ;
	plasma.write_output(to_string(Cycle)) ;
}

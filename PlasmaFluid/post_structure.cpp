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

	plasma.set_output( var->Phi.data_id ) ;
	plasma.set_output( var->NetQ.data_id ) ;
	for ( int nDim = 0 ; nDim < m->nDim ; nDim ++ ) {
		//Output <<  var->EField[ nDim ] ;
		plasma.set_output( var->EField[nDim].data_id ) ;
	}
	plasma.set_output( var->E_Mag.data_id ) ;


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
	}//End streamer

	/* Production Rate */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->ProductionRate[ iSpecies ].data_id ) ;
	

	/* Mobility */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->Mobi[ iSpecies ].data_id ) ;
	

	/* Diffusion coefficients */
	for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ )	
	plasma.set_output( var->Diff[ iSpecies ].data_id ) ;
	plasma.write_output(to_string(Step)) ;

}
void CPost::PlotDecomposition( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var )
{
	// IO Output ;
	// Output.set_output_type ( config->OutputFormat ) ;
	// Output.set_filename( "MPI_ID.dat" ) ;
	// Output <<  var->MPI_ID ;
	// Output<<endl;
}
void CPost::OutputAverageFlow( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int Cycle )
{
	// IO Output ;
	// Output.set_output_type ( config->OutputFormat ) ;
	// Output.set_filename( "1AvgFlow-"+to_string(Cycle)+".dat" ) ;

	// Output <<  var->AvgPhi ;
	// Output <<  var->AvgDebyeLength ;
	// Output <<  var->AvgCFL ;
	// Output <<  var->eAvgEnergyLoss ;
	// Output <<  var->AvgEField[0] ;
	// Output <<  var->AvgEField[1] ;
	// //Output <<  var->AvgEField[2] ;
	// for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<   var->AvgT[iSpecies] ;
	// for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<  var->AvgU0[ iSpecies ] ;
	// for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<  var->AvgU1[ iSpecies ] ;
	// for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<  var->AvgU2[ iSpecies ] ;
	// //for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<  var->AvgU3[ iSpecies ] ;
	// for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<  var->AvgU4[ iSpecies ] ;
	// for ( int iSpecies=0; iSpecies < config->TotalSpeciesNum ; iSpecies++ ) Output <<  var->AvgJouleHeating[ iSpecies ] ;

	// Output<<endl;
}

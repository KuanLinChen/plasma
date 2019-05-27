#include "chemistry.h"

//int     table_new::n = 0;
double          table_new::dT = 0.0;
int             table_new::channel_n = 0;
string          separators=" :=\t\r\n";
bool            gb_Debug = false ;

inline double interpolation(double* xmin,double* xmax,double* ymin,double* ymax,double* x)
{
	return *ymin + (*x - *xmin)*(*ymax - *ymin) / (*xmax- *xmin);
}

inline double mob_interpolation(double xmin,double xmax,double ymin,double ymax,double x)
{
	return  ymin + (x - xmin)*(ymax - ymin) / (xmax- xmin);
}
/*=============================================================
|  To getting parameters from the file of species.
|  It is called in the structor of object chemistry not
|  species_identity::species_identity
==============================================================*/
void species_identity::ReadSpeciesInfomations(string  spec_file_temp)
{
		int               i,j,k;
		int       	  wordstart;
		int       	  wordend;
		ifstream  	  namefile;
		string    	  wordtemp;
		string            wordtemp2;
		vector<string>    wordunit;
		//IsLightTurnOn=0;
		namefile.open(spec_file_temp.c_str(),ios::in);
		if (!namefile) {
			cout << "Failed to open species_name.dat" << endl ; 
			exit(1) ;
		}

	  while ( getline( namefile, wordtemp) ) {

			if ( wordtemp[0] != '#' && wordtemp[0] != '\0' && wordtemp[0] != '\n' && wordtemp[0] != '\r') {
				/*
					copy from textbook page216                                                             
					clean the separators in the sentence                                                   
					Applied string::find_first_not_of(),string::length(),string::substr()                  
					string::find_first_not_of()                                                            
				*/
				wordunit.clear();                                
				wordstart=wordtemp.find_first_not_of(separators);
				wordend=0;                                       
				while (wordstart != string::npos)                
				{                                                
					wordend=wordtemp.find_first_of(separators,wordstart+1);
						if(wordend == string::npos) wordend = wordtemp.length();
						wordunit.push_back(wordtemp.substr(wordstart,wordend-wordstart));
						wordstart=wordtemp.find_first_not_of(separators,wordend+1);
				}
				//*****************************************************************************************
				if (wordunit.size() > 0 ) {

					if (wordunit.front() == "electron" || wordunit.front() == "neutral" || wordunit.front() == "ion" || wordunit.front() == "light") {

						wordtemp2 = wordunit.front() ;

					}	else {

						if (wordtemp2 == "electron") {

					    name.push_back(wordunit.front());
					    e_id=name.size()-1;

						}	else if (wordtemp2 == "neutral") {

					    name.push_back(wordunit.front());
					    neutral_id.push_back(name.size()-1);

						}	else if (wordtemp2 == "ion") {

					    name.push_back(wordunit.front());
					    ion_id.push_back(name.size()-1);

						}	else if (wordtemp2 == "light") {
						   //IsLightTurnOn=1;
					    light_name.push_back(wordunit.front());
						}

						if (wordtemp2 == "electron" || wordtemp2 == "neutral" || wordtemp2 == "ion") {

							for (i=0 ;i < wordunit.size()-1 ; ++i ) {

								if ( wordunit[i] == "diffusion_type"     ) diff_type.push_back(atoi(wordunit[i+1].c_str()));
								if ( wordunit[i] == "mobility_type"      ) mob_type.push_back(atoi(wordunit[i+1].c_str()));
								if ( wordunit[i] == "diff_constant"      ) diff_constant.push_back(atof(wordunit[i+1].c_str()));
								if ( wordunit[i] == "mob_constant"       ) mob_constant.push_back(atof(wordunit[i+1].c_str()));
								if ( wordunit[i] == "amu"                ) amu.push_back(atof(wordunit[i+1].c_str( )));
								if ( wordunit[i] == "self_diff_constant" ) self_diff_constatn.push_back(atof(wordunit[i+1].c_str( )));
								if ( wordunit[i] == "self_mob_constant"  ) self_mob_constatn.push_back(atof(wordunit[i+1].c_str( )));
								if ( wordunit[i] == "binary_diameter"    ) binary_diameter.push_back(atof(wordunit[i+1].c_str()));
								if ( wordunit[i] == "LJ_potential"       ) LJ.push_back(atof(wordunit[i+1].c_str()));
								if ( wordunit[i] == "polarizability"     ) polarizability.push_back(atof(wordunit[i+1].c_str()));
								if ( wordunit[i] == "viscosity"          ) viscosity.push_back(atof(wordunit[i+1].c_str()));
								if ( wordunit[i] == "mob_file"           ) mob_file.push_back(wordunit[i+1]);
								if ( wordunit[i] == "diff_file"          ) diff_file.push_back(wordunit[i+1]);
								if ( wordunit[i] == "thermal_file"       ) thermal_file.push_back(wordunit[i+1]);

								if ( wordunit[i] == "interpolation"      ) {
									if (wordunit[i+1] == "n" || wordunit[i+1] == "0" || wordunit[i+1] == "N"|| wordunit[i+1] == "no") {
											interpolation.push_back(0);

									}	else if(wordunit[i+1] == "y" || wordunit[i+1] == "1" || wordunit[i+1] == "Y"|| wordunit[i+1] == "yes") {
										interpolation.push_back(1);
									} else { 
										cerr << "Could not identity the option of interpolation " << endl;
									}
								}//End interpolation
							}//End information of each species 

							//******************************************
							//  Check the lastest word in species file
							//******************************************
							 
							i=wordunit.size()-1;
							if (wordunit[i] == "diffusion_type")     diff_type.push_back(1000);
							if (wordunit[i] == "mobility_type")      mob_type.push_back(1000);
							if (wordunit[i] == "diff_constant")      diff_constant.push_back(0.0);
							if (wordunit[i] == "mob_constant")       mob_constant.push_back(0.0);
							if (wordunit[i] == "amu")                amu.push_back(0.0);
							if (wordunit[i] == "self_diff_constant") self_diff_constatn.push_back(0.0);
							if (wordunit[i] == "self_mob_constant" ) self_mob_constatn.push_back(0.0);
							if (wordunit[i] == "binary_diameter")    binary_diameter.push_back(0.0);
							if (wordunit[i] == "LJ_potential")       LJ.push_back(0.0);
							if (wordunit[i] == "polarizability")     polarizability.push_back(0.0);
							if (wordunit[i] == "viscosity")          viscosity.push_back(0.0);
							if (wordunit[i] == "mob_file")           mob_file.push_back("none");
							if (wordunit[i] == "diff_file")          diff_file.push_back("none");
							if (wordunit[i] == "interpolation")      interpolation.push_back(0);
							if (wordunit[i] == "thermal_file")       thermal_file.push_back("none");

						}	else if (wordtemp2 == "light") {

							for (i=0 ;i < wordunit.size() ; ++i ) {
								if (wordunit[i] == "wave_length")    
									wavelength.push_back( atof(wordunit[i+1].c_str() ) );
							}

						}//End individual species 
					}
				}//Bypass comment
			}//if wordunit > 0 
		}//end while
		namefile.close();
		namefile.clear();

		// initial other parameter
		table_opt temp_table_for_initial;
		temp_table_for_initial.size_n = 0;
		temp_table_for_initial.dT = 0.0;
		temp_table_for_initial.max_T = 0.0;
		temp_table_for_initial.min_T = 0.0;

		for (i=0;i<name.size();++i)
		{
			thermal_table.push_back(temp_table_for_initial);
		}

		for (i=0;i<amu.size();++i)
		{
			if ( i != e_id ) {

				mass.push_back(amu[i] / 6.022142e26);

			}	else { 
				mass.push_back (9.10938e-31);
			}
		}

		for (i=0;i<amu.size();++i)
		{
			for (j=0;j<amu.size();++j)
			{
				reduce_amu.push_back(amu[i]*amu[j]/(amu[i]+amu[j]));
			}
		}

		// check number of species
		if ( (name.size()               - diff_type.size()
		     +mob_type.size()           - diff_constant.size()
		     +mob_constant.size()       - mass.size()
		     +binary_diameter.size()    - LJ.size()
		     +polarizability.size()     - viscosity.size() ) != 0 )
		{cerr << "    Size of species is not matched. Please check. " << endl;exit(1);}



/*
for (i=0;i<name.size();++i) {
cout << "  name :   " << name[i] ;
//cout << "  diffusion_type  " << diff_type[i];
//cout << "  mobility_type  " << mob_type[i];
////cout << "  diff_constant  " << diff_constant[i];
//cout << "  mob_constant  " << mob_constant[i] ;
//cout << "  mass:  "<< mass[i] ;
//cout << "  self_diff_constant  " << self_diff_constatn[i] ;
//cout << "  binary_diameter  " << binary_diameter[i];
//cout << "  LJ  " << LJ[i] ;
//cout << "  polarizability  " << polarizability[i]<< endl;
//cout << "  viscosity       " << viscosity[i] << endl;
//cout << "  the file of mob    "    << mob_file[i] << endl;
//cout << "  the file of diff    "    << diff_file[i] << endl;
cout << "  the file of thermal    "  << thermal_file[i] << endl;
////cout << "  interpolation    "    << interpolation[i] << endl;
cout << endl;
}
cout << "e id   :  " << e_id << endl;
for (i=0;i<ion_id.size();++i) {
cout << " ion   id  " << ion_id[i] << endl;
}
for (i=0;i<neutral_id.size();++i) {
cout << " neutral   id  " << neutral_id[i] << endl;
}

exit(0);

//for (k=0;k<light_name.size();++k) cerr << light_name[k] << "   " << wavelength[k] << endl;
*/

}




chemistry::chemistry( string chemistry_file, string species_temp, int cal_size )
{
	domain_size = cal_size;
	//read_input(inputname);
	species.ReadSpeciesInfomations(SpeciesFileName);
	read_reaction_table(ChannelFileName);
	loading_each_channel_file();
	newtable_establish();
	energy_loss_coll_fre_diff_mob_init();
	OptimizationSourceSink();
}

chemistry::chemistry (string  chemistry_file,string species_temp)
{
	domain_size=1;
	//read_input(inputname);
	species.ReadSpeciesInfomations(SpeciesFileName);
	read_reaction_table(ChannelFileName);
	loading_each_channel_file();
	newtable_establish();
	energy_loss_coll_fre_diff_mob_init()	;
	OptimizationSourceSink();
}
/* Generally called in our code. */
void chemistry::init(string  inputname,int cal_size)
{

	if ( gb_Debug ) {
		cout << "chemistry::init() - Debug_1" << endl ;
	}

	read_input(inputname);
	domain_size=cal_size;
	species.ReadSpeciesInfomations(SpeciesFileName);
	read_reaction_table(ChannelFileName);

    if(gb_Debug){
      cout << "chemistry::init() - Debug_2" << endl;
    }

	loading_each_channel_file();
    if(gb_Debug){
      cout << "chemistry::init() - Debug_3" << endl;
    }
	newtable_establish();
    if(gb_Debug){
      cout << "chemistry::init() - Debug_4" << endl;
    }

    if(gb_Debug){
      cout << "chemistry::init() - Debug_5" << endl;
    }

	energy_loss_coll_fre_diff_mob_init();
	OptimizationSourceSink();

    if(gb_Debug){
      cout << "chemistry::init() - Debug_6" << endl;
    }
}

void   chemistry::read_input( string react_file )
{
	string   wordtemp ;
	ifstream namefile ;
	int wordstart ;
	int wordend ;
	vector<string> wordunit ;

	namefile.open(react_file.c_str(),ios::in) ;
	/* Check file is exist*/
	if ( !namefile ) {
		cerr << "Failed to open reaction channel file" << endl ;
		exit(1) ;
	} 

	while ( getline( namefile, wordtemp ) ) {
		if ( wordtemp[0] != '#' && wordtemp[0] != '\0' && wordtemp[0] != '\n' && wordtemp[0] != '\r') {
			//*****************************************************************************************
			// copy from textbook page216                                                             *
			// clean the separators in the sentence                                                   *
			// Applied string::find_first_not_of(),string::length(),string::substr()                  *
			// string::find_first_not_of()                                                            *
			//*****************************************************************************************
			/***/wordunit.clear();                                                                   //
			/***/wordstart = wordtemp.find_first_not_of(separators);                                 //
			/***/wordend=0;                                                                          //
			/***/while (wordstart != string::npos)                                                   //
			/***/{                                                                                   //
			/***/		wordend=wordtemp.find_first_of(separators,wordstart+1);                  //
			/***/		if(wordend == string::npos) wordend = wordtemp.length();                 //
			/***/		wordunit.push_back(wordtemp.substr(wordstart,wordend-wordstart));        //
			/***/		wordstart=wordtemp.find_first_not_of(separators,wordend+1);              //
			/***/}                                                                                   //
			//*****************************************************************************************
			if (wordunit[0] == "Electron_table_size" ) table_e_size  = atoi( wordunit[1].c_str() ) ;
			if (wordunit[0] == "Table_Te_min" )        table_Te_min  = atof( wordunit[1].c_str() ) ;
			if (wordunit[0] == "Table_Te_max" )        table_Te_max  = atof( wordunit[1].c_str() ) ;
			if (wordunit[0] == "Table_gas_min" )       table_gas_min = atof( wordunit[1].c_str() ) ;
			if (wordunit[0] == "Table_gas_max" )       table_gas_max = atof( wordunit[1].c_str() ) ;
			if (wordunit[0] == "T_Gas" )               T_Gas         = atof( wordunit[1].c_str() ) ;
			if (wordunit[0] == "SPECIES" )             SpeciesFileName = wordunit[1] ;
			if (wordunit[0] == "CHEMISTRY" )           ChannelFileName = wordunit[1] ;
			if (wordunit[0] == "Ion_table_size" )      table_ion_size 		= atoi(wordunit[1].c_str()) ;
			if (wordunit[0] == "Neutral_table_size" )  table_neutral_size = atoi(wordunit[1].c_str()) ;
		}
	}//End read file
		namefile.close();
		namefile.clear();

}



/*============================================================
| It is a function to handle reaction file to get parameter. |
| Be called by structor Object chemistry.                    |
| reaction channel file                                      |
=============================================================*/
void   chemistry::read_reaction_table(string react_file)
{
	single_reaction   *temp;
	int       	  i,j,k;
	int               temp_int;
	int       	  wordstart;
	int       	  wordend;
	int       	  s_s_flag;
	ifstream  	  reac_file;
	string            wordtemp;
	vector<int>       pow_index;
	vector<string>    wordunit;
	vector<string>::iterator  iter_i;
	vector<int>::iterator     iter_j;

	reac_file.open(react_file.c_str(),ios::in);
	if (!reac_file) {
		cerr << "Failed to open reaction_channel.list" << endl ; 
		exit(1) ;
	} else {
		while ( getline( reac_file, wordtemp) ) {
			wordunit.clear();
			if (wordtemp.length() > 3 && wordtemp[0] != '#' && wordtemp[0] != '\0' && wordtemp[0] != '\n' && wordtemp[0] != '\r') {

				wordstart=wordtemp.find_first_not_of(separators);
		 		wordend=0;

				while (wordstart != string::npos)
				{
					wordend=wordtemp.find_first_of(separators,wordstart+1);
					if(wordend == string::npos) wordend = wordtemp.length();
					wordunit.push_back(wordtemp.substr(wordstart,wordend-wordstart));
					wordstart=wordtemp.find_first_not_of(separators,wordend+1);
				}

				s_s_flag=0;
				for ( i=0 ; i < wordunit.size() ; ++i ) {
					if ( wordunit[i] == "->" ) s_s_flag = 19028 ;
				}
				if (s_s_flag == 0) {
					cerr << "The reaction of number " << i <<" has no sign '->'."<<endl ; 
					exit(1);
				}

				s_s_flag=0;
				i=0;

				/* Chech the name in species file and channel file are match ? */
				while ( s_s_flag == 0 ) {

					if ( isdigit( wordunit[i][0] ) != 0 ) {
						temp_int=0;
						for (k=0;k<species.name.size();++k)	{
							if (wordunit[i+1] == species.name[k]) temp_int=8736;
						}
						for (k=0;k<species.light_name.size();++k) {
							if (wordunit[i+1] == species.light_name[k])  temp_int = 89878;
						}

						if (temp_int == 0) {
							cout << "The species of species file and reaction channel file are not matched."<<endl;
							exit(1);
						}
					}

					if ( wordunit[i] == "->" || wordunit[i] == "ratefile" ) s_s_flag=9173;
					i=i+1;

				}//end while


				while (wordunit[i] != "ratefile" && wordunit[i] != "AAA") {

					if ( isdigit(wordunit[i][0])!=0 ) {

						temp_int=0;
						for (k=0;k<species.name.size();++k) {
							if (wordunit[i+1] ==species.name[k]) temp_int = 1000 ;
						}
							//cout << wordunit[i+1] << endl;
						for ( k=0 ; k < species.light_name.size() ; ++k ) {
						    if (wordunit[i+1] == species.light_name[k])  temp_int = 89878;
						}
						//cout <<temp_int << endl;
						if ( temp_int == 0 ) {
							cout << "The reaction channel file are not matched. In read_reaction_table"<<endl;
							cout << "The species "<< wordunit[i+1] <<" do not be found in species_name.txt"<< endl;
							exit(1);
						}
					}//end if
					i=i+1;
				}//end while
			}
		}
	}
	reac_file.close();
	reac_file.clear();


	for ( i=0 ; i < species.name.size() ; ++i ) pow_index.push_back(0);

	reac_file.open(react_file.c_str(),ios::in);

	if ( !reac_file ) {
		cerr << "Failed to open reaction_channel.list" << endl ; 
		exit(1) ;

	}	else {

		while ( getline( reac_file, wordtemp ) ) {

			temp = new single_reaction ;

			/************************************************************
			*      		Initial temp vector
			*************************************************************/
			for (i=0;i<species.name.size();++i) temp->source.push_back(0);
			for (i=0;i<species.name.size();++i) temp->sink.push_back(0);
			for (i=0;i<4;++i) temp->formula_coefficient.push_back(0.0);
			for (i=0;i<domain_size;++i) temp->Channel_SourceSink.push_back(0.0);
			temp->threshold=0.0;
			temp->rate_constant=0.0;
			temp->is_data_output=0;
			temp->is_interpolation=0;
			temp->emit_id = -1;
			/*************************************************************/
			wordunit.clear();
			if (wordtemp.length() > 3 && wordtemp[0] != '#' && wordtemp[0] != '\0' && wordtemp[0] != '\n') {
			        wordstart=wordtemp.find_first_not_of(separators);
			        wordend=0;
				while (wordstart != string::npos){
					wordend=wordtemp.find_first_of(separators,wordstart+1);
					if(wordend == string::npos) wordend = wordtemp.length();
					wordunit.push_back(wordtemp.substr(wordstart,wordend-wordstart));
					wordstart=wordtemp.find_first_not_of(separators,wordend+1);
				}
				s_s_flag=0; j=0;
				for (i=0;i<species.name.size();++i) pow_index[i]=0;
				for (iter_i=wordunit.begin();iter_i<wordunit.end();++iter_i) {

					if ((*iter_i).compare("threshold")==0)  temp->threshold=atof((*(iter_i+1)).c_str());
					if ((*iter_i).compare("->")==0) s_s_flag=1;
					if (s_s_flag == 0) {
						for (i=0;i<species.name.size();++i) {
							if ( (*iter_i).compare(species.name[i]) == 0) {
								temp->sink[i] = temp->sink[i] + atoi((*(iter_i-1)).c_str());
								j=j+1;
								pow_index[i]=j;
							}
						}
						for (i=0;i<species.light_name.size();++i)
						{
							if ((*iter_i) == species.light_name[i])
							{
								cerr << "light absorption is not include." << endl;
							}
						}
					}
					else{
						for (i=0;i<species.name.size();++i) {
							if ( (*iter_i).compare(species.name[i]) == 0)
							temp->source[i] = temp->source[i] + atoi((*(iter_i-1)).c_str());
						}
						for (i=0;i<species.light_name.size();++i)
						{
							if ((*iter_i) == species.light_name[i])
							{
								temp->emit_id = i;
//cerr << " light :  " << species.light_name[i] << endl;
							}

						}
					}
					if ((*iter_i).compare("ratefile")==0) {temp->filename = *(iter_i+1);}
					if ((*iter_i).compare("type")==0 )
					{
						if (*(iter_i+1) == "formula" || *(iter_i+1) == "f" || *(iter_i+1) == "1")  {temp->rate_type = 1;}
						else if (*(iter_i+1) == "table" || *(iter_i+1) == "t" || *(iter_i+1) == "2") {temp->rate_type = 2;}
						else if (*(iter_i+1) == "constant" || *(iter_i+1) == "c" || *(iter_i+1) == "3") {temp->rate_type = 3;}
						else if (*(iter_i+1) == "gas_table" || *(iter_i+1) == "gt" || *(iter_i+1) == "4") {temp->rate_type = 4;}
						else { cerr << "The reaction type you insert could not identify. Please check." << endl;exit(100);}
					}
					if ((*iter_i) == "rateconstant") temp->rate_constant=atof((*(iter_i+1)).c_str());
					if ((*iter_i) == "data_output")
					{
						if (*(iter_i+1) == "yes" || *(iter_i+1) == "y" || *(iter_i+1) == "1")  {temp->is_data_output = 1;}
						else if (*(iter_i+1) == "no" || *(iter_i+1) == "n" || *(iter_i+1) == "0") {temp->is_data_output = 0;}
						else { cerr << "The data_output you insert could not identify. Please check. It only yes or no." << endl;exit(100);}
					}
					if ((*iter_i) == "outname") {temp->output_name = *(iter_i+1);}
					if ((*iter_i) == "interpolation")
					{
						if (*(iter_i+1) == "yes" || *(iter_i+1) == "y" || *(iter_i+1) == "1")  {temp->is_interpolation = 1;}
						else if (*(iter_i+1) == "no" || *(iter_i+1) == "n" || *(iter_i+1) == "0") {temp->is_interpolation = 0;}
						else { cerr << "The interpolation you insert could not identify. Please check. It only yes or no." << endl;exit(100);}
					}
				}
				channel.push_back(*temp);
			}
			delete temp;
		    }
		}
		reac_file.close();
		reac_file.clear();
		wordunit.clear();
}


/**
 * @brief Read ecah channel file.
 * @date 05/27/2019
 */
void chemistry::loading_each_channel_file()
{
	ifstream        rate_file;
	string          inword;
	vector<string>  wordunit;
	int             wordstart;
	int             wordend;
	double          temp_double;

	for ( int i = 0 ; i < channel.size() ; ++i ) {

		/* rate_type = foumula */
		if ( channel[ i ].rate_type == FORMULA ) {

			wordunit.clear();
			rate_file.open(channel[i].filename.c_str(),ios::in);

			/* Check file is exist*/
			if (!rate_file) {
				cerr << "Failed to open "<< channel[i].filename<< endl ;
				exit(1) ;
			}
			/* Star read the file data. */
			while(getline(rate_file,inword)) {

				if (inword[0] != '#' && inword[0] != '\n' && inword[0] != '\0' && inword[0] != '\r') {

					wordstart = inword.find_first_not_of(separators);

					wordend = 0 ;

					while (wordstart != string::npos) {

						wordend=inword.find_first_of(separators,wordstart+1);
						if(wordend == string::npos) wordend = inword.length();
						wordunit.push_back(inword.substr(wordstart,wordend-wordstart));
						wordstart=inword.find_first_not_of(separators,wordend+1);
					}

				}//bypass the comment

			}//End read file

			/* rate = A1*(A2/Tg)^A3 exp(A4/T) */
			for ( int j = 0; j < wordunit.size() ; ++j ) {
				if(wordunit[j] == "A1") channel[i].formula_coefficient[0]=atof(wordunit[j+1].c_str());
				if(wordunit[j] == "A2") channel[i].formula_coefficient[1]=atof(wordunit[j+1].c_str());
				if(wordunit[j] == "A3") channel[i].formula_coefficient[2]=atof(wordunit[j+1].c_str());
				if(wordunit[j] == "A4") channel[i].formula_coefficient[3]=atof(wordunit[j+1].c_str());
			}

			rate_file.close();
			rate_file.clear();

		/* rate_type = table */
		} else if (channel[i].rate_type == TABLE ) {

			rate_file.open(channel[i].filename.c_str(),ios::in);
			if (!rate_file) {cerr << "Failed to open " << channel[i].filename << endl; exit(1);}
			else
			{
				while (!rate_file.eof())
				{
					rate_file >> inword ;

					if (inword == "EndDescription")
					{
						while (rate_file >> temp_double)
						{
							channel[i].origin_eV.push_back(temp_double);
							rate_file >> temp_double;
							channel[i].origin_rate.push_back(temp_double);
						}
					}
					else if (isdigit(inword[0]) != 0 || isdigit(inword[1]) != 0)
					{
						cerr << "Waring!!!  please insert the word \"EndDescription\" in front of reaction file" << endl;
						cerr << "The name of file is " << channel[i].filename << endl; exit(1);
					}
				}
			}
			rate_file.close();
			rate_file.clear();

		/* rate_type = constant */
		} else if (channel[i].rate_type == CONSTANT ) {

			//for constant rate type, we don't need to read the file.

		/* rate_type = gas table */
		}	else if ( channel[i].rate_type == GAS_TABLE )	{
			rate_file.open(channel[i].filename.c_str(),ios::in);

			if (!rate_file) {
				cerr << "Failed to open " << channel[i].filename << endl; 
				exit(1);
			} 

			while ( !rate_file.eof() ) {

				rate_file >> inword ;

				/*read data after "EndDescription" */
				if ( inword == "EndDescription") {

					while ( rate_file >> temp_double ) {

						channel[i].origin_eV.push_back(temp_double);
						rate_file >> temp_double;
						channel[i].origin_rate.push_back(temp_double);
					}

				}	else if ( isdigit(inword[0]) != 0 || isdigit(inword[1]) != 0 ) {

					cerr << "Waring!!!  please insert the word \"EndDescription\" in front of reaction file" << endl;
					cerr << "The name of file is " << channel[i].filename << endl; exit(1);
				}
			}//End while
			rate_file.close();
			rate_file.clear();
			/*
			for (j=0;j<channel[i].origin_eV.size();++j)
			{
				cout <<"dkldkldjldjdkl;jddkljldjd    " << "  T=" << channel[i].origin_eV[j] << "  rate" << channel[i].origin_rate[j] << endl;
			}
			*/
		}
		else{
			cerr << "The error is happened in reaction. Please check reaction list data." << endl;
			exit(1);
		}
	}

}



void chemistry::newtable_establish()
{
	int               i,j,k;
	int               l;
	int               counter;
	ofstream          out_file;
	vector<int>       vectemp;
	ifstream  	  difffile;
	ifstream          mobfile;
	ifstream          ThermalFile;
	string            inword;
	//int table_new::n = 0;
        //double table_new::dT = 0.0;
	//table_e_size=50000;
	//table_Te_min = 0.08;
        //table_Te_max = 17.0;
	//T_Gas = 400.0;

    if(gb_Debug){
      cout << "chemistry::newtable_establish() - Debug_1" << endl;
    }

	channel_new_table = new table_new [table_e_size];
	channel_to_table = new int [channel.size()];
	counter = 0;
	for (l=0;l<channel.size();++l)
	{

		if (channel[l].rate_type == 2)
		{
			vectemp.push_back(l);
			channel_to_table[l]=vectemp.size()-1;
		}
		else{channel_to_table[l] = channel.size();}
	}


	table_new::channel_n = vectemp.size();
	which_is_table = new int [vectemp.size()];
	for (l=0;l<vectemp.size();++l)
	{
		which_is_table[l] = vectemp[l];
	}
	for (k=0;k<table_e_size;++k)
	{
		channel_new_table[k].R = new double [table_new::channel_n];
		channel_new_table[k].dR =new double [table_new::channel_n];
	}
	table_new::dT = ( table_Te_max - table_Te_min ) /(table_e_size -1);
	channel_new_table[0].T = table_Te_min;
	for (i=1;i<table_e_size;++i)
	{
		channel_new_table[i].T = channel_new_table[i-1].T  + channel_new_table[i-1].dT ;
	}

    if(gb_Debug){
      cout << "chemistry::newtable_establish() - Debug_2" << endl;
    }

	/************************************************************
	    build new table of Te and rate
	    Form:
	    Te   channel[1]   channel[2]  ....      e_diff   e_mob
	    0.1    value        value      value    value    value
	    0.2    value        value      value    value    value
	    0.3    value        value      value    value    value

	************************************************************/
	for (l=0;l<table_new::channel_n;++l)
	{
		k = which_is_table[l];
		if (channel[k].rate_type == 1 )
		{
			channel_new_table[0].R[l] = rate_formula(& channel_new_table[0].T ,k);
			for (i=1;i<table_e_size;++i)
			{
				channel_new_table[i].R[l] = rate_formula( &channel_new_table[i].T ,k);
			}
			// partial Rate / partial Te
			for (i=0 ; i< table_e_size ; ++i)
			{
				if ( i == 0 ) { channel_new_table[i].dR[l] = (channel_new_table[i+1].R[l] - channel_new_table[i].R[l]) / table_new::dT ;  }
				else { channel_new_table[i].dR[l] = (channel_new_table[i].R[l] - channel_new_table[i-1].R[l]) / table_new::dT ;  }
			}
		}
		else if (channel[k].rate_type == 2)
		{
			j=0;
			for (i=0 ; i< table_e_size ; ++i)
			{
				if (channel_new_table[i].T <= channel[k].origin_eV[0])
				{
					channel_new_table[i].R[l]= interpolation (& channel[k].origin_eV[0],& channel[k].origin_eV[1],& channel[k].origin_rate[0],& channel[k].origin_rate[1],& channel_new_table[i].T);
					if  (channel_new_table[i].R[l] < 0.0) channel_new_table[i].R[l] = channel[k].origin_rate[0];
				}
				else if (channel_new_table[i].T > channel[k].origin_eV.front() && channel_new_table[i].T <= channel[k].origin_eV.back())
				{
					if (j < (channel[k].origin_eV.size()-2) )
					{
					    while ( channel_new_table[i].T <= channel[k].origin_eV[j] || channel_new_table[i].T > channel[k].origin_eV[j+1] )
					    {
						    if ( channel_new_table[i].T >= channel[k].origin_eV[j+1] )  j=j+1;
						    else if ( channel_new_table[i].T <= channel[k].origin_eV[j] )    j=j-1;
					    }
					}
					channel_new_table[i].R[l]= interpolation (&channel[k].origin_eV[j] , &channel[k].origin_eV[j+1] ,& channel[k].origin_rate[j] , & channel[k].origin_rate[j+1] ,& channel_new_table[i].T);
				}
				else
				{
					channel_new_table[i].R[l] = interpolation (& channel[k].origin_eV[channel[k].origin_eV.size()-2] ,& channel[k].origin_eV[channel[k].origin_eV.size()-1]
										  ,& channel[k].origin_rate[channel[k].origin_rate.size()-2], & channel[k].origin_rate[channel[k].origin_rate.size()-1]
										  ,& channel_new_table[i].T);
				}

			}
			// partial Rate / partial Te
			for (i=0 ; i< table_e_size ; ++i)
			{
				if ( i == 0 ) { channel_new_table[i].dR[l] = (channel_new_table[i+1].R[l] - channel_new_table[i].R[l]) / table_new::dT ;  }
				else { channel_new_table[i].dR[l] = (channel_new_table[i].R[l] - channel_new_table[i-1].R[l]) / table_new::dT ;  }
			}
		}
		else if (channel[k].rate_type == 3)
		{
			cerr << "   The Reaction " << channel[k].filename << "is constant" << endl;
			cerr << "  Chemistry module will get out. " << endl;
			exit(1);
		}
		else if (channel[k].rate_type == 4)
		{
			cerr << "   The Reaction " << channel[k].filename << "is gas table" << endl;
			cerr << "  Chemistry module will get out. " << endl;
			exit(1);
		}
		else
		{
			cerr << "The error is happened in reaction. Please check reaction list data."<< endl;
			exit(1);
		}
	}

	    /*
	    k = which_is_table[l];
	    for (i=0 ; i< table_e_size ; ++i)
	    {
		    cout << channel_new_table[i].T;

		    for (l=0;l<table_new::channel_n;++l)
		    {
			    cout  << "  " << channel_new_table[i].R[l]  << channel_new_table[i].dR[l];
			    cout  <<"  "<< channel_new_table[i].dR[l] << endl;
		    }
		    cout << endl;
	    }
	    */

    if(gb_Debug){
      cout << "chemistry::newtable_establish() - Debug_3" << endl;
    }

		/**********************************************************
			Rebulit diffusion in single species
			diff_type == 6
			diff_type == 7
		***********************************************************/
		double  tmp_Te,tmp_value;
		vector <double> Tee;
		vector <double> value;
		k=species.e_id;
		if (species.diff_type[k] == 5 || species.diff_type[k] == 6)
		{
			difffile.open(species.diff_file[k].c_str(),ios::in);
			if (!difffile) {cerr << "Failed to open file of electron diffusivity" << endl; exit(1);}
			else
			{
				while (!difffile.eof())
				{
					difffile >> inword ;
					if (inword == "EndDescription")
					{
						while(difffile >> tmp_Te >> tmp_value)
						{
							Tee.push_back(tmp_Te);
							value.push_back(tmp_value);
						}
					}
					else if (isdigit(inword[0]) != 0 || isdigit(inword[1]) != 0)
					{
						cerr << "Waring!!!  please insert the word \"EndDescription\" in front of reaction file" << endl;
						cerr << "The name of file is " << channel[i].filename << endl; exit(1);
					}
				}
			}
			difffile.close();
		        difffile.clear();
			j=0;
			for (i=0 ; i< table_e_size ; ++i)
			{
				if (channel_new_table[i].T <= Tee[0] )
				{
					channel_new_table[i].diff= interpolation (& Tee[0],& Tee[1],& value[0],& value[1],& channel_new_table[i].T);
					if  (channel_new_table[i].diff < 0.0) channel_new_table[i].diff = value[0];
				}
				else if (channel_new_table[i].T > Tee[0] && channel_new_table[i].T <= Tee.back())
				{
					if (j < (Tee.size()-2) )
					{
					    while ( channel_new_table[i].T <= Tee[j] || channel_new_table[i].T > Tee[j+1] )
					    {
						    if ( channel_new_table[i].T >= Tee[j+1] )  j=j+1;
						    else if ( channel_new_table[i].T <= Tee[j] )    j=j-1;
					    }
					}
					channel_new_table[i].diff = interpolation (& Tee[j] , &Tee[j+1] ,& value[j] , & value[j+1] ,& channel_new_table[i].T);
				}
				else
				{
					channel_new_table[i].diff = interpolation (& Tee[Tee.size()-2] ,& Tee[Tee.size()-1]
										  ,& value[value.size()-2], & value[value.size()-1]
										  ,& channel_new_table[i].T);
				}

			}
			// partial diffusion / partial Te
			for (i=0 ; i< table_e_size ; ++i)
			{
				if ( i == 0 ) { channel_new_table[i].ddiff = (channel_new_table[i+1].diff - channel_new_table[i].diff) / table_new::dT ;  }
				else { channel_new_table[i].ddiff = (channel_new_table[i].diff - channel_new_table[i-1].diff) / table_new::dT ;  }
			}
			Tee.clear();
			value.clear();
		}

    if(gb_Debug){
      cout << "chemistry::newtable_establish() - Debug_4" << endl;
    }

		/**********************************************************
			Rebulit Mobility in single species
			mob_type == 6
			mob_type == 7
		***********************************************************/
		k=species.e_id;
		if (species.mob_type[k] == 5 || species.mob_type[k] == 6)
		{
			mobfile.open(species.mob_file[k].c_str(),ios::in);
			if (!mobfile) {cerr << "Failed to open file of electron mobility" << endl; exit(1);}
			else
			{
				while (!mobfile.eof())
				{
					mobfile >> inword ;
					if (inword == "EndDescription")
					{
						while(mobfile >> tmp_Te >> tmp_value)
						{
							Tee.push_back(tmp_Te);
							value.push_back(tmp_value);
						}
					}
					else if (isdigit(inword[0]) != 0 || isdigit(inword[1]) != 0)
					{
						cerr << "Waring!!!  please insert the word \"EndDescription\" in front of reaction file" << endl;
						cerr << "The name of file is " << channel[i].filename << endl; exit(1);
					}
				}
			}
			mobfile.close();
		        mobfile.clear();
			j=0;
			for (i=0 ; i< table_e_size ; ++i)
			{
				if (channel_new_table[i].T <= Tee[0] )
				{
					channel_new_table[i].mob= interpolation (& Tee[0],& Tee[1],& value[0],& value[1],& channel_new_table[i].T);
					if  (channel_new_table[i].mob < 0.0) channel_new_table[i].mob = value[0];
				}
				else if (channel_new_table[i].T > Tee[0] && channel_new_table[i].T <= Tee.back())
				{
					if (j < (Tee.size()-2) )
					{
					    while ( channel_new_table[i].T <= Tee[j] || channel_new_table[i].T > Tee[j+1] )
					    {
						    if ( channel_new_table[i].T >= Tee[j+1] )  j=j+1;
						    else if ( channel_new_table[i].T <= Tee[j] )    j=j-1;
					    }
					}
					channel_new_table[i].mob = interpolation (& Tee[j] , &Tee[j+1] ,& value[j] , & value[j+1] ,& channel_new_table[i].T);
				}
				else
				{
					channel_new_table[i].mob = interpolation (& Tee[Tee.size()-2] ,& Tee[Tee.size()-1]
										  ,& value[value.size()-2], & value[value.size()-1]
										  ,& channel_new_table[i].T);
				}

			}
			// partial diffusion / partial Te
			for (i=0 ; i< table_e_size ; ++i)
			{
				if ( i == 0 ) { channel_new_table[i].dmob = (channel_new_table[i+1].mob - channel_new_table[i].mob) / table_new::dT ;  }
				else { channel_new_table[i].dmob = (channel_new_table[i].mob - channel_new_table[i-1].mob) / table_new::dT ;  }
			}
			Tee.clear();
			value.clear();
		}

    if(gb_Debug){
      cout << "chemistry::newtable_establish() - Debug_5" << endl;
    }
		/**********************************************************
			Rebulit Other in single species
			thermal table
		***********************************************************/

		for (k=0;k<species.name.size();++k){

			if(gb_Debug){cout << "chemistry::newtable_establish() - Debug_5.0" << endl;}

			if (species.thermal_file[k] != "none"){

				if(gb_Debug){cout << "chemistry::newtable_establish() - Debug_5.1" << endl;}

				ThermalFile.open(species.thermal_file[k].c_str(),ios::in);

				if (!ThermalFile){

					cerr << "Failed to open species_name.dat" << endl; exit(1);

				}else{

				if(gb_Debug){cout << "chemistry::newtable_establish() - Debug_5.2" << endl;}
					while ( !ThermalFile.eof() ){

						ThermalFile >> inword ;
						if (inword == "EndDescription"){

							while(ThermalFile >> tmp_Te >> tmp_value)
							{
								Tee.push_back(tmp_Te);
								value.push_back(tmp_value);
							}
						}else if (isdigit(inword[0]) != 0 || isdigit(inword[1]) != 0){

							cerr << "Waring!!!  please insert the word \"EndDescription\" in front of reaction file" << endl;
							cerr << "The name of file is " << channel[i].filename << endl; exit(1);

						}
					}
				}
				species.thermal_table[k].min_T 	= Tee.front() ;
				species.thermal_table[k].max_T	= Tee.back();
				species.thermal_table[k].size_n   	= table_neutral_size;
				species.thermal_table[k].dT        	= (species.thermal_table[k].max_T - species.thermal_table[k].min_T)/(species.thermal_table[k].size_n-1);
				species.thermal_table[k].T         	= new double [species.thermal_table[k].size_n];
				species.thermal_table[k].value     	= new double [species.thermal_table[k].size_n];
				species.thermal_table[k].T[0] 		= Tee.front();
				species.thermal_table[k].value[0] 	= value.front();

				for (i=1;i<species.thermal_table[k].size_n-1;++i){

					tmp_Te=species.thermal_table[k].min_T+ i*species.thermal_table[k].dT;

					for (j=0;j<Tee.size()-1;++j){

						if (tmp_Te >= Tee[j] && tmp_Te<Tee[j+1]) counter=j;
					}
					species.thermal_table[k].value[i] = mob_interpolation(Tee[counter],Tee[counter+1],value[counter],value[counter+1],tmp_Te);
					species.thermal_table[k].T[i] = tmp_Te;
				}
				species.thermal_table[k].T[species.thermal_table[k].size_n-1] = Tee.back();
				species.thermal_table[k].value[species.thermal_table[k].size_n-1] = value.back();
				/*
				for (i=0;i<species.thermal_table[k].size_n;++i){
					cout << species.thermal_table[k].T[i] << "  " <<species.thermal_table[k].value [i] << endl;
				}
				*/
			}
			Tee.clear();
			value.clear();
			ThermalFile.clear();
			ThermalFile.close();


		}

    if(gb_Debug){
      cout << "chemistry::newtable_establish() - Debug_6" << endl;
    }

}


/**
 * @brief      Calculate rate by formula, rate=A1 x (A2/Tg)^A3 exp(A4/Tg) 
 *
 * @param[in]  temperature  The temperature.
 * @param[in]  k            The channel.
 *
 * @return     rate constant of channel k.
 */
double  chemistry::rate_formula(double *temperature, int k ) 
{
	double rate=0.0 ;
	if ( channel[k].formula_coefficient[3] == 0.0 && channel[k].formula_coefficient[2] != 0.0) {

		rate=channel[k].formula_coefficient[0]*pow(channel[k].formula_coefficient[1]/(*temperature),channel[k].formula_coefficient[2]);

	}	else if (channel[k].formula_coefficient[2] == 0.0 && channel[k].formula_coefficient[3] != 0.0) {

		rate=channel[k].formula_coefficient[0] * exp(channel[k].formula_coefficient[3]/(*temperature));

	}	else if (channel[k].formula_coefficient[2] == 0.0 && channel[k].formula_coefficient[3] == 0.0) {

		rate=channel[k].formula_coefficient[0];

	} else {

		rate = channel[k].formula_coefficient[0]*pow(channel[k].formula_coefficient[1]/(*temperature),channel[k].formula_coefficient[2])
				 * exp(channel[k].formula_coefficient[3]/(*temperature));

	}

	return rate;
}

/*========================================================
|  It is the function to search rate in new table        |
|  2009-04-20                                            |
|  In: Electric temperature                              |
|  Out: A series rate of domain (double *)               |
=========================================================*/
inline double*     chemistry::search_table(const double *Te)
{
	//int index;
	if (*Te <= channel_new_table[0].T)
	{
		return & channel_new_table[0].R[0];
	}
	else if (*Te >= channel_new_table[table_e_size-1].T)
	{
		return & channel_new_table[table_e_size-1].R[0];
	}
	else
	{
		//index=static_cast<int> ( ( *Te- channel_new_table[0].T ) / table_new::dT );
		return & channel_new_table[Te_index].R[0];
	}
	cerr << " Error at serach table exit.  " << endl;
	exit(1);
}

inline double*     chemistry::search_perturb(const double *Te)
{
	//int index;
	if (*Te <= channel_new_table[0].T)
	{
		return & channel_new_table[0].dR[0];
	}
	else if (*Te >= channel_new_table[table_e_size-1].T)
	{
		return & channel_new_table[table_e_size-1].dR[0];
	}
	else
	{
		//index=static_cast<int> ( ( *Te- channel_new_table[0].T ) / table_new::dT );
		return & channel_new_table[Te_index].dR[0];
	}
	cerr << " Error at serach table exit.  " << endl;
	exit(1);
}

inline double*     chemistry::search_diffusion(const double *Te)
{
	//int index;
	if (*Te <= channel_new_table[0].T)
	{
		return & channel_new_table[0].diff;
	}
	else if (*Te >= channel_new_table[table_e_size-1].T)
	{
		return & channel_new_table[table_e_size-1].diff;
	}
	else
	{
		//index=static_cast<int> ( ( *Te- channel_new_table[0].T ) / table_new::dT );
		return & channel_new_table[Te_index].diff;
	}
	cerr << " Error at serach table exit.  " << endl;
	exit(1);
}

inline double*     chemistry::search_mobility(const double *Te)
{
	//int index;
	if (*Te <= channel_new_table[0].T)
	{
		return & channel_new_table[0].mob;
	}
	else if (*Te >= channel_new_table[table_e_size-1].T)
	{
		return & channel_new_table[table_e_size-1].mob;
	}
	else
	{
		//index=static_cast<int> ( ( *Te- channel_new_table[0].T ) / table_new::dT );
		return & channel_new_table[Te_index].mob;
	}
	cerr << " Error at serach table exit.  " << endl;
	exit(1);
}

inline double*     chemistry::search_ddiff(const double *Te)
{
	//int index;
	if (*Te <= channel_new_table[0].T)
	{
		return & channel_new_table[0].ddiff;
	}
	else if (*Te >= channel_new_table[table_e_size-1].T)
	{
		return & channel_new_table[table_e_size-1].ddiff;
	}
	else
	{
		//index=static_cast<int> ( ( *Te- channel_new_table[0].T ) / table_new::dT );
		return & channel_new_table[Te_index].ddiff;
	}
	cerr << " Error at serach table exit.  " << endl;
	exit(1);
}

inline double*     chemistry::search_dmob(const double *Te)
{
	//int index;
	if (*Te <= channel_new_table[0].T)
	{
		return & channel_new_table[0].dmob;
	}
	else if (*Te >= channel_new_table[table_e_size-1].T)
	{
		return & channel_new_table[table_e_size-1].dmob;
	}
	else
	{
		//index=static_cast<int> ( ( *Te- channel_new_table[0].T ) / table_new::dT );
		return & channel_new_table[Te_index].dmob;
	}
	cerr << " Error at serach table exit.  " << endl;
	exit(1);
}



/*********************************************************************************************************************
**********************************************************************************************************************

                    Rate constant, source-sink constant, collision frequence and e energy loss
2009-04-01:  No consider gas temperature
**********************************************************************************************************************
**********************************************************************************************************************/


void   chemistry::CalSourceSinkRate_for_Te(const double *T_species )
{
	//Mark by KL: int     i,j,
	int 	k,l;
	int     domain_int;
	int     buf_start;
	int     buf_start_channel;
	double  average_Te;
	double  *new_rate_buf;
	double  *new_rate_perturb_buf;
	double  *new_diff_buf;
	double  *new_mob_buf;
	double  *new_ddiff_buf;
	double  *new_dmob_buf;



	for (domain_int=0;domain_int<domain_size;++domain_int)
	{
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;

		average_Te=( *(T_species + species.e_id + buf_start) );
		Te_index = static_cast<int> ( ( average_Te - channel_new_table[0].T ) / table_new::dT );
		new_rate_buf = search_table(&average_Te);
		new_rate_perturb_buf = search_perturb(&average_Te);
		new_diff_buf = search_diffusion(&average_Te);
		new_mob_buf = search_mobility(&average_Te);
		new_ddiff_buf = search_ddiff(&average_Te);
	    new_dmob_buf = search_dmob(&average_Te);

		for (l=0;l<table_new::channel_n;++l)
		{
			k = which_is_table[l];
			chem_rate_buf[k+buf_start_channel] = *(new_rate_buf+l);
			chem_rate_perturb_buf[k+buf_start_channel] = *(new_rate_perturb_buf+l);
		}

		//To speed up
		if (species.mob_type[species.e_id] == 5)
		{
			mob_buf[domain_int] = *new_mob_buf / total_gas_density[domain_int];
			dmob_buf[domain_int] = *new_dmob_buf / total_gas_density[domain_int];
		}
		else if (species.mob_type[species.e_id] == 6)
		{
			mob_buf[domain_int] = *new_mob_buf ;
			dmob_buf[domain_int] = *new_dmob_buf ;
		}
		else
		{
			mob_buf[domain_int] = 0.0 ;
			dmob_buf[domain_int] = 0.0 ;
		}

		if (species.diff_type[species.e_id] == 5)
		{
			diff_buf[domain_int] = *new_diff_buf / total_gas_density[domain_int];
			ddiff_buf[domain_int] = *new_ddiff_buf / total_gas_density[domain_int];
		}
		else if (species.diff_type[species.e_id] == 6)
		{
			diff_buf[domain_int] = *new_diff_buf;
			ddiff_buf[domain_int] = *new_ddiff_buf;
		}
		else
		{
			diff_buf[domain_int] = 0.0;
			ddiff_buf[domain_int] = 0.0;
		}




	}



	/*
	for (domain_int=0;domain_int<domain_size;++domain_int){
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;
		for (k=0;k<channel_size;++k)
		{
		cout   << "domain: "<< domain_int <<" reaction:  " << k<< "  Te:  "<< *(T_species+species.e_id+buf_start) <<"  rate:   " << chem_rate_buf[k+buf_start_channel] << endl;
		}
		cout << endl;
	}


	for (domain_int=0;domain_int<domain_size;++domain_int){
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;
		for (k=0;k<channel_size;++k)
		{
		cout   << "domain: "<< domain_int <<" reaction:  " << k<< "  Te:  "<< *(T_species+species.e_id+buf_start) <<"  rate:   " << chem_rate_perturb_buf[k+buf_start_channel] << endl;
		}
		cout << endl;
	}

	for (domain_int=0;domain_int<domain_size;++domain_int){
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;
		for (k=0;k<channel_size;++k)
		{
		cout   << "domain: "<< domain_int <<" reaction:  " << k<< "  Te:  "<< *(T_species+species.e_id+buf_start) <<"  rate:   " << endl;
		}
		cout << endl;
	}
	exit(1) ;
	*/

}
/**
 * @brief      Calculate source/sink for gas. (relate to the gas temperature)
 *
 * @param[in]  T_species  The temperature of species.
 * @date 04/01/2009
 */
void   chemistry::CalSourceSinkRate_for_gas(const double *T_species )
{
	int        i, l;
	int        domain_int,coefficient;
	int        buf_start, buf_start_channel;
	double     buf_temperature, counter_duoble;

	for (domain_int=0;domain_int<domain_size;++domain_int)
	{
		buf_start         = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;

		for (l=0;l<channel_size;++l)
		{
			if (channel[l].rate_type == 4)
			{
				buf_temperature = 0.0;
				counter_duoble =0.0;
				for (i=0;i< OPT_channel_sourcesink[l].sink_n ; ++i)
				{
					coefficient = OPT_channel_sourcesink[l].reactant_ID[i];
					buf_temperature += (*(T_species+coefficient+buf_start));
					counter_duoble += 1.0;
				}
				buf_temperature = buf_temperature / counter_duoble;



					if ( buf_temperature < channel[l].origin_eV[0])
					{
						chem_rate_buf[l+buf_start_channel] = channel[l].origin_rate[0];
					}
					else if (buf_temperature >= channel[l].origin_eV.back())
					{
						chem_rate_buf[l+buf_start_channel] = channel[l].origin_rate.back();
					}
					else
					{

						i = 0;
						while ( ! (buf_temperature >= channel[l].origin_eV[i] && buf_temperature < channel[l].origin_eV[i+1] ))
						{
							i+=1;
						}
						chem_rate_buf[l+buf_start_channel] =  interpolation (& channel[l].origin_eV[i]  , & channel[l].origin_eV[i+1]
										                    ,& channel[l].origin_rate[i], & channel[l].origin_rate[i+1]
										                    ,& buf_temperature);

					}

			}
		}
	}

}
/**
 * @brief Calculate total pressure and total density. (No consider gas temperature)
 *
 * @param[in]  D_species  The density of species.
 * @param[in]  T_species  The temperature of species.
 * @date 04/01/2009
 */
void   chemistry::CalTotalPressure( const double *D_species , const double* T_species)
{
	int     k ;
	int     buf_start;

	for ( int i = 0 ; i < domain_size ; ++i ) {

		buf_start 				= i*species_size ;

		/* reset */
		total_gas_density [ i ] = 0.0 ;
		total_gas_pressure[ i ] = 0.0 ;

		for ( int iSpecies =0 ; iSpecies  < species_size ; ++iSpecies  ) {

			/* ignore partial pressure of electron */
			if ( iSpecies  != species.e_id ) {
				total_gas_density[ i ] = total_gas_density[ i ] + (*( D_species + iSpecies + buf_start ) );
			}

		}//End iSpecies

		total_gas_pressure[i]=1.03558e-25*total_gas_density[i]*(*(T_species+species_size-1+buf_start)) ;

	}//End domein loop

}




/*=======================================================================================
|========================================================================================
|          After obtaining rate, It is time to estimate source and sink
|          electron energy loss, Collision frequence
|========================================================================================
|=======================================================================================*/
//cerr << "sosi size :    " << species.sourcesink.size() << endl;
//for (i=0;i<species.sourcesink.size();++i)  cerr << "in sourcesink  in global : " << species.sourcesink[i].global.size() << endl;
//cerr << "energy_loss  size :  " << energy_loss.global.size() << endl;
//cerr << "coll_frequency size :  " << coll_frequency.global.size() << endl;

/**
 * @brief After obtaining rate, It is time to estimate source and sink, electron energy loss, Collision frequence
 *
 * @param[in]  D_species  The density of species.
 * @param[in]  T_species  The temperature of species.
 * @date 05/27/2009
 */
void   chemistry::SourceSink(const double* D_species,const double* T_species)
{
	int       domain_int;
	int       i , j , k , reactant_id ;
	int       coefficient;
	int       buf_start_channel;
	int       buf_start;
	double    temp,temp_e;

	/* initial sourcesink and light power */
	for ( i = 0 ; i < sourcesink_size ; ++i ) {
		for ( j = 0 ; j < domain_size ; ++j ) {
			species.sourcesink[ i ].global[ j ] = 0.0 ;
		}//End cell loop
	}//End iSpecies

	for (i=0;i<light_size;++i) {
		for (j=0;j<domain_size;++j) {
			species.light_power[ i ].global[ j ] = 0.0 ;
		}
	}

	reactant_id=10000;

	for ( domain_int = 0 ; domain_int < domain_size ; ++domain_int ) {

		buf_start = domain_int*species_size ;

		buf_start_channel = domain_int*channel_size;

		coll_frequency.global[domain_int]=0.0;

		energy_loss.global[domain_int]=0.0;

		for ( k = 0 ; k < channel_size ; ++k ) {

			temp=chem_rate_buf[ k+buf_start_channel ] ;

			temp_e = temp ;

			for ( i = 0 ; i < OPT_channel_sourcesink[k].sink_n ; ++i ) {

				coefficient = OPT_channel_sourcesink[k].reactant_ID[i];

				if ( OPT_channel_sourcesink[k].coefficient_reactant[coefficient] == 1 ) {

					temp = temp * (*(D_species+coefficient+buf_start));

					if ( OPT_channel_sourcesink[k].IS_momentum_Xsfer == 1 ) {
						if (coefficient != species.e_id) {
							reactant_id=coefficient;
							//cout<<" Kmt =  "<<temp_e<<endl;
							temp_e=temp_e * (*(D_species+coefficient+buf_start));
							//cout<<"N_back = "<<*(D_species+coefficient+buf_start)<<endl;//BackGround NumberDensity
							//cout<<"Kmt*N_back = "<<temp_e<<endl;
						}
					}//Is momentum transfer

				}	else if (OPT_channel_sourcesink[k].coefficient_reactant[coefficient] == 2) {

					temp = temp * pow( (*(D_species+coefficient+buf_start)) , 2);
					if (OPT_channel_sourcesink[k].IS_momentum_Xsfer == 1) {

						cerr << " Big error at momentum Xsfer. Exit Chemistry-module " << endl;
						exit(1);
						if (coefficient != species.e_id) {
						    reactant_id=coefficient;
						    temp_e=temp_e * pow( (*(D_species+coefficient+buf_start)) , 2) ;
						    //cout<<"KL2: "<<*(D_species+coefficient+buf_start)<<endl;
						}
					}

				}	else if (OPT_channel_sourcesink[k].coefficient_reactant[coefficient] == 3) {

					temp = temp * pow( (*(D_species+coefficient+buf_start)) , 3);
					if ( OPT_channel_sourcesink[k].IS_momentum_Xsfer == 1 ) {
						cerr << " Big error at momentum Xsfer. Exit Chemistry-module " << endl;
						exit(1);
						if (coefficient != species.e_id){
							reactant_id=coefficient;
							temp_e=temp_e * pow( (*(D_species+coefficient+buf_start)) , 3);
							//cout<<"KL3: "<<*(D_species+coefficient+buf_start)<<endl;
						}
					}

				}//
			}//end sink n

			/*  channel Source-Sink term */
			channel[k].Channel_SourceSink[domain_int]=temp;

			//  collisional frequency of electron  //
			if ( OPT_channel_sourcesink[k].coefficient_reactant[species.e_id] != 0 ) {
				if (OPT_channel_sourcesink[k].IS_momentum_Xsfer == 1)  coll_frequency.global[domain_int] = coll_frequency.global[domain_int] +  temp_e;
			}

			// electron energy loss //
		        if ( OPT_channel_sourcesink[k].coefficient_reactant[species.e_id] != 0 || OPT_channel_sourcesink[k].coefficient_product[species.e_id] != 0) {
				//cerr << "reactant_id" << species.name[reactant_id] << endl;
				//cerr << "klkl;lkljc,moool   " <<  energy_loss.global[domain_int] << endl;
				//
				//elastic energy loss term
				if (OPT_channel_sourcesink[k].IS_momentum_Xsfer != 1)
				{
					energy_loss.global[domain_int] = energy_loss.global[domain_int] + (temp * e_loss_energy[k]);
				//cout << k <<"  e loss  energy  " << (temp * e_loss_energy[k])<< endl;
				}
				else
				{
				    energy_loss.global[domain_int] = energy_loss.global[domain_int]
				    + 3.0 * species.mass[species.e_id]/species.mass[reactant_id] * (*(D_species+species.e_id+buf_start))*temp_e*(*(T_species+species.e_id+buf_start));
				    //cout<<"Line:1585"<<endl;
				    //cout<<" m_e = "<<species.mass[species.e_id]<<endl;
				    //cout<<" m_N = "<<species.mass[reactant_id]<<endl;
				    //cout<<"N_e = "<<(*(D_species+species.e_id+buf_start))<<endl;
				    //cout<<"Te = "<<(*(T_species+species.e_id+buf_start))<<endl;
				    //cout<<"kmt*n_back = "<<temp_e<<endl;
				    //cout<<"elastic = "<<3.0 * species.mass[species.e_id]/species.mass[reactant_id] * (*(D_species+species.e_id+buf_start))*temp_e*(*(T_species+species.e_id+buf_start))<<endl;
				    //cout<<"elastic_myself = "<<3.0*(species.mass[species.e_id]/species.mass[reactant_id])*3.218035E21*3.E-14*(*(D_species+species.e_id+buf_start))*1.602E-19*2.0<<endl;
				    //exit(1) ;
				    //cerr << "jcjdd    "  << k << endl;
				}
				//cerr << k << "  " << (temp * channel[k].threshold) << endl;
					//cerr << "elastic  "<< k << "  "  <<3.0 * species.mass[species.e_id]/species.mass[reactant_id] * (*(D_species[species.e_id]+domain_int))*temp_e*(*(T_species[species.e_id]+domain_int)) << endl;

		}


			/* sourcesink = source - sink */
			for ( i=0 ; i<species_size ;++i ) {
				if ( net_sourcesink[k].sign[i] != 0 ) {

					species.sourcesink[i].global[domain_int]
					= species.sourcesink[i].global[domain_int] + temp * net_sourcesink[k].sign[i];
					cout<<"iSpecies: "<<i<<"\t"<<"Source: "<<species.sourcesink[i].global[domain_int]<<endl;

				}//End net rate
			}//End iSpecies

			if (channel[k].emit_id >= 0) {
				species.light_power[channel[k].emit_id].global[domain_int] = species.light_power[channel[k].emit_id].global[domain_int] + (temp * 1.9864748e-16 / species.wavelength[channel[k].emit_id]);
				//cerr << "  emit id  " << channel[k].emit_id<< endl;
				//cerr << k <<"    "<<temp * 1.9864748e-16 / species.wavelength[channel[k].emit_id] << endl;
			}
		}
	}//End cell loop
}



void    chemistry::UpdateElectronDiffMob(const double* T_species,const double* D_species){
	//Mark By KL:int      j,k;
	int      i,buf_start;
	//double   total_pressure;
	//double   totalGasDensity;
	// Mark BY KL :double   TeTemp;
	// Mark BY KL :int      table_element;

	for (i=0;i<domain_size;++i) {
	        buf_start = i*species_size;

			/*===========================
			|   mobility of electron    |
			===========================*/
			switch (species.mob_type[species.e_id])
			{
			    case 0:
				    species.mobility[species.e_id].global[i] = species.mob_constant[species.e_id];
			    break;
			    case 1:
				    species.mobility[species.e_id].global[i] = species.self_mob_constatn[species.e_id]/total_gas_pressure[i];
			    break;
			    case 2:
				    species.mobility[species.e_id].global[i] = 1.6021764e-19/species.mass[species.e_id]/coll_frequency.global[i];
			    break;
			    case 3:
				    cerr << "  The mobility_type of e, we do not have frmular. Please choose another." << endl;
				    exit(1);
			    break;
			    case 4:
				    cerr << "  The mobility_type of e, we do not have data from Ward. Please choose another." << endl;
				    exit(1);
			    break;
			    case 5:
				    species.mobility[species.e_id].global[i] = mob_buf[i];
			    break;
			    case 6:
				    species.mobility[species.e_id].global[i] = mob_buf[i];
			    break;
			    case 7:
				    cerr << " The mobility_type of e, It does not be prepared.  " << endl;
			    break;
			    default:
				    cerr << "  The mobility_type of e you choose is illegal. Please check!!!!  " << endl;
				    exit(1);
			    break;
			}
//cerr << "mobility:   " << species.mobility[species.e_id].global[i] << endl;
			/*============================
			|   diffusion of electron    |
			============================*/
			switch (species.diff_type[species.e_id])
			{
			    case 0:
				    species.diffusion[species.e_id].global[i]=species.diff_constant[species.e_id];
			    break;
			    case 1:
				    species.diffusion[species.e_id].global[i]=species.self_diff_constatn[species.e_id]/total_gas_pressure[i];
			    break;
			    case 2:
				    species.diffusion[species.e_id].global[i]=(*(T_species+species.e_id+buf_start))*1.6021764e-19/species.mass[species.e_id]/coll_frequency.global[i];
			    break;
			    case 3:
				    cerr << "  The diffusion_type of e you choose is 3, we do not have formula. Please choose another." << endl;
				    exit(1);
			    break;
			    case 4:
				    cerr << "  The diffusion_type of e you choose is 4, we do not have data. Please choose another." << endl;
				    exit(1);
			    break;
			    case 5:

				    species.diffusion[species.e_id].global[i] = diff_buf[i] ;
			    break;
			    case 6:

				    species.diffusion[species.e_id].global[i] = diff_buf[i] ;

			    break;
			    case 7:
				    cerr << "  The diffusion_type of e you choose is 7, we do not prepare.   " << endl;
			    break;
			    default:
				    cerr << " The diffusion_type of e you choose is illegal. Please check!!!! " << endl;
				    exit(1);
			    break;
			}
//cerr << "ans:  " << species.diffusion[species.e_id].global[i] << endl;
	}
}


void    chemistry::UpdateNeutralDiffMob(const double* T_species,const double* D_species){
	int               i;
	int               k;
	int               kk;
	int               buf_start;
	//double            total_pressure;
	double            total_pressure_atm;
	//double            total_density;
	double            collisional_integral;
	double            dimless_T;
	double            temp_T;
	double            buffer;
	double            SaveTimeBuf;
	vector<double>    D_ab;
	vector<double>    per_pressure;
	SaveTimeBuf = 0.0;

	per_pressure.resize(species_size,0.0);
	D_ab.resize(species_size,0.0);
	for (i=0;i<domain_size;++i)
	{
	    buf_start = i*species_size;

	    for (k=0;k<species_size;++k) {
		    if (k != species.e_id) {
			    //per_pressure[k] = 1.03558e-25*(*(D_species+k+buf_start))*T_Gas ;
			    per_pressure[k] = 1.03558e-25*(*(D_species+k+buf_start)) * (*(T_species+species_size-1+buf_start));
		    }
		    else{per_pressure[k] = 0.0;}
	    }
	    total_pressure_atm = total_gas_pressure[i] * 1.3158e-3;


		    for (k=0;k<neutral_size;++k)
		    {
                        /*============================
			|   mobility of neutrals     |
			============================*/
			species.mobility[species.neutral_id[k]].global[i]=0.0;
			/*============================
			|   diffusion of electron    |
			============================*/
			switch (species.diff_type[species.neutral_id[k]])
			{
			    case 0:
				    species.diffusion[species.neutral_id[k]].global[i]=species.diff_constant[species.neutral_id[k]];
				    //cerr << "diff_constant  " << species.diff_constant[species.neutral_id[k]]<<endl;
			    break;
			    case 1:
				    species.diffusion[species.neutral_id[k]].global[i]=species.viscosity[species.neutral_id[k]]*species.self_diff_constatn[species.neutral_id[k]]/(total_gas_density[i]*species.mass[species.neutral_id[k]]);
			    break;
			    case 2:
				    cerr << "  The neutral species do not have mobility go obtain diffusion   " << endl;
				    exit(1);
			    break;
			    case 3:
				    //D_ab.clear();
				    for (kk=0 ; kk<species_size ; ++kk) {
					    if (kk != species.e_id) {
						    //temp_T = ( (*(T_species+species.neutral_id[k]+buf_start)) + (*(T_species+kk+buf_start)) )/2.0;
						    temp_T =  (*(T_species+species_size-1+buf_start));
						    //if ( abs(temp_T - SaveTimeBuf) > 1.0 ){
						    dimless_T = temp_T / sqrt(species.LJ[species.neutral_id[k]]*species.LJ[kk]);
						    collisional_integral = 1.06036 / pow(dimless_T,0.1561)
								         + 0.193   / exp(0.47635*dimless_T)
									 + 1.03587 / exp(1.52996*dimless_T)
									 + 1.76474 / exp(3.89411*dimless_T) ;
						    //SaveTimeBuf = temp_T;
						    //}
						    D_ab[kk]=(1.8583e-7 * sqrt( pow(temp_T,3.0) * (1.0/species.amu[species.neutral_id[k]] + 1.0/species.amu[kk]) )
						                   / (total_pressure_atm * pow ( (species.binary_diameter[species.neutral_id[k]]+species.binary_diameter[kk])/2.0 , 2.0 )
								   * collisional_integral)  ) ;
					    }
					    else {D_ab[kk]=0.0;}
				    }
				    buffer=0.0;
				    for (kk=0;kk<species_size;++kk) {
					    if (kk != species.e_id) {
						    buffer= buffer + per_pressure[kk] / D_ab[kk];
					    }
				    }
				    species.diffusion[species.neutral_id[k]].global[i]=total_gas_pressure[i] / buffer;
				    //cerr << "    D_ab size            "   << D_ab.size() << endl;
				    //cerr << "   pressure in atm:    " << total_pressure_atm << endl;
				    //cerr << "  diff  size   "<<species.diffusion[species.neutral_id[k]].global.size() << endl;
			    break;
			    case 4:
				    cerr << "  The diffusion_type of neutral you choose is 4, we do not have data. Please choose another." << endl;
				    exit(1);
			    break;
			    case 5:
				    cerr << "  The diffusion_type of neutral you choose is 5, we do not have data from BOLSIG." << endl;
				    cerr << "  Please choose another." << endl;
				    exit(1);
			    break;
			    case 6:
				    cerr << "  The diffusion_type of neutral you choose is 6, I do not design yet." << endl;
				    cerr << "  Please choose another or tell me." << endl;
				    exit(1);
			    break;
			    case 7:
				    cerr << "  The diffusion_type of neutral you choose is 6, I do not design yet." << endl;
				    cerr << "  Please choose another " << endl;
				    exit(1);
			    break;
			    default:
				    cerr << " The diffusion_type you choose is illegal. Please check!!!! " << endl;
				    exit(1);
			    break;
			}
		    }
	}
}


void    chemistry::UpdateIonDiffMob(const double* T_species,const double* D_species,const double * E_field){
    //Mark By KL: int               j ;
    int               i, kk, k;
    int               buf_start;
    //double            total_pressure;
    double            temp_T;
    double            buffer;

    vector<double>    mu_ij;
    vector<double>    per_pressure;
    noble_drift       gas_transpotr;
    double            EoverN;

    per_pressure.resize(species_size,0.0);
    mu_ij.resize(species_size,0.0);
	for (i=0;i<domain_size;++i)
	{
		buf_start = i*species_size;
		for (k=0;k<species_size;++k) {
			if (k != species.e_id) {
				per_pressure[k]=( 1.03558e-25*(*(D_species+k+buf_start))*T_Gas );
			}
			else { per_pressure[k]=0.0; }
		}

		for (k=0;k<ion_size;++k)
		{

		        /*======================
			|   mobility of ion    |
			======================*/
			switch (species.mob_type[species.ion_id[k]]){
			    case 0:
				    species.mobility[species.ion_id[k]].global[i]=species.mob_constant[species.ion_id[k]];
			    break;
			    case 1:
				    species.mobility[species.ion_id[k]].global[i]=species.self_mob_constatn[species.ion_id[k]]/total_gas_pressure[i];
			    break;
			    case 2:
				    cerr << "The mobility of ions can not be estimated from diffusion by Einstein relation" << endl;
			    break;
			    case 3:
				    //mu_ij.clear();
				    for (kk=0;kk<species_size;++kk) {
					    if (kk != species.e_id) {
						    temp_T = ( (*(T_species+species.ion_id[k]+buf_start)) + (*(T_species+kk+buf_start)) )/2.0;
						    mu_ij[kk]=( 0.0038553277 * temp_T /total_gas_pressure[i]
								     / sqrt(species.reduce_amu[species.ion_id[k]*species_size+kk]*species.polarizability[kk])) ;
					    }
					    else {mu_ij[kk]=0.0;}
//cerr << "ckljcjuiek       cklcjkjklje  :"<<mu_ij[kk] << endl;
				    }
				    buffer=0.0;
				    for (kk=0;kk<species_size;++kk)
				    {
					    if (kk != species.e_id)
					    {
						    buffer= buffer + per_pressure[kk] / mu_ij[kk];
					    }
				    }
				    species.mobility[species.ion_id[k]].global[i]=total_gas_pressure[i] / buffer;
			    break;
			    case 4:
				    species.mobility[species.ion_id[k]].global[i]=gas_transpotr.mobility_Ward(species.amu[species.ion_id[k]],*(E_field+i)/total_gas_pressure[i]) / *(E_field+i);
			    break;
			    case 5:
				    cerr << "  The mobility_type of ions, we do not have data from BOLSIG. Please choose another." << endl;
				    exit(1);
			    break;
			    case 6:
				    cerr << "  The diffusion_type of neutral you choose is 6, I do not design yet." << endl;
				    cerr << "  Please choose another or tell me." << endl;
				    exit(1);
			    break;
			    case 7:
				    EoverN=abs(*(E_field+i))/total_gas_density[i];
				    if (species.amu[species.ion_id[k]] >3.5 && species.amu[species.ion_id[k]] <4.5)
				    {
					    if (EoverN<10.0)   {species.mobility[species.ion_id[k]].global[i]=10.2*(760.0/total_gas_pressure[i])*(T_Gas/273.16);}
					    else if (EoverN>=10.0 && EoverN<50.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(10.0,50.0,10.2*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,8.97*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=50.0 && EoverN<100.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(50.0,100.0,8.97*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,7.67*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=100.0 && EoverN<200.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(100.0,200.0,7.67*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,6.12*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=200.0 && EoverN<300.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(200.0,300.0,6.12*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,5.19*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=300.0 && EoverN<400.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(300.0,400.0,5.19*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,4.58*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=400.0 && EoverN<500.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(400.0,500.0,4.58*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,4.17*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=500.0 && EoverN<600.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(500.0,600.0,4.17*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,3.81*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=600.0){ species.mobility[species.ion_id[k]].global[i]=3.81*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001;}
					    //cout << "He + mob : "<<species.mobility[species.ion_id[k]].global[i]<< endl;
				    }
				    else if (species.amu[species.ion_id[k]] >7.5 && species.amu[species.ion_id[k]] <9.5)
				    {
					    if (EoverN<10.0)   {species.mobility[species.ion_id[k]].global[i]=16.9*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001;}
					    else if (EoverN>=10.0 && EoverN<20.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(10.0,20.0,16.9*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,17.7*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=20.0 && EoverN<22.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(20.0,22.0,17.7*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,18.0*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=22.0 && EoverN<24.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(22.0,24.0,18.0*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,18.3*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,EoverN);}
					    else if (EoverN>=24.0 && EoverN<100.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(24.0,100.0,18.3*(760.0/total_gas_pressure[i])*(T_Gas/273.16)*0.0001,0.00379265,EoverN);}
					    else if (EoverN>=100.0 && EoverN<300.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(100.0,300.0,0.00379265,0.00672134,EoverN);}
					    else if (EoverN>=300.0 && EoverN<500.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(300.0,500.0,0.00672134,0.00965002,EoverN);}
					    else if (EoverN>=500.0){ species.mobility[species.ion_id[k]].global[i]=0.00965002;}
					    //cout << "He222 + mob : "<<species.mobility[species.ion_id[k]].global[i]<< endl;
				    }
				    else if (species.amu[species.ion_id[k]] >130.0 && species.amu[species.ion_id[k]] <135.0)
				    {
					    if (EoverN<50.0)   {species.mobility[species.ion_id[k]].global[i]=0.0000508;}
					    else if (EoverN>= 50.0 && EoverN<60.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(50.0,60.0,0.0000508,0.0000502,EoverN);}
					    else if (EoverN>= 60.0 && EoverN<70.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(60.0,70.0,0.0000502,0.0000495,EoverN);}
					    else if (EoverN>= 70.0 && EoverN<80.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(70.0,80.0,0.0000495,0.0000489,EoverN);}
					    else if (EoverN>= 80.0 && EoverN<100.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(80.0,100.0,0.0000489,0.0000476,EoverN);}
					    else if (EoverN>= 100.0 && EoverN<220.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(100.0,220.0,0.0000476,0.0000409,EoverN);}
					    else if (EoverN>= 220.0 && EoverN<300.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(220.0,300.0,0.0000409,0.0000376,EoverN);}
					    else if (EoverN>= 300.0 && EoverN<400.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(300.0,400.0,0.0000376,0.0000344,EoverN);}
					    else if (EoverN>= 400.0 && EoverN<500.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(400.0,500.0,0.0000344,0.0000320,EoverN);}
					    else if (EoverN>= 500.0 && EoverN<600.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(500.0,600.0,0.0000320,0.0000300,EoverN);}
					    else if (EoverN>= 600.0 && EoverN<700.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(600.0,700.0,0.0000300,0.0000284,EoverN);}
					    else if (EoverN>= 700.0 && EoverN<800.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(700.0,800.0,0.0000284,0.000027,EoverN);}
					    else if (EoverN>= 800.0 && EoverN<900.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(800.0,900.0,0.0000270,0.0000258,EoverN);}
					    else if (EoverN>= 900.0 && EoverN<1000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(900.0,1000.0,0.0000258,0.0000247,EoverN);}
					    else if (EoverN>= 1000.0 && EoverN<2000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(1000.0,2000.0,0.0000247,0.0000186,EoverN);}
					    else if (EoverN>= 2000.0 && EoverN<3000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(2000.0,3000.0,0.0000186,0.0000156,EoverN);}
					    else if (EoverN>= 3000.0 && EoverN<4000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(3000.0,4000.0,0.0000156,0.0000137,EoverN);}
					    else if (EoverN>= 4000.0 && EoverN<5000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(4000.0,5000.0,0.0000137,0.0000124,EoverN);}
					    else if (EoverN>= 5000.0 && EoverN<6000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(5000.0,6000.0,0.0000124,0.0000115,EoverN);}
					    else if (EoverN>= 6000.0 && EoverN<7000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(6000.0,7000.0,0.0000115,0.0000107,EoverN);}
					    else if (EoverN>= 7000.0 && EoverN<8000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(7000.0,8000.0,0.0000107,0.0000101,EoverN);}
					    else if (EoverN>= 8000.0 && EoverN<9000.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(8000.0,9000.0,0.0000101,0.0000095,EoverN);}
					    else if (EoverN>= 9000.0 && EoverN<10000.0){species.mobility[species.ion_id[k]].global[i]=mob_interpolation(9000.0,10000.0,0.0000095,0.0000091,EoverN);}
					    else if (EoverN>= 10000.0){ species.mobility[species.ion_id[k]].global[i]=0.0000091;}
					    //cout << "He222 + mob : "<<species.mobility[species.ion_id[k]].global[i]<< endl;


				    }
				    else if (species.amu[species.ion_id[k]] >260.0 && species.amu[species.ion_id[k]] <265.0)
				    {
					    if (EoverN<80.0)   {species.mobility[species.ion_id[k]].global[i]=0.0000605;}
					    else if (EoverN>= 80.0 && EoverN<100.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(80.0,100.0,0.0000605,0.0000603,EoverN);}
					    else if (EoverN>= 100.0 && EoverN<120.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(100.0,120.0,0.0000603,0.0000601,EoverN);}
					    else if (EoverN>= 120.0 && EoverN<150.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(120.0,150.0,0.0000601,0.0000595,EoverN);}
					    else if (EoverN>= 150.0 && EoverN<190.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(150.0,190.0,0.0000595,0.0000581,EoverN);}
					    else if (EoverN>= 190.0 && EoverN<220.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(190.0,220.0,0.0000581,0.0000567,EoverN);}
					    else if (EoverN>= 220.0 && EoverN<250.0){	species.mobility[species.ion_id[k]].global[i]=mob_interpolation(220.0,250.0,0.0000567,0.0000551,EoverN);}
					    else if (EoverN>= 250.0){ species.mobility[species.ion_id[k]].global[i]=0.0000551;}


				    }
				    else
				    {
					    cout << " Error    Exit   " << endl;
					    exit(1);
				    }
			    break;
			    default:
				    cerr << "  The mobility_type of ions you choose is illegal. Please check!!!!  " << endl;
				    exit(1);
			    break;
			}
//cerr << "mobility:   " << species.mobility[species.ion_id[k]].global[i] << endl;
			/*======================
			|    diffusion of ion  |
			======================*/
			switch (species.diff_type[species.ion_id[k]])
			{
			    case 0:
				    species.diffusion[species.ion_id[k]].global[i]=species.diff_constant[species.ion_id[k]];
				    //cerr << "diff_constant  " << species.diff_constant[species.ion_id[k]]<<endl;
			    break;
			    case 1:
				    species.diffusion[species.ion_id[k]].global[i]=species.self_diff_constatn[species.ion_id[k]]/total_gas_pressure[i];
				    //cerr << "constant  " << species.self_diff_constatn[species.ion_id[k]] << endl;
			    break;
			    case 2:
				    species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*species.mobility[species.ion_id[k]].global[i];
			    break;
			    case 3:
				    cerr << "  The diffusion of ions can not be estimated from Chapman-Enskog eq. " << endl;
				    exit( 1 );
			    break;
			    case 4:
				    if (species.amu[species.ion_id[k]] >= 38.0 && species.amu[species.ion_id[k]] <= 42.0) // Ar
				    {
					    cout << "I do not prepare Ar by Ward formular " << endl;
				    }
				    else if (species.amu[species.ion_id[k]] >= 129.0 && species.amu[species.ion_id[k]] <= 134.0) //Xe
				    {
					    // Formula is not simplfy
					    //((*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*5.8e-5)
		       //+ 4.72* mass  /(3.0*2.908*1.602177e-19) * pow(species.mobility[species.ion_id[k]].global[i],3)*pow( *(E_field+i),2)* 1.0e6;
					    //species.diffusion[species.ion_id[k]].global[i]=((*(T_species+species.ion_id[k]+buf_start))*4.9980837e-9)
		       species.diffusion[species.ion_id[k]].global[i]= 4.72* 2.1809E-022  /(1.397739e-18) * pow(species.mobility[species.ion_id[k]].global[i],3)*pow( *(E_field+i),2);
		    		    //cout << "Xe  "<< species.mass[species.ion_id[k]] << endl;

				    }
				    else if (species.amu[species.ion_id[k]] >= 259.0 && species.amu[species.ion_id[k]] <= 265.0) //Xe2
				    {
					    //species.mobility[species.ion_id[k]].global[i]=1.5 * gas_transpotr.mobility_Ward(131.293,*(E_field+i)/total_gas_pressure[i]) / *(E_field+i);
					    //species.diffusion[species.ion_id[k]].global[i]=((*(T_species+species.ion_id[k]+buf_start))*6.8e-9)
		       species.diffusion[species.ion_id[k]].global[i]= 8.44* 2.1809E-022  /(2.3148e-18) * pow(species.mobility[species.ion_id[k]].global[i],3)*pow( *(E_field+i),2) ;
		    		    //cout << "Xe2  "<<species.mass[species.ion_id[k]] << endl;
				    }
				    else
				    {
					    cerr << "Can not identitfy your particles " << endl;
					    cerr << " Stop it now!!!!!!!!!!!!" << endl;
					    exit(1);
				    }


			    break;
			    case 5:
				    cerr << "  The mobility_type of ions, we do not have data from BOLSIG. Please choose another." << endl;
				    exit(1);
			    break;
			    case 6:
				    cerr << "  The diffusion_type of neutral you choose is 6, I do not design yet." << endl;
				    cerr << "  Please choose another or tell me." << endl;
				    exit(1);
			    break;
			    case 7:
				    if ( abs(*(E_field+i)) != 0.0 )
				    {
					    if (species.amu[species.ion_id[k]] >3.5 && species.amu[species.ion_id[k]] <4.5)
					    {
					        species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*9.48e-4
						+ 1.0e-4*4.72/7.824*6.6464758e-24/(1.6e-19* abs(*(E_field+i)) *1.0e5)*pow(species.mobility[species.ion_id[k]].global[i]*abs(*(E_field+i))*100.0,3);
						//cout << "  mob   " << species.mobility[species.ion_id[k]].global[i] << endl;
						//cout << "  Einstein   " << (*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*9.48e-4 << endl;
						//cout << "diff  "<<species.diffusion[species.ion_id[k]].global[i] << endl;
					    }
					    else if (species.amu[species.ion_id[k]] >7.5 && species.amu[species.ion_id[k]] <9.5)
					    {
						species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*17.4e-4
						+ 1.0e-4*8.44/14.448*6.6464758e-24/(1.6e-19* abs(*(E_field+i)) *1.0e5)*pow(species.mobility[species.ion_id[k]].global[i]*abs(*(E_field+i))*100.0,3);
						//cout << "   mob   " << species.mobility[species.ion_id[k]].global[i] << endl;
						//cout << "  Einstein  " << (*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*17.4e-4 << endl;
						//cout << "diff  "<<species.diffusion[species.ion_id[k]].global[i] << endl;
					    }
					    else if (species.amu[species.ion_id[k]] >130.0 && species.amu[species.ion_id[k]] <135.0)
					    {
						    species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*0.58e-4
						+ 1.0e-4*4.72/7.824*2.1801213e-22/(1.6e-19* abs(*(E_field+i)) *1.0e5)*pow(species.mobility[species.ion_id[k]].global[i]*abs(*(E_field+i))*100.0,3);

					    }
					    else if (species.amu[species.ion_id[k]] >260.0 && species.amu[species.ion_id[k]] <265.0)
					    {
						    species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*0.79e-4
						+ 1.0e-4*8.44/14.448*4.3602426e-22/(1.6e-19* abs(*(E_field+i)) *1.0e5)*pow(species.mobility[species.ion_id[k]].global[i]*abs(*(E_field+i))*100.0,3);

					    }

					    else
					    {
						cerr << "   Errror at mob 7  "<< endl;
					    }
				    }
				    else
				    {
					    if (species.amu[species.ion_id[k]] >3.5 && species.amu[species.ion_id[k]] <4.5)
					    {
					        species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*9.48e-4 ;
						//cout << "  mob   " << species.mobility[species.ion_id[k]].global[i] << endl;
						//cout << "  Einstein   " << (*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*9.48e-4 << endl;
						//cout << "diff  "<<species.diffusion[species.ion_id[k]].global[i] << endl;
					    }
					    else if (species.amu[species.ion_id[k]] >7.5 && species.amu[species.ion_id[k]] <9.5)
					    {
						species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*17.4e-4 ;
						//cout << "   mob   " << species.mobility[species.ion_id[k]].global[i] << endl;
						//cout << "  Einstein  " << (*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*17.4e-4 << endl;
						//cout << "diff  "<<species.diffusion[species.ion_id[k]].global[i] << endl;
					    }
					    else if (species.amu[species.ion_id[k]] >130.0 && species.amu[species.ion_id[k]] <135.0)
					    {
						species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*0.58e-4;
					    }
					    else if (species.amu[species.ion_id[k]] >260.0 && species.amu[species.ion_id[k]] <265.0)
					    {
						species.diffusion[species.ion_id[k]].global[i]=(*(T_species+species.ion_id[k]+buf_start))*8.61738569225667e-05*0.79e-4;
					    }
					    else
					    {
						    cerr << "   Errror at mob 7  "<< endl;
					    }

				    }
			    break;
			    default:
				    cerr << " The diffusion_type you choose is illegal. Please check!!!! " << endl;
				    exit( 1 );
			    break;
			}
		}
	}
}








double     noble_drift::mobility_Ward(double  noble_amu , double W)
{
	mu_plus = 0.0;
	C = 0.0;
	k_plus = 0.0;
	D = 0.0;
	if (noble_amu > 12.0 && noble_amu < 30.0)  // Ne
	{
		cerr << "Sorry, we do not have   Ne   data now!!!!!!" << endl;
		exit(1);
	}
	else if (noble_amu > 30.0 && noble_amu < 60.0) // Ar
	{
		if (W <= 6000.0)
		{
			mu_plus = 0.1;
			C = 2.22e-5;
			return v_fun_1(W);
		}
		else
		{
			k_plus = 8.25;
			D = 8.625e4;
			return v_fun_2(W);
		}

	}
	else if (noble_amu > 60.0 && noble_amu < 110.0) // Kr
	{
		cerr << "Sorry, we do not have   Kr   data now!!!!!!" << endl;
		exit(1);
	}
	else if (noble_amu > 110.0 && noble_amu < 175.0) // Xe
	{
		if (W <= 10000.0)
		{
			mu_plus = 0.04;
			C = 2.25e-5;
			return v_fun_1(W);
		}
		else
		{
			k_plus = 4.0;
			D = 2.25e5;
			return v_fun_2(W);
		}
	}
	else if (noble_amu >= 259.0 && noble_amu <= 265.0) // Xe
	{
		if (W <= 10000.0)
		{
			mu_plus = 0.04;
			C = 2.25e-5;
			return 1.5 * v_fun_1(W);
		}
		else
		{
			k_plus = 4.0;
			D = 2.25e5;
			return 1.5 * v_fun_2(W);
		}
	}
	else{
		cerr << "Sorry, we do not have data now!!!!!!" << endl;
		exit(1);
	}


}
double   noble_drift::v_fun_1(double W)
{
	return mu_plus*W*(1.0-C * W);
}
double   noble_drift::v_fun_2(double W)
{
	return k_plus*sqrt(W)*(1.0-D*pow(1.0/W,1.5));
}


//double   noble_drift::diffusion_McDaniel(double K_0,double mass,double drift_v,double E,double Tem)
//{
//		return (Tem*8.61738569225667e-05*K_0)
//		       + (1.0+3.72)* mass * drift_v *drift_v*drift_v /(3.0*(1.0+1.908)*1.602177e-19*E);
//}


















void  chemistry::energy_loss_coll_fre_diff_mob_init()
{
	int i;
	for (i=0;i<domain_size;++i) energy_loss.global.push_back(0.0);
	for (i=0;i<domain_size;++i) coll_frequency.global.push_back(0.0);
	for (i=0;i<species.name.size();++i) species.sourcesink.push_back(energy_loss);
	for (i=0;i<species.name.size();++i) species.diffusion.push_back(energy_loss);
	for (i=0;i<species.name.size();++i) species.mobility.push_back(energy_loss);
	for (i=0;i<species.light_name.size();++i)  species.light_power.push_back(energy_loss);
}


//optimization of source-sink 2009-01-07
//no more apply do-loop to scan source and sink term
void   chemistry::OptimizationSourceSink()
{
	int    k,i,j;
	int    l;
	int    m;
	int    temp;
	int    buf_start;
	int    buf_start_channel;
	int    counter;
	int    domain_int;
  	  int    i_Rank;
	//vector<int>   BUFFER;

 	   MPI_Comm_rank(MPI_COMM_WORLD, &i_Rank);

	for (k=0;k<channel.size();++k)
	{
		for (i=0;i< species.name.size();++i)
		{
			if (channel[k].sink[i] != 0)
			{
				channel[k].sink_ID.push_back(i);
			}
		}
		for (i=0;i< species.name.size();++i)
		{
			if (channel[k].source[i] != 0)
			{
				channel[k].source_ID.push_back(i);
			}
		}
	}
	for (k=0;k<channel.size();++k)
	{
		temp=0;
		channel[k].momentum_reactions = 0;
		for (i=0;i< species.name.size();++i)
		{
			temp=temp+abs(channel[k].source[i]-channel[k].sink[i]);
		}
		if (temp == 0)
		{
			channel[k].momentum_reactions = 1;
		}
		else
		{
			channel[k].momentum_reactions = 0;
		}
	}


	/*
	for (k=0;k<channel.size();++k)
	{
		cout << "Reaction  " << k << ":  ";
		for (i=0;i<channel[k].sink_ID.size();++i){
		    cout << channel[k].sink[channel[k].sink_ID[i]]<< " " <<species.name[channel[k].sink_ID[i]] << " + ";
		}
		cout << " => " ;
		for (i=0;i<channel[k].source_ID.size();++i){
		    cout << channel[k].source[channel[k].source_ID[i]]<<" "<<species.name[channel[k].source_ID[i]] << " + ";
		}
		cout << endl;
	}
	*/

	species_size=species.name.size();
	channel_size=channel.size();
	ion_size=species.ion_id.size();
	neutral_size=species.neutral_id.size();
	sourcesink_size=species.sourcesink.size();
	light_size=species.light_name.size();

	OPT_channel_sourcesink = new single_sourcesink[channel_size];

	for (k=0;k<channel_size;++k)
	{
		OPT_channel_sourcesink[k].IS_momentum_Xsfer = channel[k].momentum_reactions;
		OPT_channel_sourcesink[k].sink_n   = channel[k].sink_ID.size();
		OPT_channel_sourcesink[k].source_n = channel[k].source_ID.size();
		OPT_channel_sourcesink[k].reactant_ID = new int [OPT_channel_sourcesink[k].sink_n];
		    for (i=0;i<OPT_channel_sourcesink[k].sink_n;++i)
		    {
			OPT_channel_sourcesink[k].reactant_ID[i] = channel[k].sink_ID[i] ;
		    }
		OPT_channel_sourcesink[k].product_ID  = new int [OPT_channel_sourcesink[k].source_n];
		    for (i=0;i<OPT_channel_sourcesink[k].source_n;++i)
		    {
			OPT_channel_sourcesink[k].product_ID[i] = channel[k].source_ID[i] ;
		    }
		OPT_channel_sourcesink[k].coefficient_reactant = new int [species_size];
		    for (i=0;i<species_size;++i)
		    {
			OPT_channel_sourcesink[k].coefficient_reactant[i] = channel[k].sink[i] ;
		    }
		OPT_channel_sourcesink[k].coefficient_product = new int [species_size];
		    for (i=0;i<species_size;++i)
		    {
			OPT_channel_sourcesink[k].coefficient_product[i] = channel[k].source[i] ;
		    }
	}

	// 0: electron
	// 1: ion
	// 2: neutral
	channel_species = new int[channel_size];
	for (k=0;k<channel_size;++k) channel_species[k]=100;

	int  species_int;
	for (k=0;k<channel_size;++k)
	{
		for (i=0;i<channel_size;++i)
		{
			int e_join=0;
			int ion_join=0;
			int n_join=0;
			if (OPT_channel_sourcesink[k].coefficient_reactant[species.e_id] > 0 ) {
					e_join=e_join+1;
			}
			if (e_join == 0 ) {
				for (species_int=0;species_int<ion_size;++species_int){
					if (OPT_channel_sourcesink[k].coefficient_reactant[species.ion_id[species_int]] >0 ) {
						ion_join=ion_join+1;
					}
				}
			}
			if (e_join == 0 && ion_join ==0) {
				for (species_int=0;species_int<neutral_size;++species_int){
					if (OPT_channel_sourcesink[k].coefficient_reactant[species.neutral_id[species_int]] >0 ) {
						n_join=n_join+1;
					}
				}

			}
			if (e_join != 0) channel_species[k] = 0;
			if (ion_join != 0 && e_join == 0) channel_species[k] = 1;
			if (n_join != 0 && e_join == 0 && ion_join == 0) channel_species[k] = 2;
		}
	}
	for (k=0;k<channel_size;++k)
	{
        if(i_Rank == 0){
	      cout << "reaction   " << k << "  type   " << channel_species[k] << endl;
        }
	    if (channel_species[k] >2) {cerr << "   Error  at OptimizationSourceSink() " << endl; exit(1);}
	}

	/**********************************************************************
	*  Design at 2009-04-02 for recording  source-sink term               *
	*								      *
	**********************************************************************/

	record_total_ss = new apply_channel [species_size];
	for (k=0;k<species_size;++k)
	{

		counter = 0;
		for (l=0;l<channel_size;++l)
		{
			if ((channel[l].source[k] -channel[l].sink[k]) != 0 ) counter = counter+1;
		}
		record_total_ss[k].n = counter;
		record_total_ss[k].join = new int [record_total_ss[k].n];
		record_total_ss[k].sign = new int [record_total_ss[k].n];
		counter = 0;
		for (l=0;l<channel_size;++l)
		{
			if ((channel[l].source[k] -channel[l].sink[k]) != 0 )
			{
				record_total_ss[k].join[counter] = l;
				//if (channel[l].source[k] -channel[l].sink[k] >0) record_total_ss[k].sign[counter] = 1;
				//else if (channel[l].source[k] -channel[l].sink[k] <0) record_total_ss[k].sign[counter] = -1;
				record_total_ss[k].sign[counter] = channel[l].source[k] -channel[l].sink[k];

				counter = counter+1;
			}
		}
	}

    if(i_Rank == 0){
  	  for (k=0;k<species_size;++k)
	  {
		//k=2;
		cout << species.name[k] << endl;
		for (l=0;l<record_total_ss[k].n;++l)
		{
			//cout << record_total_ss[k].join[l] <<  "\t"<< record_total_ss[k].sign[l] ;
			cout <<  "\t"<< record_total_ss[k].sign[l]  << " K_" << record_total_ss[k].join[l];
			for (i=0;i<channel[record_total_ss[k].join[l]].sink_ID.size();++i)  cout << "\t" << species.name[channel[record_total_ss[k].join[l]].sink_ID[i]] << "^" << channel[record_total_ss[k].join[l]].sink[channel[record_total_ss[k].join[l]].sink_ID[i]];
			cout << endl;
		}
		cout << endl << endl;
	  }
    }

	//exit(1);

	/**********************************************************************
	*  Design at 2009-04-02 for recording                                 *
	*								      *
	*                /  partial sourcesink  \                             *
	*		|  ------------------   |                             *
	*		 \  partial n_species   /                             *
	*        							      *
	**********************************************************************/
	for_Formfunction = new partial_s [species_size*species_size];
	for (k=0;k<species_size;++k)
	{
		buf_start=k*species_size;
		for (l=0;l<species_size;++l)
		{

			counter = 0;
			for (i=0;i<record_total_ss[l].n;++i)
			{
				temp = record_total_ss[l].join[i] ;
				for (j=0;j<channel[temp].sink_ID.size();++j)
				{
					if (channel[temp].sink_ID[j] == k) counter= counter +1;
				}
			}
			for_Formfunction[buf_start+l].n   = counter;
			for_Formfunction[buf_start+l].knn = new KNN [counter];
			counter = 0;
			for (i=0;i<record_total_ss[l].n;++i)
			{
				temp = record_total_ss[l].join[i] ;
				for (j=0;j<channel[temp].sink_ID.size();++j)
				{
					if (channel[temp].sink_ID[j] == k)
					{
						for_Formfunction[buf_start+l].knn[counter].k = temp;
						//if (channel[temp].sink[k] > 1)
						//{
							for_Formfunction[buf_start+l].knn[counter].nn_size = channel[temp].sink_ID.size();
							for_Formfunction[buf_start+l].knn[counter].sign = channel[temp].sink[k] * record_total_ss[l].sign[i];
							for_Formfunction[buf_start+l].knn[counter].nn = new int [for_Formfunction[buf_start+l].knn[counter].nn_size];
							for (m=0;m<channel[temp].sink_ID.size();++m)
							{
								for_Formfunction[buf_start+l].knn[counter].nn[m] = channel[temp].sink_ID[m];
							}
							for_Formfunction[buf_start+l].knn[counter].partial_c = channel[temp].sink[k]-1;
						//}
						//else
						//{
						/*
							for_Formfunction[buf_start+l].knn[counter].nn_size = channel[temp].sink_ID.size() - 1;
							for_Formfunction[buf_start+l].knn[counter].sign = record_total_ss[l].sign[i];
							for_Formfunction[buf_start+l].knn[counter].nn = new int [for_Formfunction[buf_start+l].knn[counter].nn_size];
							int  sss; sss=0;
							for (m=0;m<channel[temp].sink_ID.size();++m)
							{
								if (channel[temp].sink_ID[m] != k)
								{
									for_Formfunction[buf_start+l].knn[counter].nn[sss] = channel[temp].sink_ID[m];
									sss=sss+1;
								}
							}
							for_Formfunction[buf_start+l].knn[counter].partial_c = channel[temp].sink[k]-1;
						*/
						//}

						counter= counter +1;
					}
				}
			}
		}
	}



	cerr << "djkdjklddjk" << endl;
	//exit(1);

	//j=2*species_size;
	//j=4*species_size+6;
/*
	if(i_Rank == 0){
		  int kk;
		  for (kk=0;kk<species_size;++kk){

			j=kk*species_size+3 ;
			cout << "partieal with "<< species.name[kk] << endl ;

			for (k=0;k<for_Formfunction[j].n;++k){

				cout <<"\t"<<for_Formfunction[j].knn[k].sign << " K_" <<for_Formfunction[j].knn[k].k ;
		      		for (i=0;i<for_Formfunction[j].knn[k].nn_size;++i){

			    	if (for_Formfunction[j].knn[k].nn[i]==kk){

			    		if (for_Formfunction[j].knn[k].partial_c>0)
			    			cout <<" "<< species.name[for_Formfunction[j].knn[k].nn[i]]<<"^" <<for_Formfunction[j].knn[k].partial_c;
			    	}else{

			    		  cout <<" "<< species.name[for_Formfunction[j].knn[k].nn[i]] << "^" <<OPT_channel_sourcesink[for_Formfunction[j].knn[k].k].coefficient_reactant[for_Formfunction[j].knn[k].nn[i]] ;
				}
			}
		      	//cout << "\t"<<for_Formfunction[j].knn[k].partial_c<<"^" << ;
			cout << endl ;
		}
		cout << endl;

		}
	}//End  i_Rank //exit(1);
	2015/03/25 close
*/


	/**********************************************************************
	*  Design at 2009-04-14 for FormFunction                              *
	*								      *
	*                /  partial e_loss      \                             *
	*		|  ------------------   |                             *
	*		 \  partial n_species   /                             *
	* e_channel_ID: 0=no electron loss				      *
	*               1=momentum Xsfer                                      *
	*               2=electron lss                                        *
	**********************************************************************/

	e_loss_energy = new double [channel_size];
	e_channel_ID  = new int [channel_size];
	for (k=0;k<channel_size;++k)
	{
	    e_loss_energy[k]=channel[k].threshold;
	    e_channel_ID[k] =100;

	    if ( OPT_channel_sourcesink[k].coefficient_reactant[species.e_id] != 0 || OPT_channel_sourcesink[k].coefficient_product[species.e_id] != 0) {
			if (OPT_channel_sourcesink[k].IS_momentum_Xsfer != 1)
			{
				e_channel_ID[k] = 2;
			}
			else
			{
				e_channel_ID[k] = 1;
			}
	    }
	    else {e_channel_ID[k] = 0;}
	}

	e_loss_for_Formfunction = new partial_s [species_size];


	for (k=0;k<species_size;++k)
	{
	counter = 0;
	for (l=0;l<channel_size;++l)
	{
		if (e_channel_ID[l]!=0)
		{
			for (m=0;m<channel[l].sink_ID.size();++m)
			{
			    if (channel[l].sink_ID[m] == k)
			    {
				    counter = counter+1;
			    }
			}
		}
	}
	e_loss_for_Formfunction[k].n = counter;
	//cout << species.name[k] <<"   " <<e_loss_for_Formfunction[k].n << endl;
	}

	//exit(1);

	for (k=0;k<species_size;++k)
	{

		e_loss_for_Formfunction[k].knn = new KNN [e_loss_for_Formfunction[k].n];
		int ll;ll=0;
		for (l=0;l<channel_size;++l)
		{
			if (e_channel_ID[l] != 0)
			{
			    for (m=0;m<channel[l].sink_ID.size();++m)
			    {
				if (channel[l].sink_ID[m] == k)
				{
					e_loss_for_Formfunction[k].knn[ll].k = l;
					ll = ll +1;
				}
			    }
			}
		}

		if (e_loss_for_Formfunction[k].n != ll ) {cout << "It is not matched in e_loss_for_Formfunction. EXIT" << endl; exit(1);}
		for (l=0;l<e_loss_for_Formfunction[k].n;++l)
		{
		    temp=e_loss_for_Formfunction[k].knn[l].k ;

		    e_loss_for_Formfunction[k].knn[l].nn_size = channel[temp].sink_ID.size();
		    e_loss_for_Formfunction[k].knn[l].sign = 1;
		    e_loss_for_Formfunction[k].knn[l].nn = new int [e_loss_for_Formfunction[k].knn[l].nn_size];
				for (m=0;m<channel[temp].sink_ID.size();++m)
				{
					e_loss_for_Formfunction[k].knn[l].nn[m] = channel[temp].sink_ID[m];
					if (channel[temp].sink_ID[m] == k) e_loss_for_Formfunction[k].knn[l].sign = e_loss_for_Formfunction[k].knn[l].sign *channel[temp].sink[k];
				}
				e_loss_for_Formfunction[k].knn[l].partial_c = channel[temp].sink[k]-1;
		}


	}

	/*
	for (k=0;k<species_size;++k)
	{
		cout << species.name[k] << endl;
		for (l=0;l<e_loss_for_Formfunction[k].n;++l)
		{
			cout << e_loss_for_Formfunction[k].knn[l].k << endl;
		}
		cout << endl;
	}
	*/

	//exit(1);




	if(i_Rank == 0){
	  for (k=0;k<species_size;++k){

	    cout << " partial with  "<<species.name[k] << endl;
	    //cout << e_loss_for_Formfunction[k].n << endl;
	    for (l=0;l<e_loss_for_Formfunction[k].n;++l){
	  	  //cout << "test " << endl;
		  cout << e_loss_for_Formfunction[k].knn[l].sign ;
		  cout << " K_"<<e_loss_for_Formfunction[k].knn[l].k ;
		  //cout << "test " << endl;
		  for (m=0;m<e_loss_for_Formfunction[k].knn[l].nn_size;++m){
		    if (e_loss_for_Formfunction[k].knn[l].nn[m] != k){
			  cout << " "<< species.name[ e_loss_for_Formfunction[k].knn[l].nn[m] ]  <<"^" <<  OPT_channel_sourcesink[e_loss_for_Formfunction[k].knn[l].k].coefficient_reactant[e_loss_for_Formfunction[k].knn[l].nn[m]];
		    }
		    else{
			  if (e_loss_for_Formfunction[k].knn[l].partial_c>0){
			    cout << " "<< species.name[ e_loss_for_Formfunction[k].knn[l].nn[m] ]  <<"^" << e_loss_for_Formfunction[k].knn[l].partial_c ;
			  }
		    }
		  }
		  cout << endl;
	    }

    	cout << endl;

	  }
    }

	//exit(1);
	/*
	for (k=0;k<channel_size;++k)
	{
	    cout << k<<"  loss energy "<<e_loss_energy[k]<<endl;
	    cout << k <<" type  "<< e_channel_ID[k] << endl;
	}
	*/

	//exit(1);

	/**********************************************************************
	*  Design at 2009-04-16 for FormFunction                              *
	*								      *
	*                /  partial e_loss      \                             *
	*		|  ------------------   |                             *
	*		 \  partial n_species   /                             *
	* e_channel_ID: 0=no electron loss				      *
	*               1=momentum Xsfer                                      *
	*               2=electron lss                                        *
	**********************************************************************/

	vector<int>   buftemp;
		for (l=0;l<channel_size;++l)
		{
			if (channel[l].rate_type == 1 || channel[l].rate_type == 2)  buftemp.push_back(l);
		}
		e_loss_dTe_ID.n = buftemp.size();
		e_loss_dTe_ID.join = new int [e_loss_dTe_ID.n];
		//e_loss_dTe_ID.sign = new int [e_loss_dTe_ID.n];
		//counter = 0;
		for (l=0;l<buftemp.size();++l)
		{
			e_loss_dTe_ID.join[l] = buftemp[l];
			//e_loss_dTe_ID.sign[l] = channel[l].source[k] -channel[l].sink[k];
		}
		buftemp.clear();



	S_dTe_ID = new apply_channel [species_size];

	for (k=0;k<species_size;++k)
	{
		for (l=0;l<channel_size;++l)
		{
			if ((channel[l].source[k] -channel[l].sink[k]) != 0 && (channel[l].rate_type == 1 || channel[l].rate_type == 2)) buftemp.push_back(l);
		}
		S_dTe_ID[k].n = buftemp.size();
		S_dTe_ID[k].join = new int [S_dTe_ID[k].n];
		S_dTe_ID[k].sign = new int [S_dTe_ID[k].n];
		//counter = 0;
		for (l=0;l<buftemp.size();++l)
		{
			S_dTe_ID[k].join[l] = buftemp[l];
			S_dTe_ID[k].sign[l] = channel[ buftemp[l] ].source[k] -channel[ buftemp[l] ].sink[k];
		}
		buftemp.clear();
	}

    if(i_Rank == 0){
	  for (l = 0;l<e_loss_dTe_ID.n;++l)
	  {
		cout << e_loss_dTe_ID.join[l] << endl;
	  }
    }

    if(i_Rank == 0){
	  for (k=0;k<species_size;++k){
		cout << species.name[k] << endl;
		for (l = 0;l<S_dTe_ID[k].n;++l){
			cout << "reaction  "<<S_dTe_ID[k].join[l] << endl;
		}
		cout << endl;
	  }
    }


	/**********************************************************************
	*  Design at 2009-04-20 for rate constant (not table)                 *
	*								      *
	**********************************************************************/

	chem_rate_buf = new double [channel_size*domain_size];
	for (k=0;k<domain_size;++k)
	{
		buf_start=k*channel_size;
		for (j=0;j<channel_size;++j) chem_rate_buf[j+buf_start]=0.0;
	}
	chem_rate_perturb_buf = new double [channel_size*domain_size];
	for (k=0;k<domain_size;++k)
	{
		buf_start=k*channel_size;
		for (j=0;j<channel_size;++j) chem_rate_perturb_buf[j+buf_start]=0.0;
	}
	diff_buf = new double [domain_size];
	mob_buf = new double [domain_size];
	ddiff_buf = new double [domain_size];
	dmob_buf = new double [domain_size];
	total_gas_density = new double [domain_size];
	total_gas_pressure = new double [domain_size];
	series_thermal     = new double [domain_size];

	for (domain_int=0;domain_int<domain_size;++domain_int)
	{
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;

		diff_buf[domain_int] = 0.0;
		mob_buf[domain_int] = 0.0;
		ddiff_buf[domain_int] = 0.0;
		dmob_buf[domain_int] = 0.0;

		for (k=0;k<channel_size;++k)
		{
			if (channel[k].rate_type == 3)
			{
				chem_rate_buf[k+buf_start_channel]=channel[k].rate_constant;
			}
			else {chem_rate_buf[k+buf_start_channel]=0.0;}
			chem_rate_perturb_buf[k+buf_start_channel] = 0.0;
		}
	}
	/*
	for (domain_int=0;domain_int<domain_size;++domain_int){
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;
		for (k=0;k<channel_size;++k)
		{
		cout   << "domain: "<< domain_int <<" reaction:  " << k <<"  rate:   " << chem_rate_buf[k+buf_start_channel] << endl;
		}
		cout << endl;
	}
	for (domain_int=0;domain_int<domain_size;++domain_int){
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;
		for (k=0;k<channel_size;++k)
		{
		cout   << "domain: "<< domain_int <<" reaction:  " << k <<"  rate:   " << chem_rate_perturb_buf[k+buf_start_channel] << endl;
		}
		cout << endl;
	}
	exit(1);
	*/
	/**********************************************************************
	*  Design at 2009-04-21 for net source sink term                      *
	*								      *
	**********************************************************************/

	net_sourcesink = new apply_channel [channel_size];
	for (l=0;l<channel_size;++l)
	{
		net_sourcesink[l].sign = new int [species_size];
	}
	for (k=0;k<channel_size;++k)
	{
		for (i=0 ; i<species_size ;++i)
		{
			net_sourcesink[k].sign[i] = OPT_channel_sourcesink[k].coefficient_product[i] - OPT_channel_sourcesink[k].coefficient_reactant[i];
	        }
	}


	//exit(1);

///////////////////////////////////////////////////


	if(i_Rank == 0){
	  for (k=0;k<channel_size;++k){
		cout << "Reaction  " << k << ":  ";
		for (i=0;i<OPT_channel_sourcesink[k].sink_n;++i){
		  cout << OPT_channel_sourcesink[k].coefficient_reactant[OPT_channel_sourcesink[k].reactant_ID[i]]<< " " <<species.name[OPT_channel_sourcesink[k].reactant_ID[i]] << " + ";
		}
		  cout << " => " ;
		for (i=0;i<OPT_channel_sourcesink[k].source_n;++i){
		    cout << OPT_channel_sourcesink[k].coefficient_product[OPT_channel_sourcesink[k].product_ID[i]]<<" "<<species.name[OPT_channel_sourcesink[k].product_ID[i]] << " + ";
		}
		cout << endl;
	  }

	  for (k=0;k<channel_size;++k){
	    cout << "  Is it a momentum Xsfer:  " << OPT_channel_sourcesink[k].IS_momentum_Xsfer << endl;
	  }
    }

	//for (k=0;k<channel.size();++k)
	//{
	//	cerr << k << "\t" << channel[k].momentum_reactions << endl;
	//}

	//exit(1);

}







int   chemistry::ShowChemRate(const double *T_species)
{
	int    domain_int;
	int    k;
	int    buf_start;
	int    buf_start_channel;
	cout << " Below are shown every chemistry channel rate " << endl;
	for (domain_int=0;domain_int<domain_size;++domain_int)
	{
		buf_start = domain_int*species_size;
		buf_start_channel = domain_int*channel_size;
		for (k=0;k<channel_size;++k)
		{
		    cout   << "Grid_number : "<< domain_int <<" reaction:  " << k<< "  Te:  "<< *(T_species+species.e_id+buf_start) <<"  rate:   " << chem_rate_buf[k+buf_start_channel] << endl;
		}
		cout << endl;
	}
	return 0;
}



int   chemistry::ShowChemTable(const double *T_species)
{
	//Mark by KL      int    i,j,k;

	return 0;
}
double*  chemistry::ptr_single_ChannelSourceSink(int k)
{
	if (k >= channel.size() ) {cerr << "You input No of channel is larger than existed channels." << endl;exit(100);}
	return  (& channel[k].Channel_SourceSink[0]);
}



int chemistry::SpeciesID ( string  name )
{
	int  k;
	if (species.name.size() == 0) {cerr << "There is no species name " << endl;exit(1);}
	for (k=0;k<species.name.size();++k)
	{
		if (species.name[k] == name ) return k;
	}
	cerr << "You fill species name may be error in SpeciesID " << endl;
	cerr << "Program will leave now." << endl;
	exit(1);
	return 100000;
}


int chemistry::LightSize ( )
{
	return  species.light_name.size ( ) ;
}

vector<double*>   chemistry::ptr_source_sink()
{
	int                i;
	vector<double*>    ptr_s_s_temp;
	for (i=0;i<species.name.size() ;++i) ptr_s_s_temp.push_back(&species.sourcesink[i].global[0]);
	return  ptr_s_s_temp;
}

double*   chemistry::ptr_coll_frequency()
{
	return    &coll_frequency.global[0];
}

double*   chemistry::ptr_energy_loss()
{
	return    &energy_loss.global[0];
}

vector<double*>   chemistry::ptr_diffusion()
{
	int                i;
	vector<double*>    ptr_d_temp;
	for (i=0;i<species_size ;++i) ptr_d_temp.push_back(&species.diffusion[i].global[0]);
	return  ptr_d_temp;
}

vector<double*>   chemistry::ptr_mobility()
{
	int                i;
	vector<double*>    ptr_m_temp;
	for (i=0;i<species_size ;++i) ptr_m_temp.push_back(&species.mobility[i].global[0]);
	return  ptr_m_temp;
}

vector<double*>   chemistry::ptr_light_power()
{
	int                i;
	vector<double*>    ptr_m_temp;
	for (i=0;i<light_size ;++i) ptr_m_temp.push_back(&species.light_power[i].global[0]);
	return  ptr_m_temp;
}

vector<double*> chemistry::ptr_total_ChannelSourceSink()
{
	int                i;
	vector<double*>    ptr_m_temp;
	for (i=0;i<channel_size;++i) ptr_m_temp.push_back(&channel[i].Channel_SourceSink[0]);
	return  ptr_m_temp;
}



double chemistry::PS_PN(int x,int y,int d,double *D_species)
{
	int     i,j,k,l,reactant;
	double  result;
	double  temp;
	int     buf_start_channel;
	int     buf_start;
	buf_start_channel = d*channel_size;
	buf_start = d*species_size;
	result = 0.0;


	//for (kk=0;kk<species_size;++kk)
	//{

	//j=y*species_size+x;

	//cout << "partieal with "<< species.name[y] << endl;
	/*
	for (k=0;k<for_Formfunction[j].n;++k)
	{
	    cout <<"\t"<<for_Formfunction[j].knn[k].sign << " K_" <<for_Formfunction[j].knn[k].k ;
	    for (i=0;i<for_Formfunction[j].knn[k].nn_size;++i)
	    {
		if (for_Formfunction[j].knn[k].nn[i]==y)
		{
		if (for_Formfunction[j].knn[k].partial_c>0)
		cout <<" "<< species.name[for_Formfunction[j].knn[k].nn[i]]<<"^" <<for_Formfunction[j].knn[k].partial_c;
		}
		else
		{
		    cout <<" "<< species.name[for_Formfunction[j].knn[k].nn[i]] << "^" <<OPT_channel_sourcesink[for_Formfunction[j].knn[k].k].coefficient_reactant[for_Formfunction[j].knn[k].nn[i]] ;
		}
	    }
	    cout << "\t"<<for_Formfunction[j].knn[k].partial_c<<"^" << ;
	    cout << endl;

	}

	cout << endl;
	}
	*/


	j=y*species_size+x;
	result = 0.0;

	for (k=0;k<for_Formfunction[j].n;++k)
	{
		l = for_Formfunction[j].knn[k].k;
		temp = for_Formfunction[j].knn[k].sign * chem_rate_buf[l+buf_start_channel];
		//cout << " K  "<< for_Formfunction[j].knn[k].sign * chem_rate_buf[l+buf_start_channel] << endl;
		//cout << "  sign " <<  for_Formfunction[j].knn[k].sign << endl;
		//cout << "  K   " << chem_rate_buf[l+buf_start_channel] << endl;
		for (i=0;i<for_Formfunction[j].knn[k].nn_size;++i)
		{
		    reactant = for_Formfunction[j].knn[k].nn[i];
		    //cout << species.name[reactant]<< (*(D_species+reactant+buf_start))  << "^"<<OPT_channel_sourcesink[l].coefficient_reactant[reactant]<<endl;
		    if  ( reactant== y)
		    {
			    if (for_Formfunction[j].knn[k].partial_c > 0)
			    {
				if(for_Formfunction[j].knn[k].partial_c == 1)
				{
				    temp = temp * (*(D_species+reactant+buf_start))  ;
				    //cout << (*(D_species+reactant+buf_start)) <<endl;
				}
				else if (for_Formfunction[j].knn[k].partial_c > 1)
				{
				    temp = temp * pow((*(D_species+reactant+buf_start)),for_Formfunction[j].knn[k].partial_c);
				    //cout << (*(D_species+reactant+buf_start)) <<endl;
				}
				else {cout << "erroe at PS_PN" <<endl; exit(1);}
			    }

			    //cout << "1  "<< temp << endl;
		    }
		    else
		    {
			    if (OPT_channel_sourcesink[l].coefficient_reactant[reactant] == 1)
			    {
				    temp = temp * (*(D_species+reactant+buf_start));
				    //cout << (*(D_species+reactant+buf_start)) << endl;
			    }
			    else if (OPT_channel_sourcesink[l].coefficient_reactant[reactant] == 2)
			    {
				    temp = temp * pow ((*(D_species+reactant+buf_start)),2);
				    //cout << (*(D_species+reactant+buf_start)) << endl;
			    }
			    else if (OPT_channel_sourcesink[l].coefficient_reactant[reactant] == 3)
			    {
				    temp = temp * pow ((*(D_species+reactant+buf_start)),3);
				    //cout << (*(D_species+reactant+buf_start)) << endl;
			    }
			    else {cout << "erroe at PS_PN" <<endl; exit(1);}
			    //cout <<  "2  " << temp << endl;
		    }
		}
		//cout << temp << endl;
		//cout << "temp "<<temp << endl;
		//cout << "result "<< result << endl;

		result = result + temp;
	}
	//cout << "  result  01 " << result << endl;
	//exit(1);

	return result;
}




double    chemistry::PE_PN(int y,int d,double *D_species , double *T_species)
{
	int   m,l;
	double   result;
	double   temp;
	int      channell;
	int     buf_start_channel;
	int     buf_start;
	int     reactant;
	result =0.0;
	//cout << " partial with  "<<species.name[y] << endl;
	//cout << e_loss_for_Formfunction[y].n << endl;
	buf_start_channel = d*channel_size;
	buf_start         = d*species_size;
	for (l=0;l<e_loss_for_Formfunction[y].n;++l)
	{
		//cout << "test " << endl;
		//cout << e_loss_for_Formfunction[y].knn[l].sign ;
		//cout << " K_"<<e_loss_for_Formfunction[y].knn[l].k ;
		//cout << "test " << endl;
		channell = e_loss_for_Formfunction[y].knn[l].k ;
		//cout << " reaction  " << channell << endl;
		//cout << "  reaction type   " << e_channel_ID[channell] << endl;
		if (e_channel_ID[channell]==1)
		{
		    int    background;
		    temp = chem_rate_buf[channell+buf_start_channel];
		    //cout << " K  "<< temp << endl;
		    for (m=0;m<e_loss_for_Formfunction[y].knn[l].nn_size;++m)
		    {
			reactant = e_loss_for_Formfunction[y].knn[l].nn[m];
			if (reactant != species.e_id) {background = reactant;}
			if ( reactant!= y)
			{
			        //cout << species.name[reactant] << " " << *(D_species+reactant+buf_start) << endl;
				if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 1)
				{
				        temp = temp * (*(D_species+reactant+buf_start));
					//cout << (*(D_species+reactant+buf_start)) << endl;
				}
				else if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 2)
				{
					temp = temp * pow ((*(D_species+reactant+buf_start)),2);
					//cout << (*(D_species+reactant+buf_start)) << endl;
				}
				else if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 3)
				{
					temp = temp * pow ((*(D_species+reactant+buf_start)),3);
					 //cout << (*(D_species+reactant+buf_start)) << endl;
				}
				else {cout << "erroe at PE_PN" <<endl; exit(1);}
			    //cout << " "<< species.name[ e_loss_for_Formfunction[y].knn[l].nn[m] ]  <<"^" <<  OPT_channel_sourcesink[e_loss_for_Formfunction[y].knn[l].k].coefficient_reactant[e_loss_for_Formfunction[y].knn[l].nn[m]];
			}
			else
			{
				if ( e_loss_for_Formfunction[y].knn[l].partial_c> 0)
				{
				    //cout << species.name[reactant] << " " << *(D_species+reactant+buf_start) << endl;
				    if(e_loss_for_Formfunction[y].knn[l].partial_c == 1)
				    {
					temp = temp * (*(D_species+reactant+buf_start))  ;
					//cout << (*(D_species+reactant+buf_start)) <<endl;
				    }
				    else if (e_loss_for_Formfunction[y].knn[l].partial_c > 1)
				    {
					temp = temp * pow((*(D_species+reactant+buf_start)),e_loss_for_Formfunction[y].knn[l].partial_c);
					//cout << (*(D_species+reactant+buf_start)) <<endl;
				    }
				    else {cout << "erroe at PS_PN" <<endl; exit(1);}
				}
				//if (e_loss_for_Formfunction[y].knn[l].partial_c>0)
				//{
				//    cout << " "<< species.name[ e_loss_for_Formfunction[y].knn[l].nn[m] ]  <<"^" << e_loss_for_Formfunction[y].knn[l].partial_c ;
				//}
			}
		    }
				//cout << " Te  " << (*(T_species+species.e_id+buf_start)) << endl;
				//cout << "momentum  energy threshold  " << 3.0 * species.mass[species.e_id]/species.mass[background] * (*(T_species+species.e_id+buf_start))<<endl;
				temp = 3.0 * species.mass[species.e_id]/species.mass[background] *temp*(*(T_species+species.e_id+buf_start));
		}
		else if (e_channel_ID[channell]==2)
		{
		temp = e_loss_for_Formfunction[y].knn[l].sign * chem_rate_buf[channell+buf_start_channel];
		//cout << "K  "<<temp << endl;
		for (m=0;m<e_loss_for_Formfunction[y].knn[l].nn_size;++m)
		{
		    reactant = e_loss_for_Formfunction[y].knn[l].nn[m];

		    if ( reactant!= y)
		    {
			    //cout << species.name[reactant] << " " << *(D_species+reactant+buf_start) << endl;
			    if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 1)
			    {
				    temp = temp * (*(D_species+reactant+buf_start));
				    //cout << (*(D_species+reactant+buf_start)) << endl;
			    }
			    else if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 2)
			    {
				    temp = temp * pow ((*(D_species+reactant+buf_start)),2);
				    //cout << (*(D_species+reactant+buf_start)) << endl;
			    }
			    else if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 3)
			    {
				    temp = temp * pow ((*(D_species+reactant+buf_start)),3);
				    //cout << (*(D_species+reactant+buf_start)) << endl;
			    }
			    else {cout << "erroe at PE_PN" <<endl; exit(1);}
			//cout << " "<< species.name[ e_loss_for_Formfunction[y].knn[l].nn[m] ]  <<"^" <<  OPT_channel_sourcesink[e_loss_for_Formfunction[y].knn[l].k].coefficient_reactant[e_loss_for_Formfunction[y].knn[l].nn[m]];
		    }
		    else
		    {
			    if ( e_loss_for_Formfunction[y].knn[l].partial_c> 0)
			    {
				//cout << species.name[reactant] << " " << *(D_species+reactant+buf_start) << endl;
				if(e_loss_for_Formfunction[y].knn[l].partial_c == 1)
				{
				    temp = temp * (*(D_species+reactant+buf_start))  ;
				    //cout << (*(D_species+reactant+buf_start)) <<endl;
				}
				else if (e_loss_for_Formfunction[y].knn[l].partial_c > 1)
				{
				    temp = temp * pow((*(D_species+reactant+buf_start)),e_loss_for_Formfunction[y].knn[l].partial_c);
				    //cout << (*(D_species+reactant+buf_start)) <<endl;
				}
				else {cout << "erroe at PS_PN" <<endl; exit(1);}
			    }
			//if (e_loss_for_Formfunction[y].knn[l].partial_c>0)
			//{
			//    cout << " "<< species.name[ e_loss_for_Formfunction[y].knn[l].nn[m] ]  <<"^" << e_loss_for_Formfunction[y].knn[l].partial_c ;
			//}
		    }
		}
		//cout << "  threshold  e  " << e_loss_energy[channell] << endl;
		temp = temp *  e_loss_energy[channell];

		}

		//cout <<"temp  " <<temp << endl;
		result = result + temp;
		//cout << " result  " <<  result << endl;
		//cout << endl;
	}
	return  result;
}

int    chemistry::set_Te_perturbation(double dTe)
{
	perturb_Te=dTe;
	return 0;
}


double  chemistry::PS_PTe(int x,int d,double *D_species)
{
	int     l,m;
	double  temp;
	double  result;
	int     channell;
	int     reactant;
	int     buf_start;
	int     buf_start_channel;
	result = 0.0;
	buf_start_channel = d*channel_size;
	buf_start         = d*species_size;

		//cout << species.name[x] << endl;
	for (l = 0;l<S_dTe_ID[x].n;++l)
	{
		channell = S_dTe_ID[x].join[l];
		//cout << " reaction "<< channell;
		temp = S_dTe_ID[x].sign[l] * chem_rate_perturb_buf[channell+buf_start_channel];
		//cout << "  sign " <<  S_dTe_ID[x].sign[l] << " "<<chem_rate_perturb_buf[channell+buf_start_channel] << " "<< chem_rate_buf[channell+buf_start_channel] << endl;
		for (m=0;m< OPT_channel_sourcesink[channell].sink_n ; ++m)
		{
			reactant = OPT_channel_sourcesink[channell].reactant_ID[m];
			//cout << " "<<species.name[reactant] << " "<< (*(D_species+reactant+buf_start)) ;
			if ( OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 1 )
			{
				temp = temp*(*(D_species+reactant+buf_start));
			}
			else if ( OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 2 )
			{
				temp = temp*pow(*(D_species+reactant+buf_start),2);
			}
			else if ( OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 3 )
			{
				temp = temp*pow(*(D_species+reactant+buf_start),3);
			}
			else {cerr << "error at PS_PTe. exit" <<endl; exit(1);}
		}
		//cout << endl;
		//cout << "  temp  " << temp << endl;

		result = result + temp;
	}
	//exit (1);
	return result;

}

double   chemistry::PE_PTe(int d ,double*D_species,double*T_species)
{
	int  	channell;
	int  	l,m;
	int     reactant;
	int     background;
	double  temp;
	double  result;
	int     buf_start;
	int     buf_start_channel;
	result = 0.0;
	buf_start_channel = d*channel_size;
	buf_start         = d*species_size;

	for (l = 0;l<e_loss_dTe_ID.n;++l)
	{
		channell = e_loss_dTe_ID.join[l];
		if (e_channel_ID[channell]==1)
		{
			background = 100;
			temp = ((*(T_species+species.e_id+buf_start))*chem_rate_perturb_buf[channell+buf_start_channel] + chem_rate_buf[channell+buf_start_channel]);
			for (m=0;m< OPT_channel_sourcesink[channell].sink_n ; ++m)
			{
			    reactant = OPT_channel_sourcesink[channell].reactant_ID[m];
			    if (reactant != species.e_id) background = reactant;
			    if ( OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 1 )
			    {
				    temp = temp*(*(D_species+reactant+buf_start));
			    }
			    else if ( OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 2 )
			    {
				    temp = temp*pow(*(D_species+reactant+buf_start),2);
			    }
			    else if ( OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 3 )
			    {
				    temp = temp*pow(*(D_species+reactant+buf_start),3);
			    }
			    else {cerr << "error at PE_PTe. exit" <<endl; exit(1);}
			}
			if (background != 100)
			{
			    temp = temp * 3.0 * species.mass[species.e_id]/species.mass[background] ;
			}
			else {cerr << "error at PE_PTe. exit " << endl; exit(1); }
		}
		else if (e_channel_ID[channell]==2)
		{
			temp = chem_rate_perturb_buf[channell+buf_start_channel];
			for (m=0;m<OPT_channel_sourcesink[channell].sink_n;++m)
			{
				reactant = OPT_channel_sourcesink[channell].reactant_ID[m];
				if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 1)
				{
					temp = temp * (*(D_species+reactant+buf_start));
					//cout << (*(D_species+reactant+buf_start)) << endl;
				}
				else if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 2)
				{
					temp = temp * pow ((*(D_species+reactant+buf_start)),2);
					//cout << (*(D_species+reactant+buf_start)) << endl;
				}
				else if (OPT_channel_sourcesink[channell].coefficient_reactant[reactant] == 3)
				{
					temp = temp * pow ((*(D_species+reactant+buf_start)),3);
					//cout << (*(D_species+reactant+buf_start)) << endl;
				}
				else {cout << "erroe at PE_PN" <<endl; exit(1);}
			}
			temp = temp * e_loss_energy[channell];
		}
		result = result + temp;
		//cout << e_loss_dTe_ID.join[l] << endl;
	}

	return result;

}



double chemistry::PD_PTe( int d )
{
	return *(ddiff_buf+d);
}

double chemistry::PMu_PTe(int d )
{
	return *(dmob_buf+d);
}



double* chemistry::ptr_thermal(const double * T_species)
{
	int     buf_start;
	int     domain_int;
	int     background_id;
	int     location;
	double  T_neutral;
	background_id = species.name.size()-1;
	for (domain_int=0;domain_int<domain_size;++domain_int)
	{
		buf_start                  = domain_int*species_size;
		T_neutral                  = *(T_species + background_id + buf_start);
		if (T_neutral<=species.thermal_table[background_id].min_T)
		{
			series_thermal[domain_int] = species.thermal_table[background_id].value[0];
		}
		else if (T_neutral>=species.thermal_table[background_id].max_T)
		{
			series_thermal[domain_int] = species.thermal_table[background_id].value[ species.thermal_table[background_id].size_n-1 ];
		}
		else
		{
			location                   = (T_neutral-species.thermal_table[background_id].min_T)/species.thermal_table[background_id].dT;
			series_thermal[domain_int] = species.thermal_table[background_id].value[location];
		}
	}
	return series_thermal;
}






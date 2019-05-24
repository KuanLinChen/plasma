#include <vector>
#include <boost/fusion/tuple.hpp>
#include "datamanager.h"
#include "surface.h"
#include "sys_log.h"

using namespace std ;

DataManager::DataManager()
{
	map_iterator = 0 ;
}

DataManager::~DataManager()
{

}

void DataManager::release()
{
	int i ;
	for ( i = 0 ; i < map_iterator ; i++ )
	{
		if ( vec_dataexchanger_ptr[i]->comm != MPI_COMM_NULL ) MPI_Comm_free ( & ( vec_dataexchanger_ptr[i]->comm) ) ;
	}
}

void DataManager::clean()
{
	int i, j, counter ;
	int map_index;
	string sub_target, sub_source ;
	std::map<string, int>::iterator dit ;

	for ( i = vec_dataexchanger_ptr.size() - 1 ; i >= 0 ; i-- )
	{
		counter = 0;
		for ( std::set<string>::iterator sit = Domain::ID_bank.begin() ; sit != Domain::ID_bank.end(); sit++ )
		{
			sub_target = vec_dataexchanger_ptr[i]->DataExchanger_tag.substr( 5, 36 );
			if ( sub_target ==  *sit )
			{
				// domain target still exists
				counter ++ ;
			}
		}

		for ( std::set<string>::iterator sit = Domain::ID_bank.begin() ; sit != Domain::ID_bank.end(); sit++ )
		{
			sub_source = vec_dataexchanger_ptr[i]->DataExchanger_tag.substr( 42, 36 );
			if ( sub_source ==  *sit )
			{
				// domain source still exists
				counter ++ ;
			}
		}

		if ( counter != 2 )
		{
			
			dit = dataexchanger_index.find ( vec_dataexchanger_ptr[i]->DataExchanger_tag  ) ;
			
			for ( std::map<string, int>::iterator it=dataexchanger_index.begin(); it!=dataexchanger_index.end(); ++it )
			{
				if ( it->second > dit->second )
				{
					it->second = it->second -1 ;
				}
				
				if ( it == dit ) 
				{
					it->second = -1 ;
				}
			}
			
			dataexchanger_index.erase ( dit ) ;
			vec_dataexchanger_ptr.erase( vec_dataexchanger_ptr.begin() + i ) ;
			map_iterator -- ;
			
			//if ( mpi_id == 0 ) cout << endl;
		}
	}
}


int DataManager::if_register ( Domain *target, Domain *source, int type )
{
	string tag ;

	if ( type == CELL )
		tag = "CELL_" + target->ID + '_' + source->ID  ;
	else if ( type == FACE )
		tag = "FACE_" + target->ID + '_' + source->ID  ;

	// Checking the useless DataManagers and remove them.
	//Log( target->comm ).MPITagDump ( logLEVEL4 ) << "if_reg: " << tag ;
	clean();

	if (  !dataexchanger_index.empty()  && ( dataexchanger_index.count ( tag ) == 1 ) )
	{
		//cout << tag << endl;
		return dataexchanger_index[ tag ] ;
	} else
	{
		return -1 ;
	}
	
}


int DataManager::if_register ( Domain *target, Domain *source, string source_key , int type )
{
	string tag ;

	if ( type == FACE )
		tag = "FACE_" + target->ID + '_' + source->ID + '_' + source_key   ;
	else
		Log( target->comm ).MPITagDump ( logLEVEL4 ) << "Wrong type given in if_register() " << tag ;	

	// Checking the useless DataManagers and remove them.
	//Log( target->comm ).MPITagDump ( logLEVEL4 ) << "if_reg: " << tag ;
	clean();

	if (  !dataexchanger_index.empty()  && ( dataexchanger_index.count ( tag ) == 1 ) )
	{
		//cout << tag << endl;
		return dataexchanger_index[ tag ] ;
	} else
	{
		return -1 ;
	}
	
}

int DataManager::create_cell_exchanger ( Domain *target, Domain *source )
{
	int i, return_value ;

	string	tag ( "CELL_" + target->ID + '_' + source->ID ) ;

	if ( target->ID == source->ID )
	{
		Log().TagDump( logLEVEL4 ) << "Initializing GHOST exchanger (Ghost): " << tag  ;
		vec_dataexchanger_ptr.emplace_back(  new DataExchanger ( target ) ) ;
	} else if ( target->meshfile == source->meshfile )
	{
		Log().TagDump( logLEVEL4 ) << "Initializing homogeneous exchanger (H): " << tag  ;
		vec_dataexchanger_ptr.emplace_back(  new DataExchanger ( target, source ) ) ;
	} else
	{
		Log().TagDump( logLEVEL4 ) << "Initializing Non-homogeneous exchanger (NH): " << tag  ;
		vec_dataexchanger_ptr.emplace_back( new DataExchanger ( target, source ) ) ;
	}

	dataexchanger_index[ tag ] = map_iterator ;
	return_value = map_iterator ;
	map_iterator =  map_iterator + 1 ;
	return return_value  ;
}


void DataManager::cell_interpolation ( Domain *domain_target, Domain *domain_source, boost::shared_array<double> target, boost::shared_array<double> source ) 
{
	int id ;

	id = DataManager::if_register ( domain_target, domain_source, CELL ) ;

	if ( id < 0 ) 
	{
		id = DataManager::create_cell_exchanger ( domain_target, domain_source ) ;
	}
	vec_dataexchanger_ptr[id]->interpolationFunc( domain_target, domain_source, target, source, vec_dataexchanger_ptr[id].get() ) ;

	// Updating ghost cells
	if ( domain_target->ID != domain_source->ID )
	{
		id = DataManager::if_register ( domain_target, domain_target, CELL ) ;
		if ( id < 0 )
		{
			id = DataManager::create_cell_exchanger ( domain_target, domain_target ) ;
		} 
		vec_dataexchanger_ptr[id]->interpolationFunc( domain_target, domain_target, target, target, vec_dataexchanger_ptr[id].get() ); 
	}
}

int DataManager::create_face_exchanger ( Domain *target, Domain *source , string source_key, const vector<int> &source_face_id  )
{
	int i, return_value ;

	string	tag ( "FACE_" + target->ID + '_' + source->ID + '_' + source_key ) ;

	if ( target->ID == source->ID )
	{
		Log().TagDump( logLEVEL4 ) << "Initializing face exchanger (identical): " << tag  ;
		vec_dataexchanger_ptr.emplace_back(  new DataExchanger ( target, source_face_id, FACE, source_key  ) ) ;
	} else if ( target->meshfile == source->meshfile )
	{
		Log().TagDump( logLEVEL4 ) << "Initializing homogeneous exchanger (H): " << tag  ;
		vec_dataexchanger_ptr.emplace_back(  new DataExchanger ( target, source, source_face_id, FACE, source_key  ) ) ;
	} else
	{
		Log().TagDump( logLEVEL4 ) << "Non-homogeneous face exchanger is not ready yet. " << tag  ;
		//vec_dataexchanger_ptr.emplace_back( new DataExchanger ( target, source ) ) ;
	}

	dataexchanger_index[ tag ] = map_iterator ;
	return_value = map_iterator ;
	map_iterator =  map_iterator + 1 ;
	return return_value  ;
}

void DataManager::face_interpolation ( Domain *domain_target, Domain *domain_source, string source_key, boost::shared_array<double> target, boost::shared_array<double> source, const vector<int> &source_face_id  ) 
{
	int id ;

	id = DataManager::if_register ( domain_target, domain_source, CELL ) ;

	if ( id < 0 ) 
	{
		id = DataManager::create_face_exchanger ( domain_target, domain_source , source_key , source_face_id ) ;
	}
	vec_dataexchanger_ptr[id]->interpolationFunc( domain_target, domain_source, target, source, vec_dataexchanger_ptr[id].get() ) ;

	// Updating ghost cells
	if ( domain_target->ID != domain_source->ID )
	{
		id = DataManager::if_register ( domain_target, domain_target, CELL ) ;
		if ( id < 0 )
		{
			id = DataManager::create_face_exchanger ( domain_target, domain_target, source_key , source_face_id  ) ;
		} 
		vec_dataexchanger_ptr[id]->interpolationFunc( domain_target, domain_target, target, target, vec_dataexchanger_ptr[id].get() ); 
	}	
	
}

void DataManager::face_interpolation ( Domain *domain_target, Domain *domain_source, boost::shared_array<double> target, boost::shared_array<double> source ) 
{

	Log().TagDump( logLEVEL4 ) << "FULL face interpolation - Not plane to do it yet " ;
}

void DataManager::remove_relation ( Domain * target, Domain * source, int  i) 
{
	int id ;
	map<string, int>::iterator it;
	string	tag ( target->meshfile + source->meshfile )  ;
	if ( i == CELL ) tag = "CELL" + tag ;
	it = dataexchanger_index.find ( tag ) ;
	id = if_register ( target,  source, CELL ) ;


	vec_dataexchanger_ptr[id].reset() ;
	//vec_dataexchanger_ptr[id] = NULL ;

	dataexchanger_index.erase( it ) ;

}

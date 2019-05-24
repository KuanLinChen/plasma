#include <cstdio>
#include "sys_log.h"
#include <mpi.h>
#include <cstring>

using namespace std ;

LogLevel Log::reportingLevel = logLEVEL4 ;
ofstream* Log::_outfile = nullptr;

//INITIALIZE_EASYLOGGINGPP

string NowTime()
{
	char buffer[11];
	time_t t;

	time ( &t );
	tm r = {0};

	strftime ( buffer, sizeof(buffer), "%X", localtime_r( &t, &r ) );
	struct timeval tv;
	gettimeofday(&tv, 0);

	char result[100] = {0};
	sprintf (result, "%s.%03ld", buffer, (long) tv.tv_usec / 1000 );

	return result;
}

Log::Log()
{
	MPI_flag = false ;
	comm = MPI_COMM_WORLD ;
	int i ;
	MPI_Finalized( &i ) ;

	if ( i )
	{
		MPI_SIZE = 0 ;
		MPI_ID = 0 ;
	} else
	{
		MPI_Comm_size( comm, &MPI_SIZE ) ;
		MPI_Comm_rank ( comm, &MPI_ID ) ;
	}
}

Log::Log( MPI_Comm _comm )
{
	MPI_flag = true ;
	comm = _comm ;
	MPI_Comm_size( comm, &MPI_SIZE ) ;
	MPI_Comm_rank ( comm, &MPI_ID ) ;
}

Log::Log(ofstream *outfile)
{
	_outfile = outfile;

	MPI_flag = false ;
	comm = MPI_COMM_WORLD ;
	int i ;
	MPI_Finalized( &i ) ;

	if ( i )
	{
		MPI_SIZE = 0 ;
		MPI_ID = 0 ;
	} else
	{
		MPI_Comm_size( comm, &MPI_SIZE ) ;
		MPI_Comm_rank ( comm, &MPI_ID ) ;
	};
};

Log::Log(ofstream *outfile, MPI_Comm _comm)
{
	_outfile = outfile;

	MPI_flag = true ;
	comm = _comm ;
	MPI_Comm_size( comm, &MPI_SIZE ) ;
	MPI_Comm_rank ( comm, &MPI_ID ) ;

};



Log::~Log()
{
	int i, j ;
    char *recv_s  ;

	if ( comm == MPI_COMM_NULL ) return;

	ostream *gbl_out;

	if(_outfile != nullptr)
	{
		gbl_out = _outfile;
	}
	else
	{
		gbl_out = &cout;
	}


	if ( recent_level <= reportingLevel  )
	{
		if ( MPI_flag )
		{
			if ( MPI_ID == 0 )
			{
				os << std::endl;
				//fprintf ( stdout, "%s", os.str().c_str() );
				*gbl_out << os.str();

				for ( i = 1 ; i < MPI_SIZE ; i++ )
				{

					MPI_Recv ( &j, 1, MPI_INT, i, 0, comm, MPI_STATUS_IGNORE ) ;
					recv_s = new char [j] ;
					MPI_Recv ( recv_s, j, MPI_CHAR, i, 0, comm, MPI_STATUS_IGNORE ) ;
					//fprintf ( stdout, "%s\n", recv_s );
					*gbl_out << recv_s << endl;
					delete [] recv_s ;
				}

				//fflush ( stdout );
			} else
			{
				int s_size = os.str().size() + 1 ;

				char *s = new char[ s_size ];
				strcpy( s, os.str().c_str() );

				MPI_Send ( &s_size, 1, MPI_INT, 0, 0, comm  ) ;
				MPI_Send ( s, s_size, MPI_CHAR, 0, 0, comm ) ;

				delete [] s;
			}
		} else
		{
			if ( MPI_ID == 0 )
			{
				os << std::endl;
				//fprintf ( stdout, "%s", os.str().c_str() );
				//fflush ( stdout );
				*gbl_out << os.str();
			}
		}
	}
}

ostringstream& Log::Dump ( LogLevel level )
{
	recent_level = level ;

	if ( MPI_flag )
	{
		if ( MPI_ID == 0)
			os << " " << std::string( level > logLEVEL0 ? level - logLEVEL0 : 0, '\t');
		else
			os.str( "" ) ;
	} else
	{
		os << " " << std::string( level > logLEVEL0 ? level - logLEVEL0 : 0, '\t');
	}

	return os;
}

ostringstream& Log::TagDump ( LogLevel level )
{
	recent_level = level ;

	if ( MPI_flag )
	{
		if ( MPI_ID == 0 )
		{
			os << " " << NowTime();
			os << "[" << ToString(level) << "]: ";
			os << std::string( level > logLEVEL0 ? level - logLEVEL0 : 0, '\t');
		} else
		{
			os.str( "" ) ;
		}
	} else
	{
		os << " " << NowTime();
		os << "[" << ToString(level) << "]: ";
		os << std::string( level > logLEVEL0 ? level - logLEVEL0 : 0, '\t');
	}


	return os;
}

ostringstream& Log::MPIDump ( LogLevel level)
{
	if ( comm == MPI_COMM_NULL ) return os ;

	recent_level = level ;
	os << "[" << MPI_ID << "]:" << " " << std::string( level > logLEVEL0 ? level - logLEVEL0 : 0, '\t');

	return os;
}

ostringstream& Log::MPITagDump ( LogLevel level )
{
	if ( comm == MPI_COMM_NULL ) return os ;

	recent_level = level ;
	os << " " << NowTime();
	os << "[" << ToString(level) << "][" << MPI_ID << "]: "  ;
	os << std::string( level > logLEVEL0 ? level - logLEVEL0 : 0, '\t');

	return os;
}

LogLevel& Log::ReportingLevel()
{
	return reportingLevel;
}

void Log::ReportingLevel( LogLevel L )
{
	reportingLevel = L ;
}

string Log::ToString( LogLevel level )
{
	static const char* const buffer[] = {"E", "W", "I", "L", "L", "L", "L", "L"};
	return buffer[level];
}


LogLevel Log::FromString( const string& level )
{
	if (level == "LEVEL4")
		return logLEVEL4;
	if (level == "LEVEL3")
		return logLEVEL3;
	if (level == "LEVEL2")
		return logLEVEL2;
	if (level == "LEVEL1")
		return logLEVEL1;
	if (level == "LEVEL0" || level == "LEVEL" )
		return logLEVEL0;
	if (level == "INFO")
		return logINFO;
	if (level == "WARNING")
		return logWARNING;
	if (level == "ERROR")
		return logERROR;

	Log().TagDump ( logWARNING ) << "Unknown logging level '" << level << "'. Using INFO level as default.";

	return logINFO;
}


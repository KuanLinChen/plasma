#include <sstream>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include <mpi.h>

using namespace std;

#if !defined(__SYS_LOG_H)
#define __SYS_LOG_H

string NowTime() ;

enum LogLevel
{
	logERROR,
	logWARNING,
	logINFO,
	logLEVEL0,
	logLEVEL1,
	logLEVEL2,
	logLEVEL3,
	logLEVEL4
};

class Log
{
	public:
		Log();
		Log( MPI_Comm );
		Log( ofstream *outfile);
		Log( ofstream *outfile, MPI_Comm _comm);
		virtual ~Log();
		ostringstream & Dump( LogLevel level = logINFO );
		ostringstream & TagDump( LogLevel level = logINFO );

		ostringstream & MPIDump( LogLevel level = logINFO );
		ostringstream & MPITagDump( LogLevel level = logINFO );

		bool MPI_flag ;
		int MPI_ID, MPI_SIZE ;
		MPI_Comm comm;

		LogLevel recent_level ;

		static LogLevel& ReportingLevel();
		static void ReportingLevel( LogLevel );
		static LogLevel reportingLevel  ;
		static string ToString( LogLevel level);
		static LogLevel FromString(const std::string& level);
		static ofstream *_outfile;

	protected:
		std::ostringstream os;

	private:
		Log (const Log&);
		Log& operator =(const Log&);
};


#endif

#include <cstdarg>
#include <iostream>

using namespace std;

// this function will take the number of values to average
// followed by all of the numbers to average
void DEBUG_STD ( int num, ... )
{
	int 	i; 
	va_list	arguments; 
	double 	sum = 0;
	
	va_start ( arguments, num ); 

	for ( i = 0; i < num; i++ )
    	sum += va_arg ( arguments, double ); 
  
  	va_end ( arguments );

	//return 0;   
}


void DEBUG_printf( const char *fmt, ... )
{
	va_list args;
	char *ptr ;
	
	va_start (args, fmt);

	ptr = const_cast<char *> ( fmt ) ;

	while( *ptr )
	{
		switch( *ptr )
		{
			case '\\':
			if( *++ptr )
				ptr++;
			continue;

			case '%':
			switch( *++ptr )
			{
				case NULL:
				continue;

				case 's':
					printf("%s", va_arg ( args, char *) );
				break;

				case 'd':
					printf("%d", va_arg( args, int) );
				break;

				case 'f':
					printf("%f", va_arg( args, double) );
				break;
			}
			ptr++;

			default:
				putchar( *ptr++ );
		}
	}
	

	va_end (args);
}
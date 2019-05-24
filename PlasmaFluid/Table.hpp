#pragma once
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
class CTable{

	private:
		double *x, *y ;
	public:
	CTable() {
		ReadTable = false ;
	} ;
	CTable( string FileName ) {
		Name = FileName ;
		DataNum = FileCounter() ;
		printf("Creat Table Container: %s, Data Number: %d\n",Name.c_str(), DataNum ) ;
		Alloc() ;
		Read() ;
	};
	string Name ;
	int DataNum ;
	bool ReadTable ;
	
	void Init( string FileName ){

		Name = FileName ;
		DataNum = FileCounter() ;
		printf("Creat Table Container: %s, Data Number: %d\n", Name.c_str(), DataNum ) ;
		Alloc();
		Read() ;
	};


	void Alloc(){
		x =  new double [ DataNum ] ;
		y =  new double [ DataNum ] ;
	}
	
	int FileCounter(){
        int counter = 0 ;
        fstream file ;
        string  line ;

        file.open( Name.c_str(), ios::in ) ;
        if (!file){printf(" File Counter Error: File Name => %s\n",Name.c_str() ) ; exit(1);}

        while( getline( file, line ) ){
			line = analysor( line ) ;
			if( line.size() > 0 and line == "EndDescription" ) break ;
		}
		while( getline( file, line ) ){
			if( line.size() > 0 ) counter++ ;
		}
		file.close() ;
		file.clear() ;
		if( counter == 0 ){printf(" File No data : File Name => %s\n",Name.c_str() ) ; exit(1);}
        return counter ;
	};

	string analysor( string s ){
        s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );
        return s;
	}
	void Read(){
        string  line ;
        fstream file ;
        file.open( Name.c_str(), ios::in ) ;
        if (!file){printf(" File Counter Error: File Name => %s ",Name.c_str() ) ; exit(1);}

        while( getline( file, line ) ){
			line = analysor( line ) ;
			if( line.size() > 0 and line == "EndDescription" ) break ;
		}
		for(int i=0 ; i < DataNum ; i++ ){
			file>> x[ i ]>> y[ i ] ;
			cout<<x[ i ]<<"\t"<<y[ i ]<<endl;
		}
		file.close() ;
		file.clear() ;
		ReadTable = true ;
	};
	double GetValue( double T )
	{
		double value = 0.0 ;
		if( T <= x[0] ){
			value = y[0] ;
		}else if( T >= x[DataNum-1] ){
			value = y[DataNum-1] ;
		}else{
			for(int i = 0 ; i < DataNum ; i++ ){
				if( T > x[i] and T <= x[i+1]){	
					value = y[i] + ( T-x[i] )*( y[i+1] - y[i] )/(  x[i+1]-x[i] )  ;
					break;
				}
			}
		}
		return value ;
	}
	double GetValueLog( double T )
	{
		double value = 0.0 ;
		if( T <= x[0] ){
			value = y[0] ;
		}else if( T >= x[DataNum-1] ){
			value = y[DataNum-1] ;
		}else{
			for(int i = 0 ; i < DataNum ; i++ ){
				if( T > x[i] and T <= x[i+1]){	
					value = pow(10.0, log10(y[i]) + ( log10(T)-log10(x[i]) )*( log10(y[i+1]) - log10(y[i]) )/(  log10(x[i+1])-log10(x[i]) )  );
					break;
				}
			}
		}
		return value ;
	}
};

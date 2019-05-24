#include <boost/shared_array.hpp>
#include "datapool.h"
#include "datamanager.h"
#include "domain.h"

#if !defined(__SCALAR_H)
#define __SCALAR_H

class Scalar
{
	public:
		Scalar	() ;
		~Scalar	() ;
		Scalar	( int ) ;
		Scalar	( boost::shared_ptr<Domain>, int ) ;
		Scalar	( boost::shared_ptr<Domain>, int, string ) ;
	
		string ID ;
		static set<string> ID_bank ;


		bool	share_memory_flag ;
		int		type;
		int		data_number ;
		int		local_data_number  ;
		boost::shared_ptr<Domain>	domain_ptr ;
		boost::shared_array<double> data ;
		string	name ;

		int 	*ref_local_id,  *ref_global_id ;
		bool	ref_flag ;

		//void	copy ( const Scalar &  ) ;
		//void	copy_value ( const Scalar &  ) ;
		void	share_memory_from ( const Scalar & ) ;
		void	initial	( boost::shared_ptr<Domain>, int ) ;
		void	initial	( boost::shared_ptr<Domain>, int, string ) ;
		//void	initial_as_reference ( Scalar * , string  ) ;
		void	zero () ;
		void	remove_relation (  const Scalar & ) ;

		vector<int> collected_faces ;		
		void    collecting_face ( string ) ;
		void    collecting_face ( ) ;
		void    collecting_face ( double * ) ;
		void    reset_collecting_face () ;


		operator int () { return local_data_number ; } 
		operator double * () { return data.get() ; } 
		operator Domain * () { return domain_ptr.get() ; } 
		operator boost::shared_array<double> () { return data; } 
		
		double & operator[] ( int i ) { return data[i] ;}  // For Scalar[i] = value ; or value = Scalar[i] ;
		const double & operator[] ( int i ) const { return data[i] ;}

		Scalar & operator= ( const Scalar & ) ; // This should take care the domain interpolation. 
		Scalar & operator= ( const double*  ) ; // This should take care the domain interpolation. 

		friend bool operator< ( int &, Scalar &) ;
		friend bool operator< ( Scalar &, int &) ;
		friend bool operator> ( int &, Scalar &) ;
		friend bool operator> ( Scalar &, int &) ;
		friend bool operator<= ( int &, Scalar &) ;
		friend bool operator<= ( Scalar &, int &) ;
		friend bool operator>= ( int &, Scalar &) ;
		friend bool operator>= ( Scalar &, int &) ;
		friend bool operator== ( int &, Scalar &) ;
		friend bool operator== ( Scalar &, int &) ;
		friend bool operator!= ( int &, Scalar &) ;
		friend bool operator!= ( Scalar &, int &) ;

		static DataManager datamanager ;
	private:

		// For face scalar
		bool flag_collecting_face ;
};



#endif

#include <boost/shared_array.hpp>
#include <boost/multi_array.hpp>

using namespace std ;

#if !defined(__ELEMENT_H)
#define __ELEMENT_H

class Node ;
class Face ;
class Cell ;
//class Domain;


class Node
{
	public:
		Node() ;
		~Node() ;
		int		id, local_id, type ;
		int		face_number ;
		int		cell_number ;
		double	x, y, z ;

		boost::shared_array< Face* > face ;
		boost::shared_array< Cell* > cell ;
		boost::shared_array< double > alpha ;

		void set_id ( int g, int l ) ; //( global, local )
		void set_position ( double p1, double p2, double p3 ) ;
		void set_face_number ( int size ) ;
		void set_cell_number ( int size ) ;

		// r[ 0 ] = x, r[ 1 ] = y, r[ 2 ] = z
		double	r[ 3 ] ;
} ;

class Face
{
	public:
		Face() ;
		~Face() ;
		int		mpi_id, id, local_id, type ;
		string	Typename ;
		double	x, y, z, dx, dy, dz ;

		boost::shared_array < double > 	cell_sign ;
		boost::shared_array < double > 	alpha ;
		boost::shared_array < double > 	dr ;
		boost::shared_array < double > 	dr_x ;
		boost::shared_array < double > 	dr_y ;

		/* dl is the distance between two cell centers cross the face.*/
		/* dA is the face area */
		double	dl, dA, angle;
		//double	dl, dA, angle, Ax, Ay, Az, Ad, nx, ny, nz, nd ;

		/* (cross_x, cross_y) is the coordinate of cross point of the face and line bwtween cell centers */
		double	pdl[ 2 ], cross_x, cross_y, cross_z, dl0_cos, dl1_cos ;

		int	node_number ; /* How many nodes connected. */
		int	cell_number ; /* How many cells connected. */
		boost::shared_array < Node* > node ;
		boost::shared_array < Cell* > cell ;

		void set_id ( int g, int l ) ; //( global, local )
		void set_type ( int t, string tn ) ;
		void set_position ( double p1, double p2, double p3 ) ;
		void set_node_number ( int size ) ;
		void set_cell_number ( int size ) ;
		// r[ 3 ]		=	( x, y, z )
		// A[ 3 ]		=	( Ax, Ay, Az )
		// nA[ 3 ]		=	( nx, ny, nz )
		// xi[ 3 ]		=	( x_xi, y_xi, z_xi )
		// cross_r[ 3 ]		=	( cross_x, cross_y, cross_z )
		// r_cf[ 0 ][ 3 ]	=	( r_f(x) - r_c0(x), r_f(y) - r_c0(y), r_f(z) - r_c0(z) ); r_cf[ 1 ][ 3 ] = ( r_f(x) - r_c1(x), r_f(y) - r_c1(y), r_f(z) - r_c1(z) )
		// dr_cf[ 0 ]		=	dist( r_f, r_c0 ), ...
		// dl_cos[ 0 ]		=	dl0_cos, dl_cos[ 1 ]	=	dl1_cos

		double	r[ 3 ], A[ 3 ], nA[ 3 ], cross_r[ 3 ], dl_cos[ 2 ], r_fof[ 3 ], nAd ;
		double	r_cf[ 2 ][ 3 ], dr_cf[ 2 ], cf_alpha[ 2 ], cc_alpha[ 2 ], diffusion_coefficient_1st, diffusion_coefficient_2nd[ 2 ][ 3 ], r_cc[ 3 ], nr_cc[ 3 ] ;

		//double  _xi[ 3 ], _eta[ 3 ], _zeta[ 3 ] ;
		//double  x_xi, y_xi, z_xi, x_eta, y_eta, xi_coefficient, eta_coefficient ;
		//double	A_dot_xi ;

	private:
} ;


/*! \brief class Cell

@param id		The global ID of cell
@param local_id		The local ID of cell
@param volume		The cell volume
@param face_sign	The direction of the face in one cell, which might be \f$+1\f$ or \f$-1\f$. The sign is opposite to the neighbor cell cross this face.
@param face_index	The index is for the convenience to give an index to neighoor cell of a face. If the face_sign=1, the face_index will be 0, otherwise will be 1.
*/

class Cell
{
	public:
		Cell() ;
		~Cell() ;

		int		mpi_id, id, local_id, mesh_id, form, type ;
		double	x, y, z, volume ;
		int		node_number, face_number, cell_number ;
		string	Typename ;

		boost::shared_array < double > face_sign ;
		boost::shared_array < double > dr ;
		boost::shared_array < int > face_index_c0 ;
		boost::shared_array < int > face_index_c1 ;
		//boost::shared_array < double > Jacobian ;
		//boost::shared_array < double > x_xi, y_xi, z_xi ;
		//boost::shared_array < double > x_eta, y_eta, z_eta ;
		//boost::shared_array < double > x_zeta, y_zeta, z_zeta ;
		boost::shared_array < double > Ax, Ay, Az ;
		//boost::shared_array < double > A_xi, A_eta, A_zeta ;
		boost::shared_array < double > nx, ny, nz ;
		//boost::shared_array < double > xi_coefficient ;
		boost::shared_array < Node* > node ;
		boost::shared_array < Face* > face ;
		boost::shared_array < Cell* > cell ;

		// New version
		double	r[ 3 ] ;

		//boost::multi_array < double , 2 > _xi ;
		//boost::multi_array < double , 2 > _eta ;
		//boost::multi_array < double , 2 > _zeta ;
		boost::shared_array < vector<double> > A  ;
		boost::shared_array < vector<double> > nA ;
		//boost::shared_array < double > A_dot_xi;

		void set_id ( int g, int l, int m ) ;  //( global, local, mesh )
		void set_type ( int t , string tn ) ;
		void set_position ( double p1, double p2, double p3 ) ;
		void set_volume ( double v ) ;

		void set_face_number ( int size ) ;
		void set_node_number ( int size ) ;
		void set_cell_number ( int size ) ;
} ;

double dist( double x1, double y1, double z1, double x2, double y2, double z2 ) ;
void check_neighbor_data( vector<int> *meshdata, int *table_mesh_local, vector<int> *neighbordata ) ;

#endif

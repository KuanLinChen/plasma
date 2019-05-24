#include "domain.h"

void Domain::CalculateFacedAdL( )
{
	int		i, j, k ;
	double	dx, dy, dz ;
	double  r, dr, r_cf_dot_A ;
	double	value[ 2 ][ 3 ] ;

	if ( cylindrical_y ==  1 )
	{
		for ( int i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			r	=	0.5 * ( face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].node[ 0 ]->r[ 0 ] ) ;
			dr	=	face[ i ].node[ 1 ]->r[ 0 ] - face[ i ].node[ 0 ]->r[ 0 ] ;
			dz	=	face[ i ].node[ 1 ]->r[ 1 ] - face[ i ].node[ 0 ]->r[ 1 ] ;

			if ( fabs ( r ) < 1.e-25 )
			{
				// boundary face at center line
				r					=	0 ;

				face[ i ].A[ 0 ]	=	0. ;
				face[ i ].A[ 1 ]	=	0. ;
				face[ i ].A[ 2 ]	=	0. ;
				face[ i ].dA		=	0. ;

				face[ i ].nA[ 0 ]	=	dz / fabs( dz ) ;
				face[ i ].nA[ 1 ]	=	0. ;
				face[ i ].nA[ 2 ]	=	0. ;
			} else
			{
				face[ i ].A[ 0 ]	=	dz * 2. * M_PI * r ;
				face[ i ].A[ 1 ]	=	- dr * 2. * M_PI * r ;
				face[ i ].A[ 2 ]	=	0. ;
				face[ i ].dA		=	sqrt ( dr * dr + dz * dz ) * 2. * M_PI * r ;

				for ( j = 0 ; j < 2 ; j++ )
				{
					face[ i ].nA[ j ]	=	face[ i ].A[ j ] / face[ i ].dA ;
					if ( abs( face[ i ].nA[ j ] ) < 1e-10 )
						face[ i ].nA[ j ]	=	0. ;
				}
				face[ i ].nA[ 2 ]	=	0. ;
			}

			face[ i ].nAd    =   face[ i ].nA[ 0 ] * face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].nA[ 1 ] * face[ i ].node[ 1 ]->r[ 1 ] ;
			//face[ i ].Ad    =   face[ i ].A[ 0 ] * face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].A[ 1 ] * face[ i ].node[ 1 ]->r[ 1 ] ;
			//face[ i ].nd    =   face[ i ].Ad / face[ i ].dA ;

			face[ i ].dl	=	0. ;
			for ( j = 0 ; j < 2 ; j++ )
			{
				if ( face[ i ].cell_number == 1 )
					face[ i ].r_cc[ j ]	=	face[ i ].r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;
				else
					face[ i ].r_cc[ j ]	=	face[ i ].cell[ 1 ]->r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;

				face[ i ].dl	+=	face[ i ].r_cc[ j ] * face[ i ].r_cc[ j ] ;
			}
			face[ i ].dl	=	sqrt ( face[ i ].dl ) ;

			for ( j = 0 ; j < 2 ; j++ )
				face[ i ].nr_cc[ j ]	=	face[ i ].r_cc[ j ] / face[ i ].dl ;

			for ( j = 0 ; j < face[ i ].cell_number ; j++ )
			{
				face[ i ].dr_cf[ j ]	=	0. ;
				for ( k = 0 ; k < 2 ; k++ )
				{
					face[ i ].r_cf[ j ][ k ]	=	face[ i ].r[ k ] - face[ i ].cell[ j ]->r[ k ] ;
					face[ i ].dr_cf[ j ]		+=	face[ i ].r_cf[ j ][ k ] * face[ i ].r_cf[ j ][ k ] ;
				}
				face[ i ].dr_cf[ j ]	=	sqrt( face[ i ].dr_cf[ j ] ) ;
			}
		}
	} else if ( cylindrical_x ==  1 )
	{
		for ( int i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			r	=	0.5 * ( face[ i ].node[ 1 ]->r[ 1 ] + face[ i ].node[ 0 ]->r[ 1 ] ) ;
			dr	=	face[ i ].node[ 1 ]->r[ 1 ] - face[ i ].node[ 0 ]->r[ 1 ] ;
			dz	=	face[ i ].node[ 1 ]->r[ 0 ] - face[ i ].node[ 0 ]->r[ 0 ] ;

			if ( fabs ( r ) < 1.e-25 )
			{
				// boundary face at center line
				r					=	0 ;

				face[ i ].A[ 0 ]	=	0. ;
				face[ i ].A[ 1 ]	=	0. ;
				face[ i ].A[ 2 ]	=	0. ;
				face[ i ].dA		=	0. ;

				face[ i ].nA[ 0 ]	=	0. ;
				face[ i ].nA[ 1 ]	=	dz / fabs( dz ) ;
				face[ i ].nA[ 2 ]	=	0. ;
			} else
			{
				face[ i ].A[ 0 ]	=	dr * 2. * M_PI * r ;
				face[ i ].A[ 1 ]	=	- dz * 2. * M_PI * r ;
				face[ i ].A[ 2 ]	=	0. ;
				face[ i ].dA		=	sqrt ( dr * dr + dz * dz ) * 2. * M_PI * r ;

				for ( j = 0 ; j < 2 ; j++ )
				{
					face[ i ].nA[ j ]	=	face[ i ].A[ j ] / face[ i ].dA ;
					if ( abs( face[ i ].nA[ j ] ) < 1e-10 )
						face[ i ].nA[ j ]	=	0. ;
				}
				face[ i ].nA[ 2 ]	=	0. ;
			}

			face[ i ].nAd    =   face[ i ].nA[ 0 ] * face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].nA[ 1 ] * face[ i ].node[ 1 ]->r[ 1 ] ;
			//face[ i ].Ad    =   face[ i ].A[ 0 ] * face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].A[ 1 ] * face[ i ].node[ 1 ]->r[ 1 ] ;
			//face[ i ].nd    =   face[ i ].Ad / face[ i ].dA ;

			face[ i ].dl	=	0. ;
			for ( j = 0 ; j < 2 ; j++ )
			{
				if ( face[ i ].cell_number == 1 )
					face[ i ].r_cc[ j ]	=	face[ i ].r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;
				else
					face[ i ].r_cc[ j ]	=	face[ i ].cell[ 1 ]->r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;

				face[ i ].dl	+=	face[ i ].r_cc[ j ] * face[ i ].r_cc[ j ] ;
			}
			face[ i ].dl	=	sqrt ( face[ i ].dl ) ;

			for ( j = 0 ; j < 2 ; j++ )
				face[ i ].nr_cc[ j ]	=	face[ i ].r_cc[ j ] / face[ i ].dl ;

			for ( j = 0 ; j < face[ i ].cell_number ; j++ )
			{
				face[ i ].dr_cf[ j ]	=	0. ;
				for ( k = 0 ; k < 2 ; k++ )
				{
					face[ i ].r_cf[ j ][ k ]	=	face[ i ].r[ k ] - face[ i ].cell[ j ]->r[ k ] ;
					face[ i ].dr_cf[ j ]		+=	face[ i ].r_cf[ j ][ k ] * face[ i ].r_cf[ j ][ k ] ;
				}
				face[ i ].dr_cf[ j ]	=	sqrt( face[ i ].dr_cf[ j ] ) ;
			}
		}
	} else if ( dimension == 2 )
	{
		for ( int i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			dx	=	face[ i ].node[ 1 ]->r[ 0 ] - face[ i ].node[ 0 ]->r[ 0 ] ;
			dy	=	face[ i ].node[ 1 ]->r[ 1 ] - face[ i ].node[ 0 ]->r[ 1 ] ;

			//face[ i ].Ax	=	dy ;
			//face[ i ].Ay	=	- dx ;
			//face[ i ].Az	=	0. ;
			//face[ i ].Ad    =	face[ i ].Ax * face[ i ].node[ 1 ]->x + face[ i ].Ay * face[ i ].node[ 1 ]->y ;
			face[ i ].dA	=	sqrt( dx * dx + dy * dy ) ;
			//face[ i ].nx	=	face[ i ].Ax / face[ i ].dA ;
			//if ( abs( face[ i ].nx ) < 1e-10 ) face[ i ].nx	=	0 ;
			//face[ i ].ny	=	face[ i ].Ay / face[ i ].dA ;
			//if ( abs( face[ i ].ny ) < 1e-10 ) face[ i ].ny	=	0 ;
			//face[ i ].nz	=	0. ;
			//face[ i ].nd   =   face[ i ].Ad / face[ i ].dA ;

			if ( face[ i ].cell_number > 1 )
			{
				dx	=	face[ i ].cell[ 0 ]->r[ 0 ] - face[ i ].cell[ 1 ]->r[ 0 ] ;
				dy	=	face[ i ].cell[ 0 ]->r[ 1 ] - face[ i ].cell[ 1 ]->r[ 1 ] ;

				face[ i ].dl		=	sqrt( dx * dx + dy * dy ) ;

				face[ i ].dr_x[ 0 ]	=	face[ i ].r[ 0 ] - face[ i ].cell[ 0 ]->r[ 0 ] ;
				face[ i ].dr_y[ 0 ]	=	face[ i ].r[ 1 ] - face[ i ].cell[ 0 ]->r[ 1 ] ;
				face[ i ].dr[ 0 ]	=	sqrt( face[ i ].dr_x[ 0 ] * face[ i ].dr_x[ 0 ] + face[ i ].dr_y[ 0 ] * face[ i ].dr_y[ 0 ] ) ;

				face[ i ].dr_x[ 1 ]	=	face[ i ].r[ 0 ] - face[ i ].cell[ 1 ]->r[ 0 ] ;
				face[ i ].dr_y[ 1 ]	=	face[ i ].r[ 1 ] - face[ i ].cell[ 1 ]->r[ 1 ] ;
				face[ i ].dr[ 1 ]	=	sqrt( face[ i ].dr_x[ 1 ] * face[ i ].dr_x[ 1 ] + face[ i ].dr_y[ 1 ] * face[ i ].dr_y[ 1 ] ) ;
			} else
			{
				dx	=	face[ i ].cell[ 0 ]->r[ 0 ] - face[ i ].r[ 0 ] ;
				dy	=	face[ i ].cell[ 0 ]->r[ 1 ] - face[ i ].r[ 1 ] ;

				face[ i ].dl		=	sqrt( dx * dx + dy * dy ) ;

				face[ i ].dr_x[ 0 ]	=	- dx ;
				face[ i ].dr_y[ 0 ]	=	- dy ;
				face[ i ].dr[ 0 ]	=	sqrt( face[ i ].dr_x[ 0 ] * face[ i ].dr_x[ 0 ] + face[ i ].dr_y[ 0 ] * face[ i ].dr_y[ 0 ] ) ;
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// New version
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			face[ i ].A[ 0 ]	=	face[ i ].node[ 1 ]->r[ 1 ] - face[ i ].node[ 0 ]->r[ 1 ]  ;
			face[ i ].A[ 1 ]	=	- ( face[ i ].node[ 1 ]->r[ 0 ] - face[ i ].node[ 0 ]->r[ 0 ] ) ;
			face[ i ].A[ 2 ]	=	0. ;

			face[ i ].dA		=	sqrt( face[ i ].A[ 0 ] * face[ i ].A[ 0 ] + face[ i ].A[ 1 ] * face[ i ].A[ 1 ] ) ;

			for ( j = 0 ; j < 2 ; j++ )
			{
				face[ i ].nA[ j ]	=	face[ i ].A[ j ] / face[ i ].dA ;
				if ( abs( face[ i ].nA[ j ] ) < 1e-10 )
					face[ i ].nA[ j ]	=	0. ;
			}
			face[ i ].nA[ 2 ]	=	0. ;

			face[ i ].nAd 		=	face[ i ].nA[ 0 ] * face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].nA[ 1 ] * face[ i ].node[ 1 ]->r[ 1 ] ;
			//face[ i ].nAd 	=	face[ i ].nAd / face[ i ].dA ;

			face[ i ].dl	=	0. ;
			for ( j = 0 ; j < 2 ; j++ )
			{
				if ( face[ i ].cell_number == 1 )
					face[ i ].r_cc[ j ]	=	face[ i ].r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;
				else
					face[ i ].r_cc[ j ]	=	face[ i ].cell[ 1 ]->r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;

				face[ i ].dl	+=	face[ i ].r_cc[ j ] * face[ i ].r_cc[ j ] ;
			}
			face[ i ].dl 	=	sqrt ( face[ i ].dl ) ;

			for ( j = 0 ; j < 2 ; j++ )
				face[ i ].nr_cc[ j ]	=	face[ i ].r_cc[ j ] / face[ i ].dl ;

			for ( j = 0 ; j < face[ i ].cell_number ; j++ )
			{
				face[ i ].dr_cf[ j ]	=	0. ;
				for ( k = 0 ; k < 2 ; k++ )
				{
					face[ i ].r_cf[ j ][ k ]	=	face[ i ].r[ k ] - face[ i ].cell[ j ]->r[ k ] ;
					face[ i ].dr_cf[ j ]		+=	face[ i ].r_cf[ j ][ k ] * face[ i ].r_cf[ j ][ k ] ;
				}
				face[ i ].dr_cf[ j ]	=	sqrt( face[ i ].dr_cf[ j ] ) ;
			}
		}
	} else if ( dimension == 3 )
	{
		for ( i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			face[ i ].A[ 0 ]	=	0. ;
			face[ i ].A[ 1 ]	=	0. ;
			face[ i ].A[ 2 ]	=	0. ;
			for ( j = 2 ; j < face[ i ].node_number ; j++ )
			{
				for ( k = 0 ; k < 3 ; k++ )
				{
					value[ 0 ][ k ]	=	face[ i ].node[ j - 1 ]->r[ k ] - face[ i ].node[ 0 ]->r[ k ] ;
					value[ 1 ][ k ]	=	face[ i ].node[ j ]->r[ k ] - face[ i ].node[ 0 ]->r[ k ] ;
				}

				face[ i ].A[ 0 ]	+=	value[ 0 ][ 1 ] * value[ 1 ][ 2 ] - value[ 0 ][ 2 ] * value[ 1 ][ 1 ] ;
				face[ i ].A[ 1 ]	+=	value[ 0 ][ 2 ] * value[ 1 ][ 0 ] - value[ 0 ][ 0 ] * value[ 1 ][ 2 ] ;
				face[ i ].A[ 2 ]	+=	value[ 0 ][ 0 ] * value[ 1 ][ 1 ] - value[ 0 ][ 1 ] * value[ 1 ][ 0 ] ;
			}
			face[ i ].A[ 0 ]	/=	2.0 ;
			face[ i ].A[ 1 ]	/=	2.0 ;
			face[ i ].A[ 2 ]	/=	2.0 ;

			face[ i ].dA		=	sqrt( face[ i ].A[ 0 ] * face[ i ].A[ 0 ] + face[ i ].A[ 1 ] * face[ i ].A[ 1 ] + face[ i ].A[ 2 ] * face[ i ].A[ 2 ] ) ;

			for ( j = 0 ; j < 3 ; j++ )
			{
				face[ i ].nA[ j ]	=	face[ i ].A[ j ] / face[ i ].dA ;
				if ( abs( face[ i ].nA[ j ] ) < 1e-10 )
					face[ i ].nA[ j ]	=	0. ;
			}

        	face[ i ].nAd 	=	face[ i ].nA[ 0 ] * face[ i ].node[ 1 ]->r[ 0 ] + face[ i ].nA[ 1 ] * face[ i ].node[ 1 ]->r[ 1 ] + face[ i ].nA[ 2 ] * face[ i ].node[ 1 ]->r[ 2 ] ;

			face[ i ].dl 	=	0. ;
			for ( j = 0 ; j < 3 ; j++ )
			{
				if ( face[ i ].cell_number == 1 )
					face[ i ].r_cc[ j ]	=	face[ i ].r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;
				else
					face[ i ].r_cc[ j ]	=	face[ i ].cell[ 1 ]->r[ j ] - face[ i ].cell[ 0 ]->r[ j ] ;

				face[ i ].dl	+=	face[ i ].r_cc[ j ] * face[ i ].r_cc[ j ] ;
			}
			face[ i ].dl	=	sqrt ( face[ i ].dl ) ;

			for ( j = 0 ; j < 3 ; j++ )
				face[ i ].nr_cc[ j ]	=	face[ i ].r_cc[ j ] / face[ i ].dl ;

			for ( j = 0 ; j < face[ i ].cell_number ; j++ )
			{
				face[ i ].dr_cf[ j ]	=	0. ;
				for ( k = 0 ; k < 3 ; k++ )
				{
					face[ i ].r_cf[ j ][ k ]	=	face[ i ].r[ k ] - face[ i ].cell[ j ]->r[ k ] ;
					face[ i ].dr_cf[ j ]		+=	face[ i ].r_cf[ j ][ k ] * face[ i ].r_cf[ j ][ k ] ;
				}
				face[ i ].dr_cf[ j ]	=	sqrt( face[ i ].dr_cf[ j ] ) ;
			}
		}
	}
}

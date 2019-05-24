#include "domain.h"


/*! \brief Calculate Face Coefficients

To obtain the cross point coordinate of face and line between cell centers, we can write the face (2D line) equation and line equation between cell centers. The line equation between cell centers is
\f[
	y - y_{ C_{0} } = m_{0} ( x - x_{ C_{0} } )
\f], and the line equation of 2D face is
\f[
	y - y_{ a } = m_{1} ( x - x_{ a } )
\f], where the slopes of lines are
\f[
	m_{0} = \frac{  y_{C_{1}} - y_{C_{0}}  }{  x_{C_{1}} - x_{C_{0}}  }
\f]
\f[
	m_{1} = \frac{  y_{b} - y_{a}  }{  x_{b} - x_{a} }
\f].
For the cases of \f$ x_{C_{1}} = x_{C_{0}} \f$ or \f$ x_{b} = x_{a} \f$, the slopes do not exsist and the line equation will be
\f[
	x = x_{C_{0}}
\f], or
\f[
	x = x_{a}
\f].
Once we have line equations and both \f$ m_{0} \f$ and \f$ m_{1} \f$ exist, we can have the coordinate of the cross point.
\f[
	x = \frac{ y_{a} - y_{C_{0}} + m_{0} x_{C_{0}} - m_{1} x_{a} } { m_{0} - m_{1} }
\f]
\f[
	y = \frac{ m_{0} y_{a} - m_{1} y_{C_{0}} + m_{0}m_{1} ( x_{C_{0}} - x_{a} ) } { m_{0} - m_{1} }
\f]



*/


void Domain::CalculateFaceCoefficient( )
{
	int	 	i, j ;
	double	x0, y0, z0, x1, y1, z1, xa, ya, xb, yb ;
	double	buffer_m0, buffer_m1 ;

	double	r_c[ 2 ][ 3 ], r_n[ 2 ][ 3 ] ;
	double	sum_dr_cf, L_cc[ 3 ], t, d ;

	if ( cylindrical_x == 1 || cylindrical_y == 1 )
	{
		for ( i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			if ( face[ i ].cell_number == 1 )
			{
				r_c[ 1 ][ 0 ]	=	face[ i ].r[ 0 ] ;
				r_c[ 1 ][ 1 ]	=	face[ i ].r[ 1 ] ;

				face[ i ].cf_alpha[ 0 ]	=	1. ;
				face[ i ].cf_alpha[ 1 ]	=	0. ;
			} else
			{
				r_c[ 1 ][ 0 ]	=	face[ i ].cell[ 1 ]->r[ 0 ] ;
				r_c[ 1 ][ 1 ]	=	face[ i ].cell[ 1 ]->r[ 1 ] ;

				sum_dr_cf				=	1. / face[ i ].dr_cf[ 0 ] + 1. / face[ i ].dr_cf[ 1 ] ;
				face[ i ].cf_alpha[ 0 ]	=	1. / face[ i ].dr_cf[ 0 ] / sum_dr_cf ;
				face[ i ].cf_alpha[ 1 ]	=	1. / face[ i ].dr_cf[ 1 ] / sum_dr_cf ;
			}

			r_c[ 0 ][ 0 ]	=	face[ i ].cell[ 0 ]->r[ 0 ] ;
			r_c[ 0 ][ 1 ]	=	face[ i ].cell[ 0 ]->r[ 1 ] ;

			r_n[ 0 ][ 0 ]	=	face[ i ].node[ 0 ]->r[ 0 ] ;
			r_n[ 0 ][ 1 ]	=	face[ i ].node[ 0 ]->r[ 1 ] ;
			r_n[ 1 ][ 0 ]	=	face[ i ].node[ 1 ]->r[ 0 ] ;
			r_n[ 1 ][ 1 ]	=	face[ i ].node[ 1 ]->r[ 1 ] ;

			// cross point and dl0 dl1
			if ( fabs ( r_c[ 0 ][ 0 ] - r_c[ 1 ][ 0 ] ) < 1.e-10 )
			{
				// m0 is vertical
				buffer_m0	=	0. ;
				buffer_m1	=	( r_n[ 1 ][ 1 ] - r_n[ 0 ][ 1 ] ) / ( r_n[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;

				face[ i ].cross_r[ 0 ]	=	r_c[ 0 ][ 0 ] ;
				face[ i ].cross_r[ 1 ]	=	r_n[ 0 ][ 1 ] +  buffer_m1 * ( face[ i ].cross_r[ 0 ] - r_n[ 0 ][ 0 ] ) ;
			} else if ( fabs ( r_n[ 0 ][ 0 ] - r_n[ 1 ][ 0 ] ) < 1.e-10   )
			{
				// m1 is vertical
				buffer_m0	= 	( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / ( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) ;
				buffer_m1	= 	0. ;

				face[ i ].cross_r[ 0 ]	= 	r_n[ 0 ][ 0 ] ;
				face[ i ].cross_r[ 1 ]	=	r_c[ 0 ][ 1 ] +  buffer_m0 * ( face[ i ].cross_r[ 0 ] -  r_c[ 0 ][ 0 ] ) ;
			} else
			{
				buffer_m0	=	( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / ( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) ;
				buffer_m1	=	( r_n[ 1 ][ 1 ] - r_n[ 0 ][ 1 ] ) / ( r_n[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;

				face[ i ].cross_r[ 0 ]	=	( r_n[ 0 ][ 1 ] - r_c[ 0 ][ 1 ] + buffer_m0 * r_c[ 0 ][ 0 ] - buffer_m1 * r_n[ 0 ][ 0 ] ) / ( buffer_m0 - buffer_m1 ) ;
				face[ i ].cross_r[ 1 ]	=	(  buffer_m0 * r_n[ 0 ][ 1 ] - buffer_m1 * r_c[ 0 ][ 1 ] + buffer_m0 * buffer_m1 * ( r_c[ 0 ][ 0 ] - r_n[ 0 ][ 0 ] ) ) / ( buffer_m0 - buffer_m1 ) ;
			}

			face[ i ].r_fof[ 0 ]	=	face[ i ].cross_r[ 0 ] - face[ i ].r[ 0 ] ;
			face[ i ].r_fof[ 1 ]	=	face[ i ].cross_r[ 1 ] - face[ i ].r[ 1 ] ;

			face[ i ].pdl[ 0 ]		=	sqrt( ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) * ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) + ( face[ i ].cross_r[ 1 ] - r_c[ 0 ][ 1 ] ) * ( face[ i ].cross_r[ 1 ] - r_c[ 0 ][ 1 ] ) ) ;
			face[ i ].pdl[ 1 ]		=	sqrt( ( face[ i ].cross_r[ 0 ] - r_c[ 1 ][ 0 ] ) * ( face[ i ].cross_r[ 0 ] - r_c[ 1 ][ 0 ] ) + ( face[ i ].cross_r[ 1 ] - r_c[ 1 ][ 1 ] ) * ( face[ i ].cross_r[ 1 ] - r_c[ 1 ][ 1 ] ) ) ;

			// dl_cos[ 0 ] is the shortest distance from cell 0 to the face.
			// dl_cos[ 1 ] is the shortest distance from cell 1 to the face.
			// The line equation of face is y - ya = buffer_m1 * ( x - x0 )
			// y - buffer_m1 * x - ( ya - buffer_m1 * x0 ) = 0
			// The distance from a point (x, y) to the line is d = abs ( y - buffer_m1 * x - (ya - buffer_m1 * x0 )  ) / sqrt ( buffer_m1^2 + 1 )

			if ( fabs  ( r_n[ 0 ][ 0 ] - r_n[ 1 ][ 0 ] )  < 1.e-10 )
			{
				face[ i ].dl_cos[ 0 ]	=	fabs ( r_c[ 0 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;
			} else
			{
				face[ i ].dl_cos[ 0 ]	=	fabs ( r_c[ 0 ][ 1 ] - buffer_m1 * r_c[ 0 ][ 0 ] - ( r_n[ 0 ][ 1 ] - buffer_m1 * r_n[ 0 ][ 0 ] ) ) / sqrt ( buffer_m1 * buffer_m1 + 1. ) ;
			}

			if ( face[ i ].cell_number > 1 )
			{
				if ( buffer_m1 == 0. )
					face[ i ].dl_cos[ 1 ]	=	fabs ( r_c[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;
				else
					face[ i ].dl_cos[ 1 ]	=	fabs ( r_c[ 1 ][ 1 ] - buffer_m1 * r_c[ 1 ][ 0 ] - ( r_n[ 0 ][ 1 ] - buffer_m1 * r_n[ 0 ][ 0 ] ) ) / sqrt ( buffer_m1 * buffer_m1 + 1. ) ;
			}

			if ( face[ i ].cell_number == 1 )
			{
				face[ i ].cc_alpha[ 0 ]	=	1.0 ;
				face[ i ].cc_alpha[ 1 ]	=	0. ;
			} else
			{
				face[ i ].cc_alpha[ 0 ]	=	face[ i ].pdl[ 1 ] / face[ i ].dl ;
				face[ i ].cc_alpha[ 1 ]	=	1.0 - face[ i ].cc_alpha[ 0 ] ;
			}

			if ( face[ i ].cell_number == 1 )
			{
				face[ i ].diffusion_coefficient_1st	=	face[ i ].dA / face[ i ].dl_cos[ 0 ] ;
			} else
			{
				face[ i ].diffusion_coefficient_1st	=	face[ i ].dA / face[ i ].dl ;
			}

			face[ i ].diffusion_coefficient_2nd[ 0 ][ 0 ]	=	face[ i ].dA * ( face[ i ].cell_sign[ 0 ] * face[ i ].nA[ 0 ] - ( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) / face[ i ].dl ) ;
			face[ i ].diffusion_coefficient_2nd[ 0 ][ 1 ]	=	face[ i ].dA * ( face[ i ].cell_sign[ 0 ] * face[ i ].nA[ 1 ] - ( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / face[ i ].dl ) ;
			face[ i ].diffusion_coefficient_2nd[ 0 ][ 2 ]	=	0. ;

			face[ i ].diffusion_coefficient_2nd[ 1 ][ 0 ]	=	- face[ i ].diffusion_coefficient_2nd[ 0 ][ 0 ] ;
			face[ i ].diffusion_coefficient_2nd[ 1 ][ 1 ]	=	- face[ i ].diffusion_coefficient_2nd[ 0 ][ 1 ] ;
			face[ i ].diffusion_coefficient_2nd[ 1 ][ 2 ]	=	0. ;
		}
	} else if ( dimension == 2 )
	{
		for ( int i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			if ( face[ i ].cell_number == 1 )
			{
				x1	=	face[ i ].x ;
				y1	=	face[ i ].y ;
			} else
			{
				x1	=	face[ i ].cell[ 1 ]->x ;
				y1	=	face[ i ].cell[ 1 ]->y ;
			}

			x0	=	face[ i ].cell[ 0 ]->x ;
			y0	=	face[ i ].cell[ 0 ]->y ;

			xa 	=	face[ i ].node[ 0 ]->x ;
			ya	=	face[ i ].node[ 0 ]->y ;
			xb	=	face[ i ].node[ 1 ]->x ;
			yb	=	face[ i ].node[ 1 ]->y ;

			//face[ i ].x_xi			=	( x1 - x0 ) / face[ i ].dl ;
			//face[ i ].y_xi			=	( y1 - y0 ) / face[ i ].dl ;

			//face[ i ].x_eta			=	( xb - xa ) / face[ i ].dA ;
			//face[ i ].y_eta			=	( yb - ya ) / face[ i ].dA ;

			//face[ i ].xi_coefficient	=	( face[ i ].Ax * face[ i ].Ax + face[ i ].Ay * face[ i ].Ay ) / ( face[ i ].Ax * face[ i ].x_xi + face[ i ].Ay * face[ i ].y_xi ) ;

			//face[ i ].eta_coefficient	=	face[ i ].xi_coefficient * ( face[ i ].x_xi * face[ i ].x_eta + face[ i ].y_xi * face[ i ].y_eta ) ;

			//face[ i ].A_dot_xi		=	face[ i ].Ax * face[ i ].x_xi + face[ i ].Ay * face[ i ].y_xi ;

			//if ( fabs( face[ i ].xi_coefficient ) < 1.e-10 )
			//	face[ i ].xi_coefficient	=	0. ;
			//if ( fabs( face[ i ].eta_coefficient ) < 1.e-10 )
			//	face[ i ].eta_coefficient	=	0. ;

			// cross point and dl0 dl1
			if (  fabs ( x0 - x1 )  < 1.e-10 )
			{
				// m0 is vertical
				buffer_m0			=	0. ;
				buffer_m1			=	( yb - ya ) / ( xb - xa ) ;

				face[ i ].cross_x	=	x0 ;
				face[ i ].cross_y	=	ya +  buffer_m1 * ( face[ i ].cross_x - xa ) ;
			} else if ( fabs  ( xa - xb )  < 1.e-10 )
			{
				// m1 is vertical
				buffer_m0			=	( y1 - y0 ) / ( x1 - x0 ) ;
				buffer_m1			=	0. ;

				face[ i ].cross_x	=	xa ;
				face[ i ].cross_y	=	y0 +  buffer_m0 * ( face[ i ].cross_x - x0 ) ;
			} else
			{
				buffer_m0			=	( y1 - y0 ) / ( x1 - x0 ) ;
				buffer_m1			=	( yb - ya ) / ( xb - xa ) ;

				face[ i ].cross_x	=	( ya - y0 + buffer_m0 * x0 - buffer_m1 * xa) / ( buffer_m0 - buffer_m1 ) ;
				face[ i ].cross_y	=	( buffer_m0 * ya - buffer_m1 * y0 + buffer_m0 * buffer_m1 * ( x0 - xa) ) / ( buffer_m0 - buffer_m1 ) ;
			}
			face[ i ].pdl[ 0 ]		=	sqrt( ( face[ i ].cross_x - x0 ) * ( face[ i ].cross_x - x0 ) + ( face[ i ].cross_y - y0 ) * ( face[ i ].cross_y - y0 ) ) ;
			//cout << "FID " << i << "\t" << setprecision(15) << face[ i ].cross_x << "\t" << setprecision(15) << face[ i ].cross_y << "\t" << setprecision(15) << x0 << "\t" << setprecision(15) << y0 <<  endl;
			face[ i ].pdl[ 1 ]		=	sqrt( ( face[ i ].cross_x - x1 ) * ( face[ i ].cross_x - x1 ) + ( face[ i ].cross_y - y1 ) * ( face[ i ].cross_y - y1 ) ) ;

			// dl0_cos is the shortest distance from cell 0 to the face.
			// dl1_cos is the shortest distance from cell 1 to the face.
			// The line equation of face is y - ya = buffer_m1 * ( x - xa )
			// y - buffer_m1 * x - (ya - buffer_m1 * xa ) = 0
			// The distance from a point (x, y) to the line is d = abs ( y - buffer_m1 * x - (ya - buffer_m1 * xa )  ) / sqrt ( buffer_m1^2 + 1 )

			if ( fabs  ( xa - xb )  < 1.e-10 )
				face[ i ].dl0_cos	=	fabs ( x0 - xa ) ;
			else
				face[ i ].dl0_cos	=	fabs ( y0 - buffer_m1 * x0 - ( ya - buffer_m1 * xa ) ) / sqrt ( buffer_m1 * buffer_m1 + 1. ) ;

			if ( face[ i ].cell_number > 1 )
			{
				if ( buffer_m1 == 0. )
					face[ i ].dl1_cos	=	fabs ( x1 - xa ) ;
				else
					face[ i ].dl1_cos	=	fabs ( y1 - buffer_m1 * x1 - ( ya - buffer_m1 * xa ) ) / sqrt ( buffer_m1 * buffer_m1 + 1. ) ;
			}

			//cout << face[ i ].dl0 << "\t" << face[ i ].dl0_cos << "\t" << face[ i ].dl0 - face[ i ].dl0_cos << "\t" << buffer_m1 << "\t" << ya << "\t" << y0 << endl ;
			//if ( face[ i ].cell_number > 1 )
			//	cout << face[ i ].dl1 << "\t" << face[ i ].dl1_cos << "\t" << face[ i ].dl1 - face[ i ].dl1_cos   << endl ;

			//cout << "plot \"< echo \'" << x0 << "\t" << y0 << "\'\" pt 7 ps 2 notitle " ;
			//cout << ", \"< echo \'" << x1 << "\t" << y1 << "\'\" pt 7 ps 2 notitle " ;
			//cout << ", \"< echo \'" << xa << "\t" << ya << "\'\" pt 7 ps 2 notitle " ;
			//cout << ", \"< echo \'" << xb << "\t" << yb << "\'\" pt 7 ps 2 notitle " ;
			//cout << ", \"< echo \'" << face[ i ].cross_x << "\t" << face[ i ].cross_y << "\'\" pt 7 ps 2 notitle " << endl ;
			//cout << endl;

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// New version
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if ( face[ i ].cell_number == 1 )
			{
				r_c[ 1 ][ 0 ]	=	face[ i ].r[ 0 ] ;
				r_c[ 1 ][ 1 ]	=	face[ i ].r[ 1 ] ;

				face[ i ].cf_alpha[ 0 ]	=	1. ;
				face[ i ].cf_alpha[ 1 ]	=	0. ;
			} else
			{
				r_c[ 1 ][ 0 ]	=	face[ i ].cell[ 1 ]->r[ 0 ] ;
				r_c[ 1 ][ 1 ]	=	face[ i ].cell[ 1 ]->r[ 1 ] ;

				sum_dr_cf					=	1. / face[ i ].dr_cf[ 0 ] + 1. / face[ i ].dr_cf[ 1 ] ;
				face[ i ].cf_alpha[ 0 ]	=	1. / face[ i ].dr_cf[ 0 ] / sum_dr_cf ;
				face[ i ].cf_alpha[ 1 ]	=	1. / face[ i ].dr_cf[ 1 ] / sum_dr_cf ;
			}

			r_c[ 0 ][ 0 ]	=	face[ i ].cell[ 0 ]->r[ 0 ] ;
			r_c[ 0 ][ 1 ]	=	face[ i ].cell[ 0 ]->r[ 1 ] ;

			r_n[ 0 ][ 0 ]	=	face[ i ].node[ 0 ]->r[ 0 ] ;
			r_n[ 0 ][ 1 ]	=	face[ i ].node[ 0 ]->r[ 1 ] ;
			r_n[ 1 ][ 0 ]	=	face[ i ].node[ 1 ]->r[ 0 ] ;
			r_n[ 1 ][ 1 ]	=	face[ i ].node[ 1 ]->r[ 1 ] ;

			//face[ i ]._xi[ 0 ]		=	( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) / face[ i ].dl ;
			//face[ i ]._xi[ 1 ]		=	( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / face[ i ].dl ;

			//face[ i ]._eta[ 0 ]		=	( r_n[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) / face[ i ].dA ;
			//face[ i ]._eta[ 1 ]		=	( r_n[ 1 ][ 1 ] - r_n[ 0 ][ 1 ] ) / face[ i ].dA ;

			//face[ i ].A_dot_xi		=	face[ i ].A[ 0 ] * face[ i ]._xi[ 0 ] + face[ i ].A[ 1 ] * face[ i ]._xi[ 1 ] ;
			//face[ i ].xi_coefficient	=	( face[ i ].A[ 0 ] * face[ i ].A[ 0 ] + face[ i ].A[ 1 ] * face[ i ].A[ 1 ] ) / face[ i ].A_dot_xi ;
			//face[ i ].eta_coefficient	=	face[ i ].xi_coefficient * ( face[ i ]._xi[ 0 ] * face[ i ]._eta[ 0 ] + face[ i ]._xi[ 1 ] * face[ i ]._eta[ 1 ] ) ;

			//if ( fabs( face[ i ].xi_coefficient ) < 1.e-10 )
			//	face[ i ].xi_coefficient	=	0. ;
			//if ( fabs( face[ i ].eta_coefficient ) < 1.e-10 )
			//	face[ i ].eta_coefficient	=	0. ;

			// cross point and dl0 dl1
			if (  fabs ( r_c[ 0 ][ 0 ] - r_c[ 1 ][ 0 ] )  < 1.e-10 )
			{
				// m0 is vertical
				buffer_m0	=	0. ;
				buffer_m1	=	( r_n[ 1 ][ 1 ] - r_n[ 0 ][ 1 ] ) / ( r_n[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;

				face[ i ].cross_r[ 0 ]	=	r_c[ 0 ][ 0 ] ;
				face[ i ].cross_r[ 1 ]	=	r_n[ 0 ][ 1 ] +  buffer_m1 * ( face[ i ].cross_r[ 0 ] - r_n[ 0 ][ 0 ] ) ;
			} else if ( fabs  ( r_n[ 0 ][ 0 ] - r_n[ 1 ][ 0 ] )  < 1.e-10   )
			{
				// m1 is vertical
				buffer_m0	= 	( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / ( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) ;
				buffer_m1	= 	0. ;

				face[ i ].cross_r[ 0 ]	= 	r_n[ 0 ][ 0 ] ;
				face[ i ].cross_r[ 1 ]	=	r_c[ 0 ][ 1 ] +  buffer_m0 * ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) ;
			} else
			{
				buffer_m0	=	( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / ( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) ;
				buffer_m1	=	( r_n[ 1 ][ 1 ] - r_n[ 0 ][ 1 ] ) / ( r_n[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;

				face[ i ].cross_r[ 0 ]	=	( r_n[ 0 ][ 1 ] - r_c[ 0 ][ 1 ] + buffer_m0 * r_c[ 0 ][ 0 ] - buffer_m1 * r_n[ 0 ][ 0 ] ) / ( buffer_m0 - buffer_m1 ) ;
				face[ i ].cross_r[ 1 ]	=	(  buffer_m0 * r_n[ 0 ][ 1 ] - buffer_m1 * r_c[ 0 ][ 1 ] + buffer_m0 * buffer_m1 * ( r_c[ 0 ][ 0 ] - r_n[ 0 ][ 0 ] ) ) / ( buffer_m0 - buffer_m1 ) ;
			}

			face[ i ].r_fof[ 0 ]	=	face[ i ].cross_r[ 0 ] - face[ i ].r[ 0 ] ;
			face[ i ].r_fof[ 1 ]	=	face[ i ].cross_r[ 1 ] - face[ i ].r[ 1 ] ;

			face[ i ].pdl[ 0 ]		=	sqrt( ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) * ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) + ( face[ i ].cross_r[ 1 ] - r_c[ 0 ][ 1 ] ) * ( face[ i ].cross_r[ 1 ] - r_c[ 0 ][ 1 ] ) ) ;
			face[ i ].pdl[ 1 ]		=	sqrt( ( face[ i ].cross_r[ 0 ] - r_c[ 1 ][ 0 ] ) * ( face[ i ].cross_r[ 0 ] - r_c[ 1 ][ 0 ] ) + ( face[ i ].cross_r[ 1 ] - r_c[ 1 ][ 1 ] ) * ( face[ i ].cross_r[ 1 ] - r_c[ 1 ][ 1 ] ) ) ;

			// dl_cos[ 0 ] is the shortest distance from cell 0 to the face.
			// dl_cos[ 1 ] is the shortest distance from cell 1 to the face.
			// The line equation of face is y - ya = buffer_m1 * ( x - x0 )
			// y - buffer_m1 * x - ( ya - buffer_m1 * x0 ) = 0
			// The distance from a point (x, y) to the line is d = abs ( y - buffer_m1 * x - (ya - buffer_m1 * x0 )  ) / sqrt ( buffer_m1^2 + 1 )

			if ( fabs  ( r_n[ 0 ][ 0 ] - r_n[ 1 ][ 0 ] )  < 1.e-10 )
				face[ i ].dl_cos[ 0 ]	=	fabs ( r_c[ 0 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;
			else
				face[ i ].dl_cos[ 0 ]	=	fabs ( r_c[ 0 ][ 1 ] - buffer_m1 * r_c[ 0 ][ 0 ] - ( r_n[ 0 ][ 1 ] - buffer_m1 * r_n[ 0 ][ 0 ] )  ) / sqrt ( buffer_m1 * buffer_m1 + 1. ) ;

			if ( face[ i ].cell_number > 1 )
			{
				if ( buffer_m1 == 0. )
					face[ i ].dl_cos[ 1 ]	=	fabs ( r_c[ 1 ][ 0 ] - r_n[ 0 ][ 0 ] ) ;
				else
					face[ i ].dl_cos[ 1 ]	=	fabs ( r_c[ 1 ][ 1 ] - buffer_m1 * r_c[ 1 ][ 0 ] - ( r_n[ 0 ][ 1 ] - buffer_m1 * r_n[ 0 ][ 0 ] )  ) / sqrt ( buffer_m1 * buffer_m1 + 1. ) ;
			}

			if ( face[ i ].cell_number == 1 )
			{
				face[ i ].cc_alpha[ 0 ]	=	1.0 ;
				face[ i ].cc_alpha[ 1 ]	=	0. ;
			} else
			{
				face[ i ].cc_alpha[ 0 ]	=	face[ i ].pdl[ 1 ] / face[ i ].dl ;
				face[ i ].cc_alpha[ 1 ]	=	1.0 - face[ i ].cc_alpha[ 0 ] ;
			}

			face[ i ].diffusion_coefficient_1st				=	face[ i ].dA / face[ i ].dl ;
			face[ i ].diffusion_coefficient_2nd[ 0 ][ 0 ]	=	face[ i ].dA * ( face[ i ].cell_sign[ 0 ] * face[ i ].nA[ 0 ] - ( r_c[ 1 ][ 0 ] - r_c[ 0 ][ 0 ] ) / face[ i ].dl ) ;
			face[ i ].diffusion_coefficient_2nd[ 0 ][ 1 ]	=	face[ i ].dA * ( face[ i ].cell_sign[ 0 ] * face[ i ].nA[ 1 ] - ( r_c[ 1 ][ 1 ] - r_c[ 0 ][ 1 ] ) / face[ i ].dl ) ;
			face[ i ].diffusion_coefficient_2nd[ 0 ][ 2 ]	=	0. ;

			face[ i ].diffusion_coefficient_2nd[ 1 ][ 0 ]	=	- face[ i ].diffusion_coefficient_2nd[ 0 ][ 0 ] ;
			face[ i ].diffusion_coefficient_2nd[ 1 ][ 1 ]	=	- face[ i ].diffusion_coefficient_2nd[ 0 ][ 1 ] ;
			face[ i ].diffusion_coefficient_2nd[ 1 ][ 2 ]	=	0. ;
		}
	} else
	{
		for ( int i = 0 ; i < local_face_number + ghost_face_number_level_1 ; i++ )
		{
			if ( face[ i ].cell_number == 1 )
			{
				r_c[ 1 ][ 0 ]	=	face[ i ].r[ 0 ] ;
				r_c[ 1 ][ 1 ]	=	face[ i ].r[ 1 ] ;
				r_c[ 1 ][ 2 ]	=	face[ i ].r[ 2 ] ;

				face[ i ].cf_alpha[ 0 ]	=	1. ;
				face[ i ].cf_alpha[ 1 ]	=	0. ;
			} else
			{
				r_c[ 1 ][ 0 ]	=	face[ i ].cell[ 1 ]->r[ 0 ] ;
				r_c[ 1 ][ 1 ]	=	face[ i ].cell[ 1 ]->r[ 1 ] ;
				r_c[ 1 ][ 2 ]	=	face[ i ].cell[ 1 ]->r[ 2 ] ;

				sum_dr_cf				=	( 1. / face[ i ].dr_cf[ 0 ] + 1. / face[ i ].dr_cf[ 1 ] ) ;
				face[ i ].cf_alpha[ 0 ]	=	1. / face[ i ].dr_cf[ 0 ] / sum_dr_cf ;
				face[ i ].cf_alpha[ 1 ]	=	1. / face[ i ].dr_cf[ 1 ] / sum_dr_cf ;
			}

			//face[ i ].A_dot_xi		=	0. ;
			for ( j = 0 ; j < 3 ; j++ )
			{
				r_c[ 0 ][ j ]	=	face[ i ].cell[ 0 ]->r[ j ] ;
				L_cc[ j ]		=	r_c[ 1 ][ j ] - r_c[ 0 ][ j ] ;
				//face[ i ]._xi[ j ]	=	L_cc[ j ] / face[ i ].dl ;
				//face[ i ].A_dot_xi	+=	face[ i ].A[ j ] * face[ i ]._xi[ j ] ;
			}
			//face[ i ].xi_coefficient	=	face[ i ].dA * face[ i ].dA / face[ i ].A_dot_xi ;

			//if ( fabs( face[ i ].xi_coefficient ) < 1.e-10 )
			//	face[ i ].xi_coefficient	=	0. ;

			// cross point and dl0 dl1
			d	=	- ( face[ i ].A[ 0 ] * face[ i ].r[ 0 ] + face[ i ].A[ 1 ] * face[ i ].r[ 1 ] + face[ i ].A[ 2 ] * face[ i ].r[ 2 ] ) ;
			t	=	- ( d + ( face[ i ].A[ 0 ] * r_c[ 0 ][ 0 ]  + face[ i ].A[ 1 ] * r_c[ 0 ][ 1 ] + face[ i ].A[ 2 ] * r_c[ 0 ][ 2 ] ) ) / ( face[ i ].A[ 0 ] * L_cc[ 0 ] + face[ i ].A[ 1 ] * L_cc[ 1 ] + face[ i ].A[ 2 ] * L_cc[ 2 ] ) ;

			face[ i ].cross_r[ 0 ]	=	r_c[ 0 ][ 0 ] + L_cc[ 0 ] * t ;
			face[ i ].cross_r[ 1 ]	=	r_c[ 0 ][ 1 ] + L_cc[ 1 ] * t ;
			face[ i ].cross_r[ 2 ]	=	r_c[ 0 ][ 2 ] + L_cc[ 2 ] * t ;

			face[ i ].r_fof[ 0 ]	=	face[ i ].cross_r[ 0 ] - face[ i ].r[ 0 ] ;
			face[ i ].r_fof[ 1 ]	=	face[ i ].cross_r[ 1 ] - face[ i ].r[ 1 ] ;
			face[ i ].r_fof[ 2 ]	=	face[ i ].cross_r[ 2 ] - face[ i ].r[ 2 ] ;

			face[ i ].pdl[ 0 ]		=	sqrt( ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) * ( face[ i ].cross_r[ 0 ] - r_c[ 0 ][ 0 ] ) + ( face[ i ].cross_r[ 1 ] - r_c[ 0 ][ 1 ] ) * ( face[ i ].cross_r[ 1 ] - r_c[ 0 ][ 1 ] ) + ( face[ i ].cross_r[ 2 ] - r_c[ 0 ][ 2 ] ) * ( face[ i ].cross_r[ 2 ] - r_c[ 0 ][ 2 ] ) ) ;
			face[ i ].pdl[ 1 ]		=	sqrt( ( face[ i ].cross_r[ 0 ] - r_c[ 1 ][ 0 ] ) * ( face[ i ].cross_r[ 0 ] - r_c[ 1 ][ 0 ] ) + ( face[ i ].cross_r[ 1 ] - r_c[ 1 ][ 1 ] ) * ( face[ i ].cross_r[ 1 ] - r_c[ 1 ][ 1 ] ) + ( face[ i ].cross_r[ 2 ] - r_c[ 1 ][ 2 ] ) * ( face[ i ].cross_r[ 2 ] - r_c[ 1 ][ 2 ] ) ) ;

			face[ i ].dl_cos[ 0 ]	=	fabs( face[ i ].A[ 0 ] * r_c[ 0 ][ 0 ]  + face[ i ].A[ 1 ] * r_c[ 0 ][ 1 ] + face[ i ].A[ 2 ] * r_c[ 0 ][ 2 ] + d ) / face[ i ].dA ;
			face[ i ].dl_cos[ 1 ]	=	fabs( face[ i ].A[ 0 ] * r_c[ 1 ][ 0 ]  + face[ i ].A[ 1 ] * r_c[ 1 ][ 1 ] + face[ i ].A[ 2 ] * r_c[ 1 ][ 2 ] + d ) / face[ i ].dA ;

			if ( face[ i ].cell_number == 1 )
			{
				face[ i ].cc_alpha[ 0 ]	=	1. ;
				face[ i ].cc_alpha[ 1 ]	=	0. ;
			} else
			{
				face[ i ].cc_alpha[ 0 ]	=	face[ i ].pdl[ 1 ] / face[ i ].dl ;
				face[ i ].cc_alpha[ 1 ]	=	1.0 - face[ i ].cc_alpha[ 0 ] ;
			}

			face[ i ].diffusion_coefficient_1st	=	face[ i ].dA / face[ i ].dl ;

			for ( j = 0 ; j < 3 ; j++ )
			{
				face[ i ].diffusion_coefficient_2nd[ 0 ][ j ]	=	face[ i ].dA * ( face[ i ].cell_sign[ 0 ] * face[ i ].nA[ j ] - L_cc[ j ] / face[ i ].dl ) ;
				face[ i ].diffusion_coefficient_2nd[ 1 ][ j ]	=	- face[ i ].diffusion_coefficient_2nd[ 0 ][ j ] ;
			}
		}
	}
}
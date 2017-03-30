#include "spe_groundsurface.h"
#include "spe_triangularization.h"
#include "gauss_legendre.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CONVERGENCE_PRECISION 1e-6

#ifndef PI
#define PI 3.14159265358979323846
#endif

double ground_surface_integral_over_r( double r_y, void *params ) {
	unsigned int finalorder;
	double *d_params = (double *)params;
	double new_params[] = {
		r_y,
		d_params[ 0 ],
		d_params[ 1 ],
		d_params[ 2 ],
		d_params[ 3 ],
		d_params[ 4 ],
		d_params[ 5 ],
		d_params[ 6 ],
		d_params[ 7 ],
		d_params[ 8 ],
		d_params[ 9 ],
		d_params[ 10 ]
	};
	
	return gauss_legendre_eval_converge_ab( &ground_surface_integral_over_theta, new_params, 0, 2*PI, 
		CONVERGENCE_PRECISION, &finalorder );
}

double ground_surface_integral_over_theta( double theta_y, void *params ) {
	
	double r_y, z_y, r_x, theta_x, phi_x, k, R, dhdy1, dhdy2, phi_1, phi_2, phi_3;
	double sph_r_y, sph_theta_y, sph_phi_y;
	double *d_params, phi, G;
	
	/* Cast the params to the proper type */
	d_params = (double *)params;
	
	/* put the variables in their proper places for readability */
	r_y = d_params[ 0 ];
	z_y = d_params[ 1 ];
	r_x = d_params[ 2 ];
	theta_x = d_params[ 3 ];
	phi_x = d_params[ 4 ];
	k = d_params[ 5 ];
	R = d_params[ 6 ];
	dhdy1 = d_params[ 7 ];
	dhdy2 = d_params[ 8 ];
	phi_1 = d_params[ 9 ];
	phi_2 = d_params[ 10 ];
	phi_3 = d_params[ 11 ];
	
	
	/* Convert y from cylindrical to spherical coordinates for script_G calculation */
	sph_r_y = sqrt( r_y*r_y + z_y*z_y );
	sph_phi_y = theta_y;
	sph_theta_y = atan2( r_y / z_y );
	
	/* calculate phi(r, theta) */
	phi = phi_1 + phi_2 * r_y + phi_3 * theta_y;
	
	/* Calculate script_G */
	G = script_G( r_x, theta_x, phi_x, sph_r_y, sph_theta_y, sph_phi_y, k, R, dhdy1, dhdy2 );
	
	/* return the final product */
	return G * r_y * phi;
}
#include "gauss-legendre.h"
#include <math.h>
#include "utilities.h"

double f_internal( double r, void *params ) {
	/* function is r^3 * sin(2*theta).  Here theta is fixed. */
	double *dp = (double *)params;
	return r*r*r*sin( 2*dp[ 0 ] );
}

double f_external( double theta, void *params ) {
	unsigned int order = 2;
	double *dp = (double *)params;
	double new_params[] = { theta };
	return gauss_legendre_eval_converge_ab( &f_internal, new_params, dp[ 0 ], dp[ 1 ], 1e-8, &order );
}

int main( int argc, char **argv ) {
	unsigned int order = 2;
	double params[] = { 2.0, 5.0 };
	
	printf( "integral is %.5f\n", gauss_legendre_eval_converge_ab( &f_external, params, 0.0, PI/2, 1e-8, &order) );
	printf( "Actual value is %.5f\n", 609.0 / 4.0 );
		
}
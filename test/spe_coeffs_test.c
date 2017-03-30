/**
Constructs the a,b,c coefficient vector for a given set of r and theta, and compares
the resulting u equations to those calculated by hand for a given j,k vertex.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilities.h"
#include "spe_triangularization.h"

typedef enum { false, true } bool;

bool isCloseEnough( double a, double b );
bool checkABC( double a1, double b1, double c1, double a2, double b2, double b3 );

int main( int argc, char **argv ) {
	
	// Set up for r = { 0, 1, 2, 3, 4 } and theta = {0, pi/2, pi, 3pi/2, 2pi }
	unsigned int nr, ntheta;
	double rmax, psi[ 3 ];
	coeffVector vector;
	int j, k, pm, m, abc;
	
	nr = 5;
	ntheta = 5;
	rmax = 4.0;
	
	if (!init_coefficient_storage( nr, ntheta, &vector )) {
		printf("Can't initialize coefficient storage!\n");
		return 1;
	}
	calculate_coefficients( rmax, &vector );
	
	/* Test a couple of points.  Start with j=0,k=0 => r=0,theta=0 */
	
	/* For t_00+, u1+ = 1-r, u2+ = r-2theta/pi, u3+ = 2theta/pi
	   For t_00-, u1- = 1-2theta/pi, u2- = r, u3- = 2theta/pi-r */
	if (!checkABC( 1, -1, 0, 
			get_abc_coeff( &vector, 0, 0, 0, 0, 0 ),
			get_abc_coeff( &vector, 0, 0, 0, 0, 1 ),
			get_abc_coeff( &vector, 0, 0, 0, 0, 2 )
			)) {
		printf("Coefficients don't match up for u_00+;1\n");
	} else {
		printf("Coefficients OK for u_00+;1\n");
	}
	if (!checkABC( 0, 1, -2/PI, 
			get_abc_coeff( &vector, 0, 0, 0, 1, 0 ),
			get_abc_coeff( &vector, 0, 0, 0, 1, 1 ),
			get_abc_coeff( &vector, 0, 0, 0, 1, 2 )
			)) {
		printf("Coefficients don't match up for u_00+;2\n");
	} else {
		printf("Coefficients OK for u_00+;2\n");
	}
	if (!checkABC( 0, 0, 2/PI,
			get_abc_coeff( &vector, 0, 0, 0, 2, 0 ),
			get_abc_coeff( &vector, 0, 0, 0, 2, 1 ),
			get_abc_coeff( &vector, 0, 0, 0, 2, 2 ) 
			)) {
		printf("Coefficients don't match up for u_00+;3\n");
	} else {
		printf("Coefficients OK for u_00+;3\n");
	}
	
	if (!checkABC( 1, 0, -2/PI, 
			get_abc_coeff( &vector, 0, 0, 1, 0, 0 ),
			get_abc_coeff( &vector, 0, 0, 1, 0, 1 ),
			get_abc_coeff( &vector, 0, 0, 1, 0, 2 )
			)) {
		printf("Coefficients don't match up for u_00-;1\n");
	} else {
		printf("Coefficients OK for u_00-;1\n");
	}
	if (!checkABC( 0, 1, 0, 
			get_abc_coeff( &vector, 0, 0, 1, 1, 0 ),
			get_abc_coeff( &vector, 0, 0, 1, 1, 1 ),
			get_abc_coeff( &vector, 0, 0, 1, 1, 2 )
			)) {
		printf("Coefficients don't match up for u_00-;2\n");
	} else {
		printf("Coefficients OK for u_00-;2\n");
	}
	if (!checkABC( 0, -1, 2/PI,
			get_abc_coeff( &vector, 0, 0, 1, 2, 0 ),
			get_abc_coeff( &vector, 0, 0, 1, 2, 1 ),
			get_abc_coeff( &vector, 0, 0, 1, 2, 2 ) 
			)) {
		printf("Coefficients don't match up for u_00-;3\n");
	} else {
		printf("Coefficients OK for u_00-;3\n");
	}
	
	
	/* For t_32+, u1+ = 4-r, u2+ = r-1-2theta/pi, u3+ = 2theta/pi-2
	   For t_32-, u1- = 3-2theta/pi, u2- = r-3, u3- = 1-r+2theta/pi */
	if (!checkABC( 4, -1, 0,
			get_abc_coeff( &vector, 3, 2, 0, 0, 0 ),
			get_abc_coeff( &vector, 3, 2, 0, 0, 1 ),
			get_abc_coeff( &vector, 3, 2, 0, 0, 2 ) 
			)) {
		printf("Coefficients don't match up for u_32+;1\n");
	} else {
		printf("Coefficients OK for u_32+;1\n");
	}
	if (!checkABC( -1, 1, -2/PI, 
			get_abc_coeff( &vector, 3, 2, 0, 1, 0 ),
			get_abc_coeff( &vector, 3, 2, 0, 1, 1 ),
			get_abc_coeff( &vector, 3, 2, 0, 1, 2 ) 
			)) {
		printf("Coefficients don't match up for u_32+;2\n");
	} else {
		printf("Coefficients OK for u_32+;2\n");
	}
	if (!checkABC( -2, 0, 2/PI, 
			get_abc_coeff( &vector, 3, 2, 0, 2, 0 ),
			get_abc_coeff( &vector, 3, 2, 0, 2, 1 ),
			get_abc_coeff( &vector, 3, 2, 0, 2, 2 ) 
			)) {
		printf("Coefficients don't match up for u_32+;3\n");
	} else {
		printf("Coefficients OK for u_32+;3\n");
	}
	
	if (!checkABC( 3, 0, -2/PI, 
			get_abc_coeff( &vector, 3, 2, 1, 0, 0 ),
			get_abc_coeff( &vector, 3, 2, 1, 0, 1 ),
			get_abc_coeff( &vector, 3, 2, 1, 0, 2 ) 
			)) {
		printf("Coefficients don't match up for u_32-;1\n");
	} else {
		printf("Coefficients OK for u_32-;1\n");
	}
	if (!checkABC( -3, 1, 0, 
			get_abc_coeff( &vector, 3, 2, 1, 1, 0 ),
			get_abc_coeff( &vector, 3, 2, 1, 1, 1 ),
			get_abc_coeff( &vector, 3, 2, 1, 1, 2 ) 
			)) {
		printf("Coefficients don't match up for u_32-;2\n");
	} else {
		printf("Coefficients OK for u_32-;2\n");
	}
	if (!checkABC( 1, -1, 2/PI, 
			get_abc_coeff( &vector, 3, 2, 1, 2, 0 ),
			get_abc_coeff( &vector, 3, 2, 1, 2, 1 ),
			get_abc_coeff( &vector, 3, 2, 1, 2, 2 ) 
			)) {
		printf("Coefficients don't match up for u_32-;3\n");
	} else {
		printf("Coefficients OK for u_32-;3\n");
	}
	
	/* Output all of the abc values */
	for (j = 0; j < nr-1; j++) {
		for (k = 0; k < ntheta-1; k++) {
			for (pm = 0; pm < 2; pm++) {
				for (m = 0; m < 3; m++) {
					printf("u_%d%d^%d%s = %.5f + %.5fr + %.5ftheta\n",j,k,m,(pm?"-":"+"),
						get_abc_coeff( &vector, j, k, pm, m, 0 ),
						get_abc_coeff( &vector, j, k, pm, m, 1 ),
						get_abc_coeff( &vector, j, k, pm, m, 2 ));
				}
			}
		}
	}
	
	/* Do test runs of psi calculations */
	printf("\n\nCalculating Psi_00\n");
	get_psi_coefficients( &vector, 0, 0, psi );
	printf("psi_00 = %.5f %s %.5fr %s %.5ftheta\n",psi[0], psi[1]<0?"-":"+" ,fabs(psi[1]), psi[2]<0?"-":"+" ,fabs(psi[2]));
	
	printf("\nCalculating Psi_22\n");
	get_psi_coefficients( &vector, 2, 2, psi );
	printf("psi_22 = %.5f %s %.5fr %s %.5ftheta\n",psi[0], psi[1]<0?"-":"+" ,fabs(psi[1]), psi[2]<0?"-":"+" ,fabs(psi[2]));
	
	
	free_coefficient_storage( &vector );
	return 0;
}

/* check equality allowing for floating-point error */
bool isCloseEnough( double a, double b ) {
	return fabs( a - b ) <= 1e-6;
}

bool checkABC( double a1, double b1, double c1, double a2, double b2, double c2 ) {
	return isCloseEnough( a1, a2 ) && isCloseEnough( b1, b2 ) && isCloseEnough( c1, c2 );
}
#include "gauss-legendre.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265359

/* test function is -3x^5 + 5x^4 + 8x^3 - 6x^2 + 8x - 3 */
double testfunc1( double x, void *params ) {
	return -3*x*x*x*x*x + 5*x*x*x*x + 8*x*x*x - 6*x*x + 8*x - 3;
}

double testfunc2( double x, void *params ) {
	return x * exp( 2*x );
}

double testfunc3( double x, void *params ) {
	return sin( 5*x ) * cos( 2*x );
}

#define N_FUNCS 1

int main( int argc, char **argv ) {
	
	double *abscissae, *weights, solution;
	unsigned int order = 2;
	int i, converged;
	double accuracy = 1e-5, true_value;
	
	
	/* testfunc1: 6th order polynomial */
	converged = 0;
	order = 2;
	true_value = -13292.0;
		
	printf("%s\nOn interval [%f,%f]:\nTrue value = %.15f\n", "-3x^5 + 5x^4 + 8x^3 - 6x^2 + 8x - 3", 2.0, 6.0, true_value);
	gauss_legendre_eval_converge_ab( testfunc1, NULL, 2.0, 6.0, accuracy, &order );
	printf("Converged at N == %d\n\n",order);
	
	/* testfunc2: exponential */
	converged = 0;
	order = 2;
	true_value = 502.4387278411862;
		
	printf("%s\nOn interval [%f,%f]:\nTrue value = %.15f\n", "x * exp( 2x )", 1.0, 3.0, true_value);
	gauss_legendre_eval_converge_ab( testfunc2, NULL, 1.0, 3.0, accuracy, &order );
	printf("Converged at N == %d\n\n",order);
	
	/* testfunc3: trigonometric */
	converged = 0;
	order = 2;
	true_value = 10.0/21.0;
		
	printf("%s\nOn interval [%f,%f]:\nTrue value = %.15f\n", "sin( 5x )cos( 2x )", 0.0, 3*PI, true_value);
	gauss_legendre_eval_converge_ab( testfunc3, NULL, 0.0, 3.0*PI, accuracy, &order );
	printf("Converged at N == %d\n\n",order);
	
}
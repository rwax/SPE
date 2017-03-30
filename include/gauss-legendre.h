#ifndef __NCPA_GAUSS_LEGENDRE_H__
#define __NCPA_GAUSS_LEGENDRE_H__

#include <math.h>

/**
  * Typedef for a function pointer to be integrated.  Function should take as arguments:
  *      double (integration variable)
  *      void * (any other parameters that need to be passed into the function)
  */
typedef double (*integrand)(double, void *);

/**
  * Performs Gauss-Legendre quadrature for an arbitrary function over the standard interval
  * [-1,1]
  */
double gauss_legendre_eval( unsigned int order, integrand func, void *params, double *abscissae, 
	double *weights );
	
double gauss_legendre_eval_converge( integrand func, void *params, double accuracy, unsigned int *finalorder );
double gauss_legendre_eval_converge_ab( integrand func, void *params, double a, double b, 
	double accuracy, unsigned int *finalorder );
	
/**
  * Performs Gauss-Legendre quadrature for an arbitrary function over an arbitrary interval [a,b].
  * Works by creating a new set of weights and abscissae through change of variables:
  *
  * w'_i = 0.5 * w_i * (b-a)
  * a'_i = 0.5 * ((b-a)a_i + b + a)
  *
  * per https://math.okstate.edu/people/yqwang/teaching/math4513_fall14/Notes/gaussian.pdf
  */
double gauss_legendre_eval_ab( unsigned int order, integrand func, void *params, double *abscissae, 
	double *weights, double a, double b );
	
/**
  * Allocated and populates the tables of abscissae and weights for a given order.  Uses lookup tables
  * to do this.  Currently works up to order == 64.
  */
unsigned int gauss_legendre_table_alloc( unsigned int order, double **abscissae, double **weights );

/**
  * Frees the memory for the tables of weights and abscissae.
  */
void gauss_legendre_table_free( double *abscissae, double *weights );



#endif
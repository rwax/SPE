#ifndef __SPE_TRIANGULARIZATION_C__
#define __SPE_TRIANGULARIZATION_C__

//#include "gridfloat_dem.h"

/* How many Legendre polynomials and Bessel functions to iterate over */
#define L_MAX 20

typedef struct coeffVector {
	unsigned int N_r, N_theta, N_elements;
	double *coeffs, *r, *theta;
	double dr, dtheta;
} coeffVector;

/*******************************************************************************************/
/* This section is for calculating the triangularization of the ground surface.  Function  */
/* names refer to the variables in the notes (i.e. psi, script_G, etc).                    */
/*******************************************************************************************/

/**
  * void get_psi_coefficients
  * Calculates and adds up the three coefficients that specify the psi function for a ground surface triangularization.
  * psi should be a double[ 3 ] that has already been allocated.
  */
void get_psi_coefficients( const coeffVector *vector, unsigned int r_ind, unsigned int theta_ind, double *psi );

/**
  * void init_coefficient_storage
  * Initialize the storage vector for the coefficients.  Since there are 5 variables that
  * go into the calculation, and it's unnecessarily error-prone to make a 5-D array, we
  * make a 1-D vector and calculate the offsets that go in.
  *
  * We store:
  *     For each 1..N_r,
  *         For each 1..N_theta,
  *             For each +/-,
  *                 For each m in {1,2,3},
  *                     Store a, b, c as doubles
  *
  * Accepts N_r, N_theta, and a pointer to a coeffVector struct, which must have been declared already
  */
unsigned int init_coefficient_storage( unsigned int N_r, unsigned int N_theta, coeffVector *vector );
void free_coefficient_storage( coeffVector *vector );

/**
  * Given the maximum range r_max, and an initialized coeffVector, calculates all of the
  * {a,b,c} coefficients for the entire triangularization of the ground surface.
  */
void calculate_coefficients( double r_max, coeffVector *vector );

/**
  * Computes script G( x, y ) for x, y in spherical coordinates.  Calculate
  * the directional derivatives of the ground surface using a splined DEM.
  */
double script_G( double r_x, double theta_x, double phi_x, double r_y, double theta_y, double phi_y, 
		 double k, double R, double dhdy1, double dhdy2 );

/*******************************************************************************************/
/* Here follow a bunch of support functions.  In theory the user should not need to call   */
/* these.                                                                                  */
/*******************************************************************************************/


/** 
  * void calculate_interp_coefficients
  *
  * Calculate the a, b, and c coefficient vectors for a given triangle of
  *	points {r1,t1}, {r2,t2}, {r3,t3}
  * The coefficients are calculated by solving the equation
  *   | 1  r1  t1 |   | a_m |   | d_m1 |
  *   | 1  r2  t2 | x | b_m | = | d_m2 |
  *   | 1  r3  t3 |   | c_m |   | d_m3 |
  * where d_mn is the Kronecker delta function for m and n.  We solve using
  * LU decomposition.
  */
void calculate_interp_coefficients( const double *r, const double *theta, double *a, double *b, double *c );

/**
  * unsigned int indices2offset
  * Calculate the offset into the storage vector from the r, theta, +/-, m, and a/b/c indices
  */
unsigned int indices2offset( const coeffVector *vector, unsigned int r_ind, 
	unsigned int theta_ind, unsigned int plusminus, unsigned int m, unsigned int abc );

/**
  * Gets a single ABC coefficient from the coeffVector for a given r, theta, plusminus, m, and 
  * abc index.
  */	
double get_abc_coeff( const coeffVector *vector, unsigned int r_ind, 
	unsigned int theta_ind, unsigned int plusminus, unsigned int m, unsigned int abc );
/**
  * Sets a single ABC coefficient from the coeffVector for a given r, theta, plusminus, m, and 
  * abc index.
  */
void set_abc_coeff( coeffVector *vector, unsigned int r_ind, 
	unsigned int theta_ind, unsigned int plusminus, unsigned int m, unsigned int abc,
	double newVal );
	
/**
  * For a given r, theta, plusminus, and m index, gets the {a,b,c} vector making up the u function
  * for that triangle and vertex.
  */
void get_abc_vector( const coeffVector *vector, unsigned int r_ind, unsigned int theta_ind, 
		unsigned int plusminus, unsigned int m, double *abc );

/**
  * Adds two vectors of the same length together element-wise
  */
void add_vectors( int N, double *cumulative, const double *new_stuff );

/**
  * Computes the radial portion of the Green's function g_l(r,r')
  * Inputs: r1, r2 : r, r'
  *         R : Radius of hemispherical shell
  *         k : omega/c
  *         l : order of Bessel functions to use
  */
double greens_function_radial( double r1, double r2, double R, double k, double l );

/**
  * Computes the partial derivatives of \Omega wrt theta_y and phi_y
  */
double dOMEGAdtheta_y( double theta_x, double theta_y, double phi_x, double phi_y );
double dOMEGAdphi_y( double theta_x, double theta_y, double phi_x, double phi_y );


/*******************************************************************************************/
/* This section is for calculating the overpressure on the ground surface due to the       */
/* explosion.                                                                              */
/*******************************************************************************************/




#endif
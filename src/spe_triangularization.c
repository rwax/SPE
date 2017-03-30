#include <math.h>

#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_legendre.h"

#include "utilities.h"
#include "spe_triangularization.h"


/** 
  * void calculate_interp_coefficients
  * Calculate the a, b, and c coefficient vectors for a given triangle of
  *	points {r1,t1}, {r2,t2}, {r3,t3}
  * The coefficients are calculated by solving the equation
  *   | 1  r1  t1 |   | a_m |   | d_m1 |
  *   | 1  r2  t2 | x | b_m | = | d_m2 |
  *   | 1  r3  t3 |   | c_m |   | d_m3 |
  * where d_mn is the Kronecker delta function for m and n.  We solve using
  * LU decomposition.
  */
void calculate_interp_coefficients( const double *r, const double *theta, double *a, double *b, double *c ) {
	
	double **M, **LU;  /* matrix of r and t values */
	int i, j, m;
	double abc[ 3 ];
	double d[ 3 ];
	
	/* create the input matrix */
	M = real_matrix_alloc( 3, 3 );
	LU = real_matrix_alloc( 3, 3 );
	for (i = 0; i < 3; i++) {
		M[ i ][ 0 ] = 1;
		M[ i ][ 1 ] = r[ i ];
		M[ i ][ 2 ] = theta[ i ];
	}
	
#ifdef __NCPA_DEBUG__
	printf("Input r/theta matrix:\n");
	for (i = 0; i < 3; i++) {
		printf("| %.3f  %.3f  %.3f |\n",M[ i ][ 0 ],M[ i ][ 1 ],M[ i ][ 2 ]);
	}
#endif
	
	do_real_LU_decomp(3,LU,M);
	
	
	/* Solve this three times, once for each m index */
	for (m = 0; m < 3; m++) {
		
		/* Zero out the v1 and d vectors */
		for (j = 0; j < 3; j++) {
			abc[ j ] = 0;
			d[ j ] = 0;
		}
		
		/* set up delta function vector */
		d[ m ] = 1;
		
		/* solve */
		real_LU_linsolver( 3, LU, d, abc );
		
		/* distribute results */
		a[ m ] = abc[ 0 ];
		b[ m ] = abc[ 1 ];
		c[ m ] = abc[ 2 ];
#ifdef __NCPA_DEBUG__
		printf("u[%d] = %.4f + %.4fr + %.4ftheta\n",m, a[m], b[m], c[m]);
#endif
	}
	
	real_matrix_free( M, 3, 3 );
	real_matrix_free( LU, 3, 3 );
}

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
  */
unsigned int init_coefficient_storage( unsigned int N_r, unsigned int N_theta, coeffVector *vector ) {
	
	//double *M;
	//unsigned int vector_size;
	
	/* there are N_r-1 x N_theta-1 points associated with triangles */
	vector->N_r = N_r;
	vector->N_theta = N_theta;
	vector->N_elements = (N_r-1) * (N_theta-1) * 18;
	
#ifdef __NCPA_DEBUG__
	printf("Allocating %d bytes of memory to r vector\n",
		vector->N_r * (int)sizeof( double ) );
#endif
	vector->r = (double *)malloc( vector->N_r * sizeof( double ) );
	if (!(vector->r)) {
		printf("Error allocating r vector!\n");
		return 0;
	}
	
#ifdef __NCPA_DEBUG__
	printf("Allocating %d bytes of memory to theta vector\n",
		vector->N_theta * (int)sizeof( double ) );
#endif
	vector->theta = (double *)malloc( vector->N_theta * sizeof( double ) );
	if (!(vector->theta)) {
		printf("Error allocating theta vector!\n");
		return 0;
	}
	
	
#ifdef __NCPA_DEBUG__
	printf("Allocating %d bytes of memory to coefficient storage vector\n",
		vector->N_elements * (int)sizeof( double ));
#endif
	
	vector->coeffs = (double *)malloc( vector->N_elements * sizeof( double ) );
	if (!(vector->coeffs)) {
		printf("Error allocating coefficient storage vector!\n");
		return 0;
	}
	
	
	return 1;
}

void free_coefficient_storage( coeffVector *vector ) {
	free( vector->coeffs );
	free( vector->r );
	free( vector->theta );
}

/**
  * unsigned int indices2offset
  * Calculate the offset into the storage vector from the r, theta, +/-, m, and a/b/c indices
  */
unsigned int indices2offset( const coeffVector *vector, unsigned int r_ind, 
	unsigned int theta_ind, unsigned int plusminus, unsigned int m, unsigned int abc ) {
	
	return (r_ind * (vector->N_theta - 1) * 18) + (theta_ind * 18) + (plusminus * 9) + (m * 3) + abc;
}

double get_abc_coeff( const coeffVector *vector, unsigned int r_ind, 
	unsigned int theta_ind, unsigned int plusminus, unsigned int m, unsigned int abc ) {
		
	return vector->coeffs[ indices2offset( vector, r_ind, theta_ind, plusminus, m, abc ) ];		
}

void set_abc_coeff( coeffVector *vector, unsigned int r_ind, 
	unsigned int theta_ind, unsigned int plusminus, unsigned int m, unsigned int abc,
	double newVal ) {

#ifdef __NCPA_DEBUG__
	printf("Setting index %d to %.4f\n",
		indices2offset( vector, r_ind, theta_ind, plusminus, m, abc ),newVal);
#endif
		
	vector->coeffs[ indices2offset( vector, r_ind, theta_ind, plusminus, m, abc ) ] = newVal;
}

void calculate_coefficients( double r_max, coeffVector *vector ) {
	
	int j, k, m;
	double r123[ 3 ], t123[ 3 ], a[3], b[3], c[3];
	
	/* allocate and precalculate r and theta vectors */
	//r = (double *)malloc( vector->N_r * sizeof( double ) );
	//theta = (double *)malloc( vector->N_theta * sizeof( double ) );
	vector->dr = r_max / (double)(vector->N_r-1);
	vector->dtheta = 2 * PI / (double)(vector->N_theta-1);
	
	for (j = 0; j < vector->N_r; j++) {
		vector->r[ j ] = vector->dr * (double)j;
	}
	for (k = 0; k < vector->N_theta; k++) {
		vector->theta[ k ] = vector->dtheta * (double)k;
	}
	
	/* start chugging through */
	for (j = 0; j < (vector->N_r-1); j++) {
		for (k = 0; k < (vector->N_theta-1); k++) {
			
#ifdef __NCPA_DEBUG__
			printf("Working j=%d (r=%.4f), k=%d (theta=%.4f)\n", j, j*vector->dr, k, k*vector->dtheta);
#endif
			/* for + triangle, use 0=(j,k), 1=(j+1,k), 2=(j+1,k+1) */
			r123[ 0 ] = j * vector->dr;
			r123[ 1 ] = (j+1) * vector->dr;
			r123[ 2 ] = (j+1) * vector->dr;
			t123[ 0 ] = k * vector->dtheta;
			t123[ 1 ] = k * vector->dtheta;
			t123[ 2 ] = k == vector->N_theta-2 ? 0 : (k+1) * vector->dtheta;
			memset( a, 0, 3*sizeof(double) );
			memset( b, 0, 3*sizeof(double) );
			memset( c, 0, 3*sizeof(double) );
			calculate_interp_coefficients( r123, t123, a, b, c );
			for (m = 0; m < 3; m++) {

				set_abc_coeff( vector, j, k, 0, m, 0, a[ m ] );
				//coefficients[ indices2offset( N_theta, j, k, 0, m, 0 ) ] = a[ m ];
				set_abc_coeff( vector, j, k, 0, m, 1, b[ m ] );
				//coefficients[ indices2offset( N_theta, j, k, 0, m, 1 ) ] = b[ m ];
				set_abc_coeff( vector, j, k, 0, m, 2, c[ m ] );
				//coefficients[ indices2offset( N_theta, j, k, 0, m, 2 ) ] = c[ m ];
			}
			
			/* for - triangle, use 0=(j,k), 1=(j+1,k+1), 2=(j,k+1) */
			r123[ 0 ] = j * vector->dr;
			r123[ 1 ] = (j+1) * vector->dr;
			r123[ 2 ] = j * vector->dr;
			t123[ 0 ] = k * vector->dtheta;
			t123[ 1 ] = k == vector->N_theta-2 ? 0 : (k+1) * vector->dtheta;
			t123[ 2 ] = k == vector->N_theta-2 ? 0 : (k+1) * vector->dtheta;
			memset( a, 0, 3*sizeof(double) );
			memset( b, 0, 3*sizeof(double) );
			memset( c, 0, 3*sizeof(double) );
			calculate_interp_coefficients( r123, t123, a, b, c );
			for (m = 0; m < 3; m++) {
				set_abc_coeff( vector, j, k, 1, m, 0, a[ m ] );
				//coefficients[ indices2offset( N_theta, j, k, 1, m, 0 ) ] = a[ m ];
				set_abc_coeff( vector, j, k, 1, m, 1, b[ m ] );
				//coefficients[ indices2offset( N_theta, j, k, 1, m, 1 ) ] = b[ m ];
				set_abc_coeff( vector, j, k, 1, m, 2, c[ m ] );
				//coefficients[ indices2offset( N_theta, j, k, 1, m, 2 ) ] = c[ m ];
			}
		}
	}
}

void get_abc_vector( const coeffVector *vector, unsigned int r_ind, unsigned int theta_ind, 
	unsigned int plusminus, unsigned int m, double *abc ) {
	
	int i;
		
	memset(abc,0,3*sizeof(double));
	for (i = 0; i < 3; i++) {
		abc[ i ] = get_abc_coeff( vector, r_ind, theta_ind, plusminus, m, i );
	}
}

void get_psi_coefficients( const coeffVector *vector, unsigned int r_ind, unsigned int theta_ind, 
	double *psi ) {
	
	double abc[ 3 ];
	int i;
	
	/* 
	   exactly which abc coefficients get added together depends on the positions of the
	   r_index and theta_index relative to the edges of the calculation space
	*/
	memset(psi,0,3*sizeof(double));
	memset(abc,0,3*sizeof(double));
		
	/* Case 1: r_ind (j) == 0 */ 
	if (r_ind == 0) {
		if (theta_ind == 0 || theta_ind == vector->N_theta-1) {
			/* T00+: bottom left corner (index 0) */
			get_abc_vector( vector, 0, 0, 0, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",0,0,0,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* T00-: bottom left corner (index 0) */
			get_abc_vector( vector, 0, 0, 1, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",0,0,0,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* T0M-1: bottom right corner (index 2) */
			get_abc_vector( vector, 0, vector->N_theta-2, 1, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",0,vector->N_theta-2,2,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
		} else {
			/* T0k+: bottom left corner (index 0) */
			get_abc_vector( vector, 0, theta_ind, 0, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",0,theta_ind,0,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* T0k-: bottom left corner (index 0) */
			get_abc_vector( vector, 0, theta_ind, 1, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",0,theta_ind,0,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* T0k-1-: bottom right corner (index 2) */
			get_abc_vector( vector, 0, theta_ind, 1, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",0,theta_ind,2,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
		}
	/* Case 2: r_ind (j) == N_r-1 (N) */
	} else if (r_ind == vector->N_r-1) {
		if (theta_ind == 0 || theta_ind == vector->N_theta-1) {
			/* TN-1,0+: top left corner (index 1) */
			get_abc_vector( vector, r_ind-1, 0, 0, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,0,1,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* TN-1k-1-: top right corner (index 1) */
			get_abc_vector( vector, r_ind-1, vector->N_theta-2, 1, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,vector->N_theta-2,1,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* TN-1k-1+: top right corner (index 2) */
			get_abc_vector( vector, r_ind-1, vector->N_theta-2, 0, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1, vector->N_theta-2,2,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
		} else {
			/* TN-1,k+: top left corner (index 1) */
			get_abc_vector( vector, r_ind-1, theta_ind, 0, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,theta_ind,1,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* TN-1k-1-: top right corner (index 1) */
			get_abc_vector( vector, r_ind-1, theta_ind-1, 1, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,theta_ind-1,1,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
			
			/* TN-1k-1+: top right corner (index 2) */
			get_abc_vector( vector, r_ind-1, theta_ind-1, 0, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,theta_ind-1,2,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
			add_vectors( 3, psi, abc );
		}
	/* Case 3: r_ind in {1..N_r-2}, theta_ind == (0 || N_theta-1) */
	} else if (theta_ind == 0 || theta_ind == vector->N_theta-1) {
		/* Tj0+: index 0 */
		get_abc_vector( vector, r_ind, 0, 0, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind,0,0,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj0-: index 0 */
		get_abc_vector( vector, r_ind, 0, 1, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind,0,0,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj-1,0+: index 1 */
		get_abc_vector( vector, r_ind-1, 0, 0, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,0,1,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj-1M-1-: index 1 */
		get_abc_vector( vector, r_ind-1, vector->N_theta-2, 1, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,vector->N_theta-2,1,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj-1M-1+: index 2 */
		get_abc_vector( vector, r_ind-1, vector->N_theta-2, 0, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,vector->N_theta-2,2,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* TjM-1-: index 2 */
		get_abc_vector( vector, r_ind, vector->N_theta-2, 1, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind,vector->N_theta-2,2,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
	
	/* Case 4: Not at any edges */	
	} else {
		/* Tjk+: bottom left corner (index 0) */
		get_abc_vector( vector, r_ind, theta_ind, 0, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind,theta_ind,0,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tjk-: bottom left corner (index 0) */
		get_abc_vector( vector, r_ind, theta_ind, 1, 0, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind,theta_ind,0,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj-1k+: top left corner (index 1) */
		get_abc_vector( vector, r_ind-1, theta_ind, 0, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,theta_ind,1,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj-1k-1-: top right corner (index 1) */
		get_abc_vector( vector, r_ind-1, theta_ind-1, 1, 1, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,theta_ind-1,1,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tj-1k-1+: top right corner (index 2) */
		get_abc_vector( vector, r_ind-1, theta_ind-1, 0, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind-1,theta_ind-1,2,"+",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
		
		/* Tjk-1-: bottom right corner (index 2) */
		get_abc_vector( vector, r_ind, theta_ind-1, 1, 2, abc );
#ifdef __NCPA_DEBUG__
			printf("Adding u_%d%d^%d%s [%.5f,%.5f,%.5f] to [%.5f,%.5f,%.5f]\n",r_ind,theta_ind-1,2,"-",abc[0],abc[1],abc[2],psi[0],psi[1],psi[2]);
#endif
		add_vectors( 3, psi, abc );
	}
}

void add_vectors( int N, double *cumulative, const double *new_stuff ) {
	int i;
	for (i = 0; i < N; i++) {
		cumulative[ i ] += new_stuff[ i ];
	}
}

double script_G( double r_x, double theta_x, double phi_x, double r_y, double theta_y, double phi_y, 
	 	double k, double R, double dhdy1, double dhdy2 ) {
	
	/* 
	Need: 
		dh/dy1, dh/dy2
		k = w/c
	*/
	double y1, y2, OMEGA;
	double sin_theta_x, sin_theta_y, cos_theta_x, cos_theta_y,
			 sin_phi_x, sin_phi_y, cos_phi_x, cos_phi_y;
	double A, B, C;  /* coefficients for partial derivatives */
	double Pl_cosOMEGA, g_rx_ry, dgdr_y, dPdtheta_y, dPdphi_y, dPl_base;
	double summation;
	int l;
	
	sin_theta_x = sin( theta_x );
	sin_theta_y = sin( theta_y );
	cos_theta_x = cos( theta_x );
	cos_theta_y = cos( theta_y );
	sin_phi_x = sin( phi_x );
	sin_phi_y = sin( phi_y );
	cos_phi_x = cos( phi_x );
	cos_phi_y = cos( phi_y );
	
	OMEGA = acos( cos_theta_x*cos_theta_y + sin_theta_x*sin_theta_y*cos( phi_x - phi_y ) );
	
	/* compute normal derivative */
	A = cos_theta_y - dhdy1*sin_theta_y*cos_phi_y - dhdy2*sin_theta_y*sin_phi_y;
	B = -1/r_y * (sin_theta_y + dhdy1*cos_theta_y*cos_phi_y + dhdy2*cos_theta_y*sin_phi_y);
	C = 1/r_y * (dhdy1*sin_phi_y/r_y/sin_theta_y - dhdy2*cos_phi_y/r_y/sin_theta_y);
	
	summation = 0.0;
	for (l = 0; l <= L_MAX; l++) {
		
		Pl_cosOMEGA = gsl_sf_legendre_Pl( l, cos( OMEGA ) );
		
		/* d/dr_y[ Pl(cos(OMEGA))*gl(r_x,r_y) ] */
		if (r_x < r_y) {
			dgdr_y = k * k * gsl_sf_bessel_jl( l, k*r_x ) * (
				-gsl_sf_bessel_yl( l+1, k*r_y ) 
				+ (1/(k*r_y)) * gsl_sf_bessel_yl( l, k*r_y )
				- (gsl_sf_bessel_yl(l,k*R)/gsl_sf_bessel_jl(l,k*R))
					* (-gsl_sf_bessel_jl(l+1,k*r_y) + (1/(k*r_y))*gsl_sf_bessel_jl(l,k*r_y))
						);  /* @todo check R==r_y, f=j */
		} else if (r_x > r_y) {
			dgdr_y = k * k * 
				(-gsl_sf_bessel_jl(l+1,k*r_y) + 1/(k*r_y)*gsl_sf_bessel_jl(l,k*r_y)) *
				(gsl_sf_bessel_yl(l,k*r_x) - 
					gsl_sf_bessel_yl(l,k*R)/gsl_sf_bessel_jl(l,k*R)*gsl_sf_bessel_jl(l,k*r_x));
		} else {
			printf("Need to deal with case of r_x == r_y!\n");
		}
		
		/* d/dtheta_y[ Pl(cos(OMEGA))*gl(r_x,r_y) ] */
		g_rx_ry = greens_function_radial( r_x, r_y, R, k, l );
		dPl_base = 1/(sin(OMEGA)*sin(OMEGA)) * (
				-(l+1)*gsl_sf_legendre_Pl(l+1,cos(OMEGA)) +
				(l+1)*cos(OMEGA)*gsl_sf_legendre_Pl(l,cos(OMEGA))
			);
		dPdtheta_y = dPl_base * dOMEGAdtheta_y( theta_x, theta_y, phi_x, phi_y );
		
		
		/* d/dphi_y[ Pl(cos(OMEGA))*gl(r_x,r_y) ] */
		dPdphi_y = dPl_base * dOMEGAdphi_y( theta_x, theta_y, phi_x, phi_y );
		
		summation += (2*l+1)/(4*PI) * (
			A * Pl_cosOMEGA * dgdr_y +
			B * g_rx_ry * dPdtheta_y +
			C * g_rx_ry * dPdphi_y
			);
	}
	
	return summation;
	
}

double dOMEGAdtheta_y( double theta_x, double theta_y, double phi_x, double phi_y ) {
	
	double u, ddt;   /* convenience variable */
	
	u = cos(theta_x)*cos(theta_y) + sin(theta_x)*sin(theta_y)*cos(phi_x-phi_y);
	
	ddt = (-1/sqrt(1 - u*u)) * (
		sin(theta_x)*cos(theta_y)*cos(phi_x-phi_y) - cos(theta_x)*sin(theta_y)
			);
	return ddt;
}

double dOMEGAdphi_y( double theta_x, double theta_y, double phi_x, double phi_y ) {
	
	double u, ddp;
	u = cos(theta_x)*cos(theta_y) + sin(theta_x)*sin(theta_y)*cos(phi_x-phi_y);
	
	ddp = (-1/sqrt(1-u*u)) * (
		sin(theta_x)*sin(theta_y)*sin(phi_x-phi_y)
		);
	return ddp;
	
}

double greens_function_radial( double r1, double r2, double R, double k, double l ) {
	
	double r_gt, r_lt, g;
	
	if (r1 > r2) {
		r_gt = r1;
		r_lt = r2;
	} else {
		r_gt = r2;
		r_lt = r1;
	}
	
	g = k * gsl_sf_bessel_jl(l,k*r_lt) * (
			gsl_sf_bessel_yl(l,k*r_gt) 
			- gsl_sf_bessel_yl(l,k*R)/gsl_sf_bessel_jl(l,k*R)*gsl_sf_bessel_jl(l,k*r_gt)
			);
	return g;
}
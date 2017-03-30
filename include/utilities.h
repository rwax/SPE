#ifndef __RWAX_UTILITIES_H__
#define __RWAX_UTILITIES_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h> 

#define PI 3.14159265358979323846
#define d2r (M_PI/180.0)
#define ROOT_TOL 1.0e-10
#define FUNKT(x) ((*funct)(x))

/*         */
/* Filters */
/*         */
double hann( int begin, int center, int i );
double half_hann( int begin, int end, int i );
double cont_half_hann( double begin, double end, double x );

/*  Choose N as 2^power.  */
int make_N( int power );

/*               */
/* Data analysis */
/*               */
double correllation( double *vec1, double *vec2, double *corr,
					int n_fft, double f_s, double t_start,
					double t_end, double t_window, double f_b,
					double roll_b, double f_t, double roll_t );

/*               */
/* Root finders. */
/*               */
double complex cnewt_raph( double complex start_step, double complex guess,
						   double complex (*funct)(double complex) );

/*                              */
/* Memory allocation routines.  */
/*                              */
double complex **matrix_alloc( int n, int m );
void matrix_free( double complex **p, int n, int m );
double **real_matrix_alloc( int n, int m );
void real_matrix_free( double **p, int n, int m );

/*                                                         */
/* The routines below are used to read data in from files. */
/*                                                         */
int count_rows( char *filename, int num_col );
int count_rows_arbcol( char *filename );
int count_columns( char *filename );
int read_1col_datafile( char *filename, double *vec );
int read_real_data( char *filename, double *t_vec, double *vec );
int read_real_data_ncols( char *filename, int n_cols, double *t_vec, double **vec );
int read_real_data_ncols2( char *filename, int n_cols, double **vec );
int read_real_data_ncols3( char *filename, int n_cols, double **vec );
int read_data_vec( char *filename, double *q_vec, double complex *vec );
int read_data_vec_ncols( char *filename, int n_cols, double *q_vec, double complex **vec);
void read_header( FILE *f );

/*                                                      */
/* The routines below are used to read data into files. */
/*                                                      */
int save_points( char *filename, int n, int m, double **vec);
int save_complex_points( char *filename, int n, int m, double complex **vec);
  
/*                           */
/*  Linear algebra routines. */
/*                           */

/* Tridiagonal matrix solver. */
void tridiag( int, double complex *, double complex *,
			 double complex *, double complex *, double complex *);

/* LU routines */
void do_real_LU_decomp( int NNN, double **LU, double **matrix );
void real_LU_linsolver( int NNN, double **LU, double *source_vec, double *bc_vec);

/* 2 by 2 matrix utilities */
double complex det_2by2( double complex **MM );
void inverse_2by2( double complex **invMM, double complex **MM);
void eigs_2by2( double complex *eigs, double complex **eigvecs,
			   double complex **MM);

/*              */
/* ODE solvers. */
/*              */

/* Runge-Kutta solvers */
void rk_bc_wave_eq_ode( int N, double a, double b, double complex eps,
					   double (*funct)(double), double complex bc,
					   double complex **vec );

void reverse_rk_bc_wave_eq_ode( int N, double a, double b, double complex eps,
							   double (*funct)(double), double complex bc,
							   double complex **vec);


/*                          */
/* Interpolation algorithms */
/*                          */

/* Cubic spline interpolation routines. */
void cub_spline_coef1( int, double complex *, double *, double complex *);
double complex cub_spline( double complex *, double *, double complex *, double );
void cub_spline_coef2( int, double complex, double complex, double complex *,
					  double *, double complex * );
double complex d_cub_spline( double complex *, double *, double complex *, double );
double complex dd_cub_spline( double complex *, double *, double complex *, double );


/* Unwrap phase data. */
void unwrap( int, double * );

/* Haversine formulae for gps to distance */
double haversine_km_h( double lat1, double long1, double lat2, double long2 );
double haversine_km_v( double lat1, double long1, double lat2, double long2 );
double haversine_km( double lat1, double long1, double lat2, double long2 );

#endif
#ifndef __GRID_FLOAT_DEM_H__
#define __GRID_FLOAT_DEM_H__

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

typedef enum { BASIS_LATLON, BASIS_XY } map_basis;


/**
  * gridfloat_dem
  * Basically just a way of keeping the lat/lon grid associated with the DEM.  The
  * variable dem is indexed as dem[ lat ][ lon ] counting from the lower left corner
  * of the map area.  If the basis is set to USES_XY, indexing is therefore dem[ y ][ x ].
  */
typedef struct gridfloat_dem {
	unsigned int ny, nx;			/* map dimensions */
	double **dem, *yvec, *xvec;		/* the map itself */
	map_basis coordinate_basis;		/* Is this DEM on a lat/lon or x/y system? */
	double origin_lat, origin_lon;		/* only valid for BASIS_XY */
} gridfloat_dem;

/**
  * gridfloat_dem_spline
  * An object containing everything needed to do bicubic spline interpolation on
  * a DEM.  Note that the underlying DEM should not be changed once this object has been
  * created.
  */
typedef struct gridfloat_dem_spline {
	gsl_spline2d *spline;  			/* GSL spline object */
	gsl_interp_accel *xacc, *yacc;  	/* Accelerators */
	double *xa, *ya, *za;  			/* underlying x, y, and z vectors */
	size_t nx, ny;         			/* size of underlying grid */
	map_basis coordinate_basis;		/* BASIS_LATLON or BASIS_XY */
} gridfloat_dem_spline;


/**
  * Initializes the DEM from a file in the output format from the gridfloat command.
  */
int gridfloat_dem_init_from_file( gridfloat_dem *dem, const char *filename, double minlon, double maxlon, 
	double minlat, double maxlat );
	
/**
  * Free the allocated memory from the DEM struct
  */
void gridfloat_dem_free( gridfloat_dem *dem );

/**
  * Convert a lat/lon based DEM struct to X/Y based on a specified origin point
  */
void gridfloat_dem_convert_to_xy( gridfloat_dem *dem, double lat0, double lon0 );


int gridfloat_dem_spline_init( gridfloat_dem_spline *spline, gridfloat_dem *dem );
void gridfloat_dem_spline_free( gridfloat_dem_spline *spline );

/**
  * Evaluate the 2-D spline for a lat/lon point
  */
double gridfloat_dem_spline_eval_latlon( gridfloat_dem_spline *spline, double lat, double lon );
double gridfloat_dem_spline_eval_xy( gridfloat_dem_spline *spline, double x, double y );
double gridfloat_dem_spline_eval_rtheta( gridfloat_dem_spline *spline, double r, double theta );

/**
  * Get the X and Y derivatives at a particular point
  */
double gridfloat_dem_spline_eval_dx_latlon( gridfloat_dem_spline *spline, double lat, double lon );
double gridfloat_dem_spline_eval_dy_latlon( gridfloat_dem_spline *spline, double lat, double lon );
double gridfloat_dem_spline_eval_dx_xy( gridfloat_dem_spline *spline, double x, double y );
double gridfloat_dem_spline_eval_dy_xy( gridfloat_dem_spline *spline, double x, double y );
double gridfloat_dem_spline_eval_dx_rtheta( gridfloat_dem_spline *spline, double r, double theta );
double gridfloat_dem_spline_eval_dy_rtheta( gridfloat_dem_spline *spline, double r, double theta );

#endif
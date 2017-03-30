#include "gridfloat_dem.h"
#include "utilities.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/**
  * reads the header information for a USGS GridFloat DEM file and constructs
  * the metadata struct.
  */
int gridfloat_dem_init_from_file( gridfloat_dem *dem, const char *filename, double minlon, double maxlon, 
		double minlat, double maxlat ) {
		 	
	FILE *demfile;
	char buffer[ 64000 ];
	char *token;
	int nlon, nlat, latind, lonind;
	double dlat, dlon;
	
	/* open the file and count the columns and rows */
	demfile = fopen( filename, "r" );
	if (!demfile) {
		printf("Can't open %s!\n",filename);
		return 0;
	}
	memset( buffer, 0, 64000 );
	nlon = 0;
	nlat = 0;
	if (fgets( buffer, 64000, demfile ) != NULL) {
		
		/* chomp off the newline character, just to be safe */
		if (!feof( demfile )) {
			buffer[ strlen( buffer ) - 1] = 0;
		}
		
		/* Tokenize the string to count elements.  Assume map is rectangular */
		token = strtok( buffer, "[], " );
		while (token != NULL) {
			nlon++;
			token = strtok( NULL, "[], " );
		}
		
		nlat = 1;
	}
	
	/* read in the rest of the file */
	while (fgets( buffer, 64000, demfile ) != NULL) {
		nlat++;
	}
	
	
#ifdef __NCPA_DEBUG__
	printf("Counted %d longitude points, %d latitude points from %s\n", nlon, nlat, filename);
#endif
	
	/* initialize the vectors and map array */
	dem->ny = nlat;
	dem->nx = nlon;
	dem->yvec = (double *)malloc( nlat * sizeof( double ) );
	dem->xvec = (double *)malloc( nlon * sizeof( double ) );
	dem->dem = (double **)malloc( nlat * sizeof( double * ) );
	for (latind = 0; latind < nlat; latind++) {
		dem->dem[ latind ] = (double *)malloc( nlon * sizeof( double ) );
	}
	
	/* Calculate the lat/lon vectors */
	dlat = (minlat - maxlat) / (double)(nlat - 1);  /* negative dlat */
	dlon = (maxlon - minlon) / (double)(nlon - 1);
#ifdef __NCPA_DEBUG__
	printf("dlat = %.6f, dlon = %.6f\n", dlat, dlon);
#endif
	for (latind = 0; latind < nlat; latind++) {
		dem->yvec[ nlat - 1 - latind ] = maxlat + (double)latind*dlat;
	}
	for (lonind = 0; lonind < nlon; lonind++) {
		dem->xvec[ lonind ] = minlon + (double)lonind*dlon;
	}
#ifdef __NCPA_DEBUG__
	printf("Latitude Vector:\n");
	for (latind = 0; latind < nlat; latind++) {
		printf("%d: %.6f\n",latind,dem->yvec[ latind ]);
	}
	printf("\n\nLongitude Vector:\n");
	for (lonind = 0; lonind < nlon; lonind++) {
		printf("%d: %.6f\n",lonind,dem->xvec[ lonind ]);
	}
	printf("\n\n");
#endif
	
	/* read the file and put the elevation values in place */
	latind = 0;
	fseek(demfile,0,SEEK_SET);
	memset( buffer, 0, 64000 );
	while (fgets( buffer, 64000, demfile ) != NULL) {
		
		lonind = 0;
		
		/* chomp off the newline character, just to be safe */
		if (!feof( demfile )) {
			buffer[ strlen( buffer ) - 1] = 0;
		}
		
		token = strtok( buffer, "[], " );
		/* Remember to put lines in reverse order */
		while (token != NULL) {
			dem->dem[ nlat - 1 - latind ][ lonind ] = atof( token );
			lonind++;
			token = strtok( NULL, "[], " );
		}
		
		latind++;
		memset( buffer, 0, 64000 );
	}
	fclose( demfile );
	dem->coordinate_basis = BASIS_LATLON;
	
	return 1;
}

void gridfloat_dem_free( gridfloat_dem *dem ) {
	int i;
	
	free( dem->yvec );
	free( dem->xvec );
	for (i = 0; i < dem->ny; i++) {
		free( dem->dem[ i ] );
	}
	free( dem->dem );
}

void gridfloat_dem_convert_to_xy( gridfloat_dem *dem, double lat0, double lon0 ) {
	
	int i;
	
	if (dem->coordinate_basis == BASIS_XY) {
		printf("DEM is already in XY basis!\n");
		return;
	}
	
	/* convert lon vector to x coordinates */
	for (i = 0; i < dem->nx; i++) {
#ifdef __NCPA_DEBUG__
		printf("Lon %.4f --> %0.6f\n",dem->xvec[ i ], (dem->xvec[ i ] - lon0) * 111.320 * cos( lat0 * PI / 180.0 ));
#endif
		dem->xvec[ i ] = (dem->xvec[ i ] - lon0) * 111.320 * cos( lat0 * PI / 180.0 );

	}
	
	/* convert lat vector to y coordinates */
	for (i = 0; i < dem->ny; i++) {
#ifdef __NCPA_DEBUG__
		printf("Lat %.4f --> %0.6f\n",dem->yvec[ i ], (dem->yvec[ i ] - lat0) * 110.574);
#endif
		dem->yvec[ i ] = (dem->yvec[ i ] - lat0) * 110.574;
	}
	
	dem->coordinate_basis = BASIS_XY;
	dem->origin_lat = lat0;
	dem->origin_lon = lon0;
}



int gridfloat_dem_spline_init( gridfloat_dem_spline *spline, gridfloat_dem *dem ) {
	
	int lonind,latind;
	
#ifdef __NCPA_DEBUG__
	printf("Allocating vectors\n");
#endif
	spline->nx = dem->nx;
	spline->ny = dem->ny;
	spline->xa = (double *)malloc( dem->nx * sizeof( double ) );
	spline->ya = (double *)malloc( dem->ny * sizeof( double ) );
	spline->za = (double *)malloc( dem->ny * dem->nx * sizeof( double ) );
#ifdef __NCPA_DEBUG__
	printf("Done allocating vectors\n");
	printf("Initializing objects\n");
#endif
	
	spline->spline = gsl_spline2d_alloc( gsl_interp2d_bicubic, dem->nx, dem->ny );
	spline->xacc = gsl_interp_accel_alloc();
	spline->yacc = gsl_interp_accel_alloc();
#ifdef __NCPA_DEBUG__
	printf("Done initializing objects\n");
#endif
	
	memcpy( spline->xa, dem->xvec, dem->nx * sizeof( double ) );
	memcpy( spline->ya, dem->yvec, dem->ny * sizeof( double ) );
	for (latind = 0; latind < dem->ny; latind++) {
		for (lonind = 0; lonind < dem->nx; lonind++) {			
			gsl_spline2d_set(spline->spline, spline->za, 
				lonind, latind, 
				/*dem->xvec[ lonind ], dem->yvec[ latind ], */
				dem->dem[ latind ][ lonind ] );
		}
	}
	
	spline->coordinate_basis = dem->coordinate_basis;
	
#ifdef __NCPA_DEBUG__
	printf("Initializing final spline...");
#endif
      	/* initialize interpolation */
        gsl_spline2d_init(spline->spline, spline->xa, spline->ya, spline->za, spline->nx, spline->ny);
#ifdef __NCPA_DEBUG__
	printf("done!\n");
#endif
	return 1;
}

void gridfloat_dem_spline_free( gridfloat_dem_spline *spline ) {
	gsl_spline2d_free( spline->spline );
	gsl_interp_accel_free( spline->xacc );
	gsl_interp_accel_free( spline->yacc );
	free( spline->za );
	free( spline->xa );
	free( spline->ya );
}

double gridfloat_dem_spline_eval_latlon( gridfloat_dem_spline *spline, double lat, double lon ) {
	if (spline->coordinate_basis != BASIS_LATLON) {
		printf("Coordinate basis must be LATLON: returning 0.0\n");
		return 0.0;
	}
	
	return gsl_spline2d_eval( spline->spline, lon, lat, spline->xacc, spline->yacc );
}

double gridfloat_dem_spline_eval_xy( gridfloat_dem_spline *spline, double x, double y ) {
	if (spline->coordinate_basis != BASIS_XY) {
		printf("Coordinate basis must be XY: returning 0.0\n");
		return 0.0;
	}
	
	return gsl_spline2d_eval( spline->spline, x, y, spline->xacc, spline->yacc );
}

double gridfloat_dem_spline_eval_rtheta( gridfloat_dem_spline *spline, double r, double theta ) {
	
	double x, y;
	double trigtheta;
	
	if (spline->coordinate_basis != BASIS_XY ) {
		printf("Coordinate basis must be XY: returning 0.0\n");
		return 0.0;
	}
	

	
	trigtheta = 90-theta;
	x = r * cos( trigtheta * PI / 180.0);
	y = r * sin( trigtheta * PI / 180.0);
#ifdef __NCPA_DEBUG__
	printf("Evaluating (r,t) = (%.4f,%.4f) --> (x,y) = (%.4f,%.4f)\n",
		r, theta, x, y );
#endif
	
	return gridfloat_dem_spline_eval_xy( spline, x, y );
}

double gridfloat_dem_spline_eval_dx_latlon( gridfloat_dem_spline *spline, double lat, double lon ) {
	if (spline->coordinate_basis != BASIS_LATLON) {
		printf("Coordinate basis must be LATLON: returning 0.0\n");
		return 0.0;
	}
	
	return gsl_spline2d_eval_deriv_x( spline->spline, lon, lat, spline->xacc, spline->yacc );
}

double gridfloat_dem_spline_eval_dy_latlon( gridfloat_dem_spline *spline, double lat, double lon ) {
	if (spline->coordinate_basis != BASIS_LATLON) {
		printf("Coordinate basis must be LATLON: returning 0.0\n");
		return 0.0;
	}
	
	return gsl_spline2d_eval_deriv_y( spline->spline, lon, lat, spline->xacc, spline->yacc );
}

double gridfloat_dem_spline_eval_dx_xy( gridfloat_dem_spline *spline, double x, double y ) {
	if (spline->coordinate_basis != BASIS_XY) {
		printf("Coordinate basis must be XY: returning 0.0\n");
		return 0.0;
	}
	
	return gsl_spline2d_eval_deriv_x( spline->spline, x, y, spline->xacc, spline->yacc );
}

double gridfloat_dem_spline_eval_dy_xy( gridfloat_dem_spline *spline, double x, double y ) {
	if (spline->coordinate_basis != BASIS_XY) {
		printf("Coordinate basis must be XY: returning 0.0\n");
		return 0.0;
	}
	
	return gsl_spline2d_eval_deriv_y( spline->spline, x, y, spline->xacc, spline->yacc );
}

double gridfloat_dem_spline_eval_dx_rtheta( gridfloat_dem_spline *spline, double r, double theta ) {
	if (spline->coordinate_basis != BASIS_XY ) {
		printf("Coordinate basis must be XY: returning 0.0\n");
		return 0.0;
	}

	return gridfloat_dem_spline_eval_dx_xy( spline, r * cos( (90-theta) * PI / 180.0), 
		r * sin( (90-theta) * PI / 180.0) );
}

double gridfloat_dem_spline_eval_dy_rtheta( gridfloat_dem_spline *spline, double r, double theta ) {
	if (spline->coordinate_basis != BASIS_XY ) {
		printf("Coordinate basis must be XY: returning 0.0\n");
		return 0.0;
	}

	return gridfloat_dem_spline_eval_dy_xy( spline, r * cos( (90-theta) * PI / 180.0), 
		r * sin( (90-theta) * PI / 180.0) );
}
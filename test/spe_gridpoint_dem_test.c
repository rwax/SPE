#include "gridfloat_dem.h"
#include <stdio.h>
#include <math.h>

int findclosest( double *dvec, int veclength, double target ) {
	int i, t;
	double diff;
	
	diff = fabs( target - dvec[ 0 ] );
	t = 0;
	
	for (i = 1; i < veclength; i++) {
		if (fabs( target - dvec[ i ] ) < diff) {
			diff = fabs( target - dvec[ i ] );
			t = i;
		}
	}
	return t;
}


int main( int argc, char **argv ) {
	gridfloat_dem dem;
	gridfloat_dem_spline spline;
	double testlat, testlon, testval, testres;
	double minlon, maxlon, minlat, maxlat;
	FILE *LATLINE, *LONLINE, *AZLINE;
	int latind, lonind, i;
	int ncontrol = 8;
	double rcontrol[] = { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40 };
	double tcontrol[] = {    0,   45,   90,  135,  180,  225,  270,  315 };
	double zcontrol[] = { 1549, 1536, 1527, 1523, 1531, 1567, 1693, 1569 };
	
	
	minlon = -116.2;
	maxlon = -116.0;
	minlat = 37.3;
	maxlat = 37.4;
	
	gridfloat_dem_init_from_file( &dem, "output_canyon.elev", minlon, maxlon, minlat, maxlat );
	
	/* output some test points */
	printf("Last 3 points in file: %.3f  %.3f  %.3f\n",
		dem.dem[ 0 ][ dem.nx-3 ],
		dem.dem[ 0 ][ dem.nx-2 ],
		dem.dem[ 0 ][ dem.nx-1 ]);
		
	printf("First 3 points in file: %.3f  %.3f  %.3f\n\n",
		dem.dem[ dem.ny-1 ][ 0 ],
		dem.dem[ dem.ny-1 ][ 1 ],
		dem.dem[ dem.ny-1 ][ 2 ]);
		
	/* Initialize the spline */
#ifdef __NCPA_DEBUG__	
	printf("Setting up spline!\n");
#endif
	if (!gridfloat_dem_spline_init( &spline, &dem )) {
		printf("Error initializing spline!\n");
		return 1;
	}
#ifdef __NCPA_DEBUG__	
	printf("Spline done!\n");
#endif
	
	/* Extract a line for lat = 37.375 */
	testlat = 37.375;
	
	/* first the DEM points themselves */
	latind = findclosest( dem.yvec, dem.ny, testlat );
	LATLINE = fopen("test_lat_line_base.dat","w");
	if (!LATLINE) {
		printf("Error opening test_lat_line_base.dat\n");
		return 1;
	}
	for (lonind = 0; lonind < dem.nx; lonind++) {
		fprintf(LATLINE,"%.4f  %.4f\n", dem.xvec[ lonind ], dem.dem[ latind ][ lonind ]);
	}
	fclose(LATLINE);
	printf("Wrote base DEM points to test_lat_line_base.dat\n");
	
	testres = 0.0001;
	LATLINE = fopen("test_lat_line_spline.dat","w");
	if (!LATLINE) {
		printf("Error opening test_lat_line_spline.dat\n");
		return 1;
	}
	for (testlon = minlon; testlon <= maxlon; testlon += testres) {
		testval = gridfloat_dem_spline_eval_latlon( &spline, testlat, testlon );
#ifdef __NCPA_DEBUG__
		printf("(%.4f,%.4f) = %.4f\n",testlat,testlon,testval);
#endif
		fprintf(LATLINE, "%.4f  %.4f\n", testlon, testval );
	}
	fclose(LATLINE);
	printf("Wrote splined DEM points to test_lat_line_spline.dat\n");
	printf("Test in gnuplot with command:\nplot 'test_lat_line_control.dat', 'test_lat_line_base.dat', 'test_lat_line_spline.dat' with lines\n\n");
	
	
	/* Extract a line for lon = -116.066 */
	testlon = -116.066;
	lonind = findclosest( dem.xvec, dem.nx, -116.066 );
	LONLINE = fopen("test_lon_line_base.dat","w");
	if (!LONLINE) {
		printf("Error opening test_lon_line_base.dat\n");
		return 1;
	}
	for (latind = 0; latind < dem.ny; latind++) {
		fprintf(LONLINE,"%.4f  %.4f\n", dem.yvec[ latind ], dem.dem[ latind ][ lonind ]);
	}
	fclose(LONLINE);
	printf("Wrote base DEM points to test_lon_line_base.dat\n");
	
	LONLINE = fopen("test_lon_line_spline.dat","w");
	if (!LONLINE) {
		printf("Error opening test_lon_line_spline.dat\n");
		return 1;
	}
	for (testlat = minlat; testlat <= maxlat; testlat += testres) {
		testval = gridfloat_dem_spline_eval_latlon( &spline, testlat, testlon );
#ifdef __NCPA_DEBUG__
		printf("(%.4f,%.4f) = %.4f\n",testlat,testlon,testval);
#endif
		fprintf(LONLINE, "%.4f  %.4f\n", testlat, testval );
	}
	fclose(LONLINE);
	printf("Wrote splined DEM points to test_lon_line_spline.dat\n");
	printf("Test in gnuplot with command:\nplot 'test_lon_line_control.dat', 'test_lon_line_base.dat', 'test_lon_line_spline.dat' with lines\n\n");
		
		
	/* Read new DEM and convert to XY to test out r,theta calculation */
	gridfloat_dem_spline_free( &spline );
	gridfloat_dem_free( &dem );
	gridfloat_dem_init_from_file( &dem, "output_gz.elev", -116.2, -116.0, 37.0, 37.5 );
	gridfloat_dem_convert_to_xy( &dem, 37.2235, -116.0614 );
	if (!gridfloat_dem_spline_init( &spline, &dem )) {
		printf("Error initializing spline!\n");
		return 1;
	}
	
	/* Step through precalculated control values */
#ifdef __NCPA_DEBUG__
	printf("Starting (r,theta) test...\n");
#endif
	for (i = 0; i < ncontrol; i++) {
#ifdef __NCPA_DEBUG__
		printf("Testing (r,theta) = (%.4f,%.4f)\n",rcontrol[ i ], tcontrol[ i ] );
#endif
		testval = gridfloat_dem_spline_eval_rtheta( &spline, rcontrol[ i ], tcontrol[ i ] );
		printf("r = %0.2f, theta = %0.2f, z_0 = %0.2f, z_calc = %0.2f\n",
			rcontrol[ i ], tcontrol[ i ], zcontrol[ i ], testval );
	}
		
	gridfloat_dem_spline_free( &spline );
	gridfloat_dem_free( &dem );
	
	return 0;
}


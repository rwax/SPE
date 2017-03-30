CC=gcc
#SPE_LIBS=src/spe_guts.c src/utilities.c
CFLAGS=-Iinclude -I/usr/local/include -L/usr/local/lib -lfftw3 -lgsl -lm
#DEBUG=-D__NCPA_DEBUG__ -ggdb
DEBUG=


test: coeffs-test gridpoint-dem-test gauss-legendre-test


coeffs-test: test/spe_coeffs_test.c src/spe_triangularization.c src/utilities.c
	$(CC) -o $@ $^ $(CFLAGS) $(DEBUG)
	
gridpoint-dem-test: test/spe_gridpoint_dem_test.c src/gridfloat_dem.c
	$(CC) -o $@ $^ $(CFLAGS) $(DEBUG)
	
gauss-legendre-test: test/gauss-legendre-test.c src/gauss-legendre.c
	$(CC) -o $@ $^ $(CFLAGS) $(DEBUG)

surface-test: test/ground-surface-test.c src/spe_triangularization.c src/spe_groundsurface.c
	$(CC) -o $@ $^ $(CFLAGS) $(DEBUG)
	
integral-test: test/integral-test.c src/utilities.c src/gauss-legendre.c
	$(CC) -o $@ $^ $(CFLAGS) $(DEBUG)
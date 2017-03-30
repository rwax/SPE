unset key

set xlabel "Eastings [m]" offset 0,-1
set ylabel "Westings [m]" offset 1
set zlabel "Ground Altitude [m]" rotate by 90

set ztics 40
set term post enh eps mono solid 20
set out "SPE_ground_surface.eps"

splot "~/work/SPE/ground_model/NV_DEM.dat" using 2:3:1 with points

! epstopdf SPE_ground_surface.eps;rm SPE_ground_surface.eps


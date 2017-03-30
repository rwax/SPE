set autoscale
unset key

RR=50.0
x0=5.0
y0=6.0

set parametric
set trange [0:1]

set xlabel "Azimuth [rad]" 
set ylabel "Elevation [rad]" offset 1
set yrange [-0.0*pi:0.5*pi]

set size 1.8,1
set term post eps enh mono solid 22 lw 4
set out "sphere_section_elements.eps"

set multiplot 
set size 1,1
set origin 0,0
set title "Elements in the {/Symbol f}-{/Symbol q} Plane"
plot [0:0.5*pi] [0:2*pi] [0.5*pi:0] "theta_phi_grid.dat" lt 1,\
     "theta_phi_grid1.dat" lt 1,\
     2*pi, t lt 1

set origin 0.86,0
set title "Vertices on the Spherical Section"
set xlabel ""
set ylabel ""

set autoscale

unset border
unset tics
unset pm3d 
set surface
#splot "sphere_section_grid.dat" using 1:2:3
#replot "sphere_section_grid.dat" using 1:2:3 with dots lt 1 lw 3
splot "sphere_section_grid.dat" using 1:2:3 with dots lt 1 lw 3
unset multiplot

! epstopdf sphere_section_elements.eps;rm sphere_section_elements.eps
unset key
unset polar
set parametric

h_r=0.1
az_steps=15
h_az=2*pi/az_steps

line_x(xj,yj,xk,yk,x)=xj+(xk-xj)*x
line_y(xj,yj,xk,yk,x)=yj+(yk-yj)*x

set size 1.8,1
set term post eps enh mono solid 22 lw 4
set out "ground_elements.eps"

set multiplot 
set size 1,1
set origin 0,0
set xrange [0:2*pi]
set xlabel "Azimuth [rad]"
set ylabel "Scaled radius" offset 1
set title "Elements in the R-{/Symbol q} Plane"
plot [0:1] \
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(1-1)*h_r,ii*h_az,1*h_r,t),line_y((ii-1)*h_az,(1-1)*h_r,ii*h_az,1*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(1-1)*h_r,ii*h_az,1*h_r,t),line_y((ii-1)*h_az,(1-1)*h_r,ii*h_az,1*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(1-1)*h_r,(ii-1)*h_az,1*h_r,t),line_y((ii-1)*h_az,(1-1)*h_r,(ii-1)*h_az,1*h_r,t) lt 1,\
line_x((ii-1)*h_az,(1-1)*h_r,(ii-1)*h_az,1*h_r,t),line_y((ii-1)*h_az,(1-1)*h_r,(ii-1)*h_az,1*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,1*h_r,ii*h_az,1*h_r,t),line_y((ii-1)*h_az,1*h_r,ii*h_az,1*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(2-1)*h_r,ii*h_az,2*h_r,t),line_y((ii-1)*h_az,(2-1)*h_r,ii*h_az,2*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(2-1)*h_r,ii*h_az,2*h_r,t),line_y((ii-1)*h_az,(2-1)*h_r,ii*h_az,2*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(2-1)*h_r,(ii-1)*h_az,2*h_r,t),line_y((ii-1)*h_az,(2-1)*h_r,(ii-1)*h_az,2*h_r,t) lt 1,\
line_x((ii-1)*h_az,(2-1)*h_r,(ii-1)*h_az,2*h_r,t),line_y((ii-1)*h_az,(2-1)*h_r,(ii-1)*h_az,2*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,2*h_r,ii*h_az,2*h_r,t),line_y((ii-1)*h_az,2*h_r,ii*h_az,2*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(3-1)*h_r,ii*h_az,3*h_r,t),line_y((ii-1)*h_az,(3-1)*h_r,ii*h_az,3*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(3-1)*h_r,ii*h_az,3*h_r,t),line_y((ii-1)*h_az,(3-1)*h_r,ii*h_az,3*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(3-1)*h_r,(ii-1)*h_az,3*h_r,t),line_y((ii-1)*h_az,(3-1)*h_r,(ii-1)*h_az,3*h_r,t) lt 1,\
line_x((ii-1)*h_az,(3-1)*h_r,(ii-1)*h_az,3*h_r,t),line_y((ii-1)*h_az,(3-1)*h_r,(ii-1)*h_az,3*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,3*h_r,ii*h_az,3*h_r,t),line_y((ii-1)*h_az,3*h_r,ii*h_az,3*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(4-1)*h_r,ii*h_az,4*h_r,t),line_y((ii-1)*h_az,(4-1)*h_r,ii*h_az,4*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(4-1)*h_r,ii*h_az,4*h_r,t),line_y((ii-1)*h_az,(4-1)*h_r,ii*h_az,4*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(4-1)*h_r,(ii-1)*h_az,4*h_r,t),line_y((ii-1)*h_az,(4-1)*h_r,(ii-1)*h_az,4*h_r,t) lt 1,\
line_x((ii-1)*h_az,(4-1)*h_r,(ii-1)*h_az,4*h_r,t),line_y((ii-1)*h_az,(4-1)*h_r,(ii-1)*h_az,4*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,4*h_r,ii*h_az,4*h_r,t),line_y((ii-1)*h_az,4*h_r,ii*h_az,4*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(5-1)*h_r,ii*h_az,5*h_r,t),line_y((ii-1)*h_az,(5-1)*h_r,ii*h_az,5*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(5-1)*h_r,ii*h_az,5*h_r,t),line_y((ii-1)*h_az,(5-1)*h_r,ii*h_az,5*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(5-1)*h_r,(ii-1)*h_az,5*h_r,t),line_y((ii-1)*h_az,(5-1)*h_r,(ii-1)*h_az,5*h_r,t) lt 1,\
line_x((ii-1)*h_az,(5-1)*h_r,(ii-1)*h_az,5*h_r,t),line_y((ii-1)*h_az,(5-1)*h_r,(ii-1)*h_az,5*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,5*h_r,ii*h_az,5*h_r,t),line_y((ii-1)*h_az,5*h_r,ii*h_az,5*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(6-1)*h_r,ii*h_az,6*h_r,t),line_y((ii-1)*h_az,(6-1)*h_r,ii*h_az,6*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(6-1)*h_r,ii*h_az,6*h_r,t),line_y((ii-1)*h_az,(6-1)*h_r,ii*h_az,6*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(6-1)*h_r,(ii-1)*h_az,6*h_r,t),line_y((ii-1)*h_az,(6-1)*h_r,(ii-1)*h_az,6*h_r,t) lt 1,\
line_x((ii-1)*h_az,(6-1)*h_r,(ii-1)*h_az,6*h_r,t),line_y((ii-1)*h_az,(6-1)*h_r,(ii-1)*h_az,6*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,6*h_r,ii*h_az,6*h_r,t),line_y((ii-1)*h_az,6*h_r,ii*h_az,6*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(7-1)*h_r,ii*h_az,7*h_r,t),line_y((ii-1)*h_az,(7-1)*h_r,ii*h_az,7*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(7-1)*h_r,ii*h_az,7*h_r,t),line_y((ii-1)*h_az,(7-1)*h_r,ii*h_az,7*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(7-1)*h_r,(ii-1)*h_az,7*h_r,t),line_y((ii-1)*h_az,(7-1)*h_r,(ii-1)*h_az,7*h_r,t) lt 1,\
line_x((ii-1)*h_az,(7-1)*h_r,(ii-1)*h_az,7*h_r,t),line_y((ii-1)*h_az,(7-1)*h_r,(ii-1)*h_az,7*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,7*h_r,ii*h_az,7*h_r,t),line_y((ii-1)*h_az,7*h_r,ii*h_az,7*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(8-1)*h_r,ii*h_az,8*h_r,t),line_y((ii-1)*h_az,(8-1)*h_r,ii*h_az,8*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(8-1)*h_r,ii*h_az,8*h_r,t),line_y((ii-1)*h_az,(8-1)*h_r,ii*h_az,8*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(8-1)*h_r,(ii-1)*h_az,8*h_r,t),line_y((ii-1)*h_az,(8-1)*h_r,(ii-1)*h_az,8*h_r,t) lt 1,\
line_x((ii-1)*h_az,(8-1)*h_r,(ii-1)*h_az,8*h_r,t),line_y((ii-1)*h_az,(8-1)*h_r,(ii-1)*h_az,8*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,8*h_r,ii*h_az,8*h_r,t),line_y((ii-1)*h_az,8*h_r,ii*h_az,8*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(9-1)*h_r,ii*h_az,9*h_r,t),line_y((ii-1)*h_az,(9-1)*h_r,ii*h_az,9*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(9-1)*h_r,ii*h_az,9*h_r,t),line_y((ii-1)*h_az,(9-1)*h_r,ii*h_az,9*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(9-1)*h_r,(ii-1)*h_az,9*h_r,t),line_y((ii-1)*h_az,(9-1)*h_r,(ii-1)*h_az,9*h_r,t) lt 1,\
line_x((ii-1)*h_az,(9-1)*h_r,(ii-1)*h_az,9*h_r,t),line_y((ii-1)*h_az,(9-1)*h_r,(ii-1)*h_az,9*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,9*h_r,ii*h_az,9*h_r,t),line_y((ii-1)*h_az,9*h_r,ii*h_az,9*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(10-1)*h_r,ii*h_az,10*h_r,t),line_y((ii-1)*h_az,(10-1)*h_r,ii*h_az,10*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(10-1)*h_r,ii*h_az,10*h_r,t),line_y((ii-1)*h_az,(10-1)*h_r,ii*h_az,10*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,(10-1)*h_r,(ii-1)*h_az,10*h_r,t),line_y((ii-1)*h_az,(10-1)*h_r,(ii-1)*h_az,10*h_r,t) lt 1,\
line_x((ii-1)*h_az,(10-1)*h_r,(ii-1)*h_az,10*h_r,t),line_y((ii-1)*h_az,(10-1)*h_r,(ii-1)*h_az,10*h_r,t) lt 1,\
for [ii=1:az_steps:1] line_x((ii-1)*h_az,10*h_r,ii*h_az,10*h_r,t),line_y((ii-1)*h_az,10*h_r,ii*h_az,10*h_r,t) lt 1

set polar
set size square 1,1
set origin 0.88,0

set autoscale
unset tics 
unset border
set xlabel ""
set ylabel ""
set title "Mapping onto the Ground Projection"
replot

unset multiplot

! epstopdf ground_elements.eps;rm ground_elements.eps

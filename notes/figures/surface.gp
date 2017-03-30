unset grid
unset key

set isosamples 40
unset tics
unset border

f(x)= x<0 ? 0 : 2.0*(-0.5+1.0/(1.0+exp(-0.2*x*x)))
g(x)= x<0 ? 0 : 4.0*(-0.5+1.0/(1.0+exp(-0.1*x*x)))
h(x,y)=f(x)*g(y)

sph(x,y)=sqrt(50.0-(x-5)**2-(y-6)**2)

set size 1,0.6
set term post enh eps mono solid
set out "surface.eps"
splot [-6:16] [-6:16] [0:8] h(x,y),sph(x,y) lt 1

!epstopdf surface.eps;rm surface.eps

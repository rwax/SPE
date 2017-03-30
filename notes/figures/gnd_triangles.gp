unset key
unset polar
set parametric

unset tics
unset border

line_x(xj,yj,xk,yk,x)=xj+(xk-xj)*x
line_y(xj,yj,xk,yk,x)=yj+(yk-yj)*x

set label "x_{jk}" at -0.05,-0.05
set label "x_{jk+1}" at 1.02,-0.05
set label "x_{j+1k+1}" at 1.02,1.05
set label "S@_{jk}^{jk+1}" at 0.5,-0.07
set label "S@_{jk+1}^{j+1k}" at 1.02,0.5
set label "S@_{jk}^{j+1k+1}" at 0.5,0.66
set label "T@_{jk}^{&.{/Symbol -}}" at 0.6,0.3

set size 1.4,1
set term post eps enh mono solid 28 lw 6
set out "ground_triangles.eps"

set multiplot
set size 0.8,0.8
set origin 0.5,0.1

plot [0:1] line_x(0,0,1,0,t),line_y(0,0,1,0,t) lt 1,\
     line_x(0,0,1,1,t),line_y(0,0,1,1,t) lt 1,\
     line_x(1,0,1,1,t),line_y(1,0,1,1,t) lt 1

unset label

set label "x_{jk}" at -0.05,-0.05
set label "x_{j+1k}" at -0.07,1.05
set label "x_{j+1k+1}" at 1.02,1.04
set label "S@_{jk+1}^{j+1k+1}" at 0.5,1.06
set label "S@_{jk}^{jk+1}" at -0.11,0.5
set label "S@_{jk}^{j+1k+1}" at 0.5,0.45
set label "T@_{jk}^{&.+}" at 0.35,0.75

set origin 0.1,0.1

plot [0:1] line_x(0,0,0,1,t),line_y(0,0,0,1,t) lt 1,\
     line_x(0,0,1,1,t),line_y(0,0,1,1,t) lt 1,\
     line_x(0,1,1,1,t),line_y(0,1,1,1,t) lt 1

unset multiplot

! epstopdf ground_triangles.eps;rm ground_triangles.eps

unset label
set origin 0,0
set size 1,1

set label "(r_1,{/Symbol q}_1)" at -0.1,-0.1
set label "(r_2,{/Symbol q}_2)" at -0.62,1.1
set label "(r_3,{/Symbol q}_3)" at 0.95,0.9

set out "gen_triangle.eps"
plot [0:1] [-0.6:1.1] [-0.2:1.2] line_x(0,0,-0.5,1,t),line_y(0,0,-0.5,1,t) lt 1,line_x(0,0,1,0.8,t),line_y(0,0,1,0.8,t) lt 1,line_x(-0.5,1,1,0.8,t),line_y(-0.5,1,1,0.8,t) lt 1

! epstopdf gen_triangle.eps;rm gen_triangle.eps


unset label
unset tics
unset border
set size 1,1
set term post eps enh mono solid 26 lw 6
set out "ground_triangles_interp.eps"

unset label

set label "x_{jk}" at 0.05,-0.1
set label "x_{j+1k}" at -0.07,1.15
set label "x_{j+1k+1}" at 1.02,1.04
set label "x_{j-1k}" at -0.05,-1.1
set label "x_{jk-1}" at -1.2,-0.05
set label "x_{j-1k-1}" at -1.05,-1.1
set label "x_{jk+1}" at 1.05,-0.05
set label "T@_{jk}^{&.+}" at 0.35,0.75
set label "T@_{jk}^{&.-}" at 0.75,0.35
set label "T@_{j-1k}^{&.+}" at 0.3,-0.35
set label "T@_{j-1k-1}^{&.-}" at -0.35,-0.75
set label "T@_{j-1k-1}^{&.+}" at -0.75,-0.35
set label "T@_{jk-1}^{&.-}" at -0.35,0.35
set label "{/ZapfDingbats l}" at -0.038,0

plot [0:1] [-1.2:1.2] [-1.2:1.2] line_x(0,0,1,0,t),line_y(0,0,1,0,t) lt 1,\
     line_x(0,0,1,1,t),line_y(0,0,1,1,t) lt 1,\
     line_x(1,0,1,1,t),line_y(1,0,1,1,t) lt 1,\
     line_x(0,0,0,1,t),line_y(0,0,0,1,t) lt 1,\
     line_x(0,1,1,1,t),line_y(0,1,1,1,t) lt 1,\
     line_x(0,-1,0,0,t),line_y(0,-1,0,0,t) lt 1,\
     line_x(0,-1,1,0,t),line_y(0,-1,1,0,t) lt 1,\
     line_x(-1,0,0,0,t),line_y(-1,0,0,0,t) lt 1,\
     line_x(-1,0,0,1,t),line_y(-1,0,0,1,t) lt 1,\
     line_x(-1,-1,0,0,t),line_y(-1,-1,0,0,t) lt 1,\
     line_x(-1,-1,0,-1,t),line_y(-1,-1,0,-1,t) lt 1,\
     line_x(-1,-1,-1,0,t),line_y(-1,-1,-1,0,t) lt 1

! epstopdf ground_triangles_interp.eps;rm ground_triangles_interp.eps

unset label
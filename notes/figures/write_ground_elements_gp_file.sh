#!/bin/bash


for jj in 1 2 3 4 5 6 7 8 9 10
do
    echo "for [ii=1:10:1] line_x((ii-1)*h_az,($jj-1)*h_r,ii*h_az,$jj*h_r,t),line_y((ii-1)*h_az,($jj-1)*h_r,ii*h_az,$jj*h_r,t) lt 1,\\"
    echo "for [ii=1:10:1] line_x((ii-1)*h_az,($jj-1)*h_r,ii*h_az,$jj*h_r,t),line_y((ii-1)*h_az,($jj-1)*h_r,ii*h_az,$jj*h_r,t) lt 1,\\"
    echo "for [ii=1:10:1] line_x((ii-1)*h_az,($jj-1)*h_r,(ii-1)*h_az,$jj*h_r,t),line_y((ii-1)*h_az,($jj-1)*h_r,(ii-1)*h_az,$jj*h_r,t) lt 1,\\"
    echo "line_x((ii-1)*h_az,($jj-1)*h_r,(ii-1)*h_az,$jj*h_r,t),line_y((ii-1)*h_az,($jj-1)*h_r,(ii-1)*h_az,$jj*h_r,t) lt 1,\\"
    echo "for [ii=1:10:1] line_x((ii-1)*h_az,$jj*h_r,ii*h_az,$jj*h_r,t),line_y((ii-1)*h_az,$jj*h_r,ii*h_az,$jj*h_r,t) lt 1,\\"
done > temp.txt

cat gp_header.txt temp.txt gp_ender.txt > test.gp

rm temp.txt


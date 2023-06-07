reset
set terminal svg size 1600,1200 enhanced font "Helvetica,14"

set size ratio -1
set xrange [0:10]
set yrange [0:50]
set key off
set palette rgb 33,13,10; set title "y-Displ."
set colorbox vertical user origin .8, .1 size .04,.8
set output "../plot/plot2.svg"
plot '../data/data2.dat' u 1:2:3 w image, \
     '../data/data.dat' with lines linestyle 1





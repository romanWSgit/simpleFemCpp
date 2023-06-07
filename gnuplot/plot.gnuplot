reset
 set terminal svg size 1600,1200 enhanced font "Helvetica,14"

set xlabel "x"
set ylabel "f(x)"
set grid
set key outside right center
set size ratio -1
set output "../plot/plot.svg"
plot "../data/data.dat" title "Mesh"  with lines linetype 1 \
          linecolor rgb "#0000FF"  linewidth 3 ,\
     "../data/data.dat" using ($1+1*$3):($2+1*$4) \
          title "Mesh+U"  with lines linetype 1 \
          linecolor rgb "#FF0000"  linewidth 3





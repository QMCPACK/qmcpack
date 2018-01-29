#/usr/bin/gnuplot

set term png size 640,480
set output "compare.png"
set key right center

set xlabel "Number of water molecules"
set ylabel "Energy per molecule"
p 'dmc_no.dat' u 1:($2/$1):($4/$1) w errorlines t "LA", \
  'dmc_v0.dat' u 1:($2/$1):($4/$1) w errorlines t "TM0", \
  'dmc_v1.dat' u 1:($2/$1):($4/$1) w errorlines t "TM1"

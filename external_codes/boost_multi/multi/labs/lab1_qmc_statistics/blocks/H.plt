#set terminal png size 400,300 
#set output 'nblock_vs_tcpu_energy_H_STO-2G.png'

TOP = 0.98
DY  = 0.42

unset key

set multiplot 

set xrange [200:200000]
set logscale x

#Total energy
set xlabel "Number of blocks"
set yrange [-0.458:-0.447]
set tmargin at screen TOP-DY
set bmargin at screen TOP-2*DY
#set ytics  -0.455,0.001,-0.450
set ylabel "E_total\n(hartree)"
plot "H.dat" using 1:5:7 with errorbars 

#Total execution time
set xtics format ''
unset xlabel
set logscale y
set yrange [5:5000]
set tmargin at screen TOP
set bmargin at screen TOP-DY
#set ytics  0,50,300
set ylabel "t_CPU,total (s)" offset -3
plot "H.dat" using 1:12

unset multiplot
pause -1

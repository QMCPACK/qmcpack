#set terminal png size 400,300 
#set output 'timestep_vs_autocorrelation_energy_H_STO-2G.png'

TOP = 0.98
DY  = 0.42

unset key

set multiplot 

set xrange [0.00005:20]
set logscale x

#Total energy
set xlabel "Time step (1/hartree)"
set yrange [-0.495:-0.435]
set tmargin at screen TOP-DY
set bmargin at screen TOP-2*DY
set ytics  -0.49,0.01,-0.44
set ylabel "E_total\n(hartree)"
plot "H.dat" using 1:2:4 with errorbars 

#Autocorrelation
set xtics format ''
unset xlabel
set yrange [0:15]
set tmargin at screen TOP
set bmargin at screen TOP-DY
set ytics  1,2,14
set ylabel "Autocorrelation\n(1/hartree)" offset -2
plot "H.dat" using 1:5 

unset multiplot
pause -1

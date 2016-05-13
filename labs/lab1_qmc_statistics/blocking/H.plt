#set terminal png size 400,300 
#set output 'steps_per_block_vs_autocorrelation_energy_H_STO-2G.png'

TOP = 0.98
DY  = 0.42

unset key

set multiplot 

set xrange [0.5:8000]
set logscale x

#Total energy
set xlabel "Steps per block"
set yrange [-0.525:-0.375]
set tmargin at screen TOP-DY
set bmargin at screen TOP-2*DY
set ytics  -0.52,0.02,-0.38
set ylabel "E_total\n(hartree)"
plot "H.dat" using 1:7:9 with errorbars 

#Autocorrelation
set xtics format ''
unset xlabel
set yrange [0:200]
set tmargin at screen TOP
set bmargin at screen TOP-DY
set ytics  10,20,190
set ylabel "Autocorrelation\n(1/hartree)" offset -1
plot "H.dat" using 1:10

unset multiplot
pause -1

#set terminal png size 400,300 
#set output 'timestep_vs_autocorrelation_energy_H_STO-2G.png'

#set xrange [0.00005:20]
#set logscale x

f(x) = a * x**2 + b * x + c
fit f(x) "nelectron_tcpu.dat" using 1:2 via a,b,c

unset key
set xlabel "Number of electrons"
#set yrange [-0.495:-0.435]
#set ytics  -0.49,0.01,-0.44
set ylabel "t_CPU (s/block)"
plot f(x), "nelectron_tcpu.dat" using 1:2

pause -1

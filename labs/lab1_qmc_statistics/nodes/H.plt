#set terminal png size 400,300 
#set output 'nnode_vs_energy_H_STO-2G.png'

unset key
set xlabel "Number of nodes"
set logscale x
set yrange [-0.465:-0.445]
set ytics  -0.464,0.002,-0.446
set ylabel "E_total (Ha)"
plot[10:1000] "H.dat" using 1:5:7 with errorbars
pause -1

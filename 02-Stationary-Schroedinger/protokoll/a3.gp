set format x "%.0e"
set format y "%.0e"

set term png size 800,600

set xlabel "r [m]"
set ylabel "V(r) [eV]"

set output "img/wood-saxon.png"

plot "../results/pot_plot.txt" with lines title ""

pause -1

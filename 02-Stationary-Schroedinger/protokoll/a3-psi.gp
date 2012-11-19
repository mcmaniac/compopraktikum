
set ylabel "PSI(r) [bel. Einheiten]"

set xlabel "r [m]"

set term png size 800,600
set output "img/psi.png"

plot "../results/psi/0-1.txt"  using 1:2 with lines title "l=0, j=1", \
     "../results/psi/0-3.txt"  using 1:2 with lines title "l=0, j=3", \
     "../results/psi/1-2.txt"  using 1:2 with lines title "l=1, j=2", \
     "../results/psi/2-2.txt"  using 1:2 with lines title "l=2, j=2"

pause -1

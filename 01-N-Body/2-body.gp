set logscale y

set xlabel "t"
set ylabel "|E_0 - E(t)|"
plot "results/euler.txt" using 1:(abs($5)) title "Euler", \
     "results/euler_cromer.txt" using 1:(abs($5)) title "Euler-Cromer", \
     "results/leap_frog.txt" using 1:(abs($5)) title "Leap-Frog", \
     "results/verlet.txt" using 1:(abs($5)) title "Verlet", \
     "results/hermite.txt" using 1:(abs($5)) title "Hermite", \
     "results/hermite_iterated_2.txt" using 1:(abs($5)) title "2x Hermite", \
     "results/hermite_iterated_3.txt" using 1:(abs($5)) title "3x Hermite

pause -1

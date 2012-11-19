
set ylabel "PSI(r) [bel. Einheiten]"

set xlabel "r [m]"

set term png size 800,600

set title "PSI(r) für R = 5fm, l_max = 3"
set output "img/psi-var-5-bf-3.png"
plot "../results/psi-var-5-bf-3/0-0.txt"  using 1:2 with lines title "l=0, j=0", \
     "../results/psi-var-5-bf-3/0-1.txt"  using 1:2 with lines title "l=0, j=1", \
     "../results/psi-var-5-bf-3/1-0.txt"  using 1:2 with lines title "l=1, j=0", \
     "../results/psi-var-5-bf-3/2-0.txt"  using 1:2 with lines title "l=2, j=0"

pause -1

set title "PSI(r) für R = 15fm, l_max = 3"
set output "img/psi-var-15-bf-3.png"
plot "../results/psi-var-15-bf-3/0-2.txt"  using 1:2 with lines title "l=0, j=2", \
     "../results/psi-var-15-bf-3/0-6.txt"  using 1:2 with lines title "l=0, j=6", \
     "../results/psi-var-15-bf-3/1-6.txt"  using 1:2 with lines title "l=1, j=6", \
     "../results/psi-var-15-bf-3/2-5.txt"  using 1:2 with lines title "l=2, j=5"

pause -1

set title "PSI(r) für R = 5fm, l_max = 5"
set output "img/psi-var-5-bf-5.png"
plot "../results/psi-var-5-bf-5/0-0.txt"  using 1:2 with lines title "l=0, j=0", \
     "../results/psi-var-5-bf-5/0-1.txt"  using 1:2 with lines title "l=0, j=1", \
     "../results/psi-var-5-bf-5/1-0.txt"  using 1:2 with lines title "l=1, j=0", \
     "../results/psi-var-5-bf-5/2-0.txt"  using 1:2 with lines title "l=2, j=0"

pause -1

set title "PSI(r) für R = 10fm, l_max = 5"
set output "img/psi-var-10-bf-5.png"
plot "../results/psi-var-10-bf-5/0-1.txt"  using 1:2 with lines title "l=0, j=1", \
     "../results/psi-var-10-bf-5/0-3.txt"  using 1:2 with lines title "l=0, j=3", \
     "../results/psi-var-10-bf-5/1-2.txt"  using 1:2 with lines title "l=1, j=2", \
     "../results/psi-var-10-bf-5/2-2.txt"  using 1:2 with lines title "l=2, j=2"

pause -1

set title "PSI(r) für R = 15fm, l_max = 5"
set output "img/psi-var-15-bf-5.png"
plot "../results/psi-var-15-bf-5/0-2.txt"  using 1:2 with lines title "l=0, j=2", \
     "../results/psi-var-15-bf-5/0-6.txt"  using 1:2 with lines title "l=0, j=6", \
     "../results/psi-var-15-bf-5/1-6.txt"  using 1:2 with lines title "l=1, j=6", \
     "../results/psi-var-15-bf-5/2-5.txt"  using 1:2 with lines title "l=2, j=5"

pause -1

set title "PSI(r) für R = 5fm, l_max = 10"
set output "img/psi-var-5-bf-10.png"
plot "../results/psi-var-5-bf-10/0-0.txt"  using 1:2 with lines title "l=0, j=0", \
     "../results/psi-var-5-bf-10/0-1.txt"  using 1:2 with lines title "l=0, j=1", \
     "../results/psi-var-5-bf-10/1-0.txt"  using 1:2 with lines title "l=1, j=0", \
     "../results/psi-var-5-bf-10/2-0.txt"  using 1:2 with lines title "l=2, j=0"

pause -1

set title "PSI(r) für R = 10fm, l_max = 10"
set output "img/psi-var-10-bf-10.png"
plot "../results/psi-var-10-bf-10/0-1.txt"  using 1:2 with lines title "l=0, j=1", \
     "../results/psi-var-10-bf-10/0-3.txt"  using 1:2 with lines title "l=0, j=3", \
     "../results/psi-var-10-bf-10/1-2.txt"  using 1:2 with lines title "l=1, j=2", \
     "../results/psi-var-10-bf-10/2-2.txt"  using 1:2 with lines title "l=2, j=2"

pause -1

set title "PSI(r) für R = 15fm, l_max = 10"
set output "img/psi-var-15-bf-10.png"
plot "../results/psi-var-15-bf-10/0-2.txt"  using 1:2 with lines title "l=0, j=2", \
     "../results/psi-var-15-bf-10/0-6.txt"  using 1:2 with lines title "l=0, j=6", \
     "../results/psi-var-15-bf-10/1-6.txt"  using 1:2 with lines title "l=1, j=6", \
     "../results/psi-var-15-bf-10/2-5.txt"  using 1:2 with lines title "l=2, j=5"

pause -1


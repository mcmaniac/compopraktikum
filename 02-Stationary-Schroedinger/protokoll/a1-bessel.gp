set xrange [0:75]
set yrange [-0.25:1]
set key outside

set term png size 800,600

set output "img/bessel-01234.png"

plot "../results/bessel/0.txt"  with lines title "l = 0", \
     "../results/bessel/1.txt"  with lines title "l = 1", \
     "../results/bessel/2.txt"  with lines title "l = 2", \
     "../results/bessel/3.txt"  with lines title "l = 3", \
     "../results/bessel/4.txt"  with lines title "l = 4"

pause -1

set output "img/bessel-56789.png"

plot "../results/bessel/5.txt"  with lines title "l = 5", \
     "../results/bessel/6.txt"  with lines title "l = 6", \
     "../results/bessel/7.txt"  with lines title "l = 7", \
     "../results/bessel/8.txt"  with lines title "l = 8", \
     "../results/bessel/9.txt"  with lines title "l = 9"

pause -1

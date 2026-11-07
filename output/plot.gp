set terminal pngcairo size 900,600
set output 'plot.png'
set grid; set xlabel 'x'; set ylabel 'y'; set key left top
plot 'points.dat' with points pt 7 title 'y_i', 'curves.dat' using 1:2 with lines lw 2 title 'P1(x)', 'curves.dat' using 1:3 with lines lw 2 title 'P2(x)'

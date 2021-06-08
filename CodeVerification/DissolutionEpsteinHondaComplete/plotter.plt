#set terminal epslatex size 3.5,2.62 color
#set output 'nitrogen.tex'

set terminal png enhanced
set output 'nitrogen.png'


set key right bottom

set border linewidth 2
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 5
set tics scale 1.25

set ylabel 'Total dissolution time [sec]'
set xlabel 'Initial bubble radius [microns]'
set logscale xy
set grid
set mxtics
set mytics
set xrange [0.5e-6:5e-3]
set xtics ( "1" 0.000001, "100" 0.0001, "1000" 0.001)
set ytics (0.01, 100, 1000 )

plot 'Miyamoto_nitrogen.txt' u 1:2 w p pt 9 ps 1 title 'Exp: Miyamoto', 'Ljunggren_eriksson_nitrogen.txt' u 1:2 w p pt 7 ps 1 title 'Exp: Ljunggren et al.', 'Epstein_honda_nitrogen.txt' u 1:2 w p pt 5 ps 1 title 'Exp: Epstein-Honda', 'CodeResultsNitrogen.txt' u 1:2 w l title 'Predicted'

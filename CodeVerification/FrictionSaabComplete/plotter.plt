set terminal epslatex size 3.5,2.62 standalone color colortext
set output 'Saad.tex'


set key right bottom

set border linewidth 2
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 5
set tics scale 1.25

set xlabel 'Q (m^3/sec)'
set ylabel 'Friction Factor (\zeta)'



plot 'Saad_measurement.txt' u 1:($2/$4):($2/$4 * 0.05) w yerrorbars title 'Experiment', '' u 1:($3/$4):($3/$4 * 0.05) w yerrorbars title 'Predicted'

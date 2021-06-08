#set terminal epslatex size 3.5,2.62 standalone color colortext
#set output 'APdeJong.tex'

set terminal png enhanced
set output 'APdeJong.png'

set key right top

set border linewidth 2
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 5
set tics scale 1.25

set xlabel 'Q (m^3/sec)'
set ylabel 'Friction Factor'



plot 'DeJong_measure_compare.txt' u 1:2 w l title 'Experiment', '' u 1:3 w l title 'Predicted'

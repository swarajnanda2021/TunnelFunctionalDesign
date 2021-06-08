#set terminal epslatex size 3.5,2.62 standalone color colortext
#set output 'Melika.tex'

set terminal png enhanced
set output 'Melika.png'

set key left top

set border linewidth 2
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 5
set tics scale 1.25

set xlabel 'Q (m^3/sec)'
set ylabel 'Pressure drop (Pa)'



plot 'pipe_data.txt' u ($4*0.001/3600):2 w l title 'Experiment', '' u ($4*0.001/3600):5 w l title 'Predicted'

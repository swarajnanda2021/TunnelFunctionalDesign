set terminal png enhanced
set output 'Wetzel_turbmgmt_validation.png'


set key right top

set border linewidth 2
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 5
set tics scale 1.25

set xlabel 'Distance from contraction (x/L)'
set ylabel 'T.I. (%)'

set grid 

plot 'wetzel2p54.txt' u 1:2 w p pt 7 ps 2 title 'Wetzel et al. (M=2.54 cm)', '' u 1:3 w p pt 1 ps 4 title 'Predicted (M=2.54 cm)',\
'wetzel5p04.txt' u 1:2 w p pt 7 ps 2 title 'Wetzel et al. (M=5.04 cm)', '' u 1:3 w p pt 1 ps 4 title 'Predicted (M=5.04 cm)'
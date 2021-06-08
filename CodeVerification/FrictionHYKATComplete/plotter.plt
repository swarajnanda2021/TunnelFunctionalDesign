#set terminal epslatex size 3.5,2.62 standalone color colortext
set term png enhanced
set output 'hykat.png'

set style data histogram
set style fill solid border


set key left top

set border linewidth 2
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 5
set tics scale 1.25

set ylabel 'Pressure drop (m)'



plot 'HYKAT_PressureDrop_Contribution.txt' u 2:xticlabels(1) title 'Wetzel et al.', '' u 3:xticlabels(1) title 'Predicted'

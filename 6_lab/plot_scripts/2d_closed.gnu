set terminal png size 700,700
set size square

set output "plots/".outputName.".png"

set xlabel "x"
set ylabel "y"
set title "Rozkład cząstek"
set xrange [-6:6]
set yrange [-6:6]

set object circle at 3,0 size Ra fillstyle solid fillcolor rgb "red"
set object circle at -4.5,0 size 0.17 fillstyle solid fillcolor rgb "blue"
set object circle at 0,0 size 5 fc rgb "navy"

plot "tmp_filename" u 1:2 lw 1 pt 7 lt rgb "black" notitle


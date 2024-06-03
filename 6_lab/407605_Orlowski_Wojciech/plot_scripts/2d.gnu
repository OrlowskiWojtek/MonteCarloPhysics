set terminal png size 700,700
set size square

set output "plots/".outputName.".png"

set xlabel "x"
set ylabel "y"
set title "Rozkład cząstek"

plot "tmp_filename" u 1:2 lw 1 pt 7 lt rgb "black" notitle


set terminal png size 700,700
set size square

set output "plots/".outputName.".png"

set xlabel "x"
set ylabel "y"
set title "Rozkład cząstek"
set xrange [-30:30]
set yrange [-30:30]

plot "tmp_filename" u 5:6 lw 1 pt 7 title "t = 20.0", "" u 3:4 lw 1 pt 7 title "t = 5.0" ,\
    "" u 1:2 lw 1 pt 7 title "t = 1.0";


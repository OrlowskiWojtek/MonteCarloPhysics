set terminal png size 1024,768

set output "plots/".outputName.".png"

set xlabel "t"
set ylabel "d(t)"
set title "Współczynnik dyfuzji"

plot "tmp_filename" u 1:2 lw 2 w l title "D_{xx}", "tmp_filename" u 1:3 lw 2 w l title "D_{yy}", "tmp_filename" u 1:4 lw 2 w l title "D_{xy}", \
    "" u 1:(0) w l notitle lt rgb "black", "" u 1:(1) w l notitle lt rgb "black" 


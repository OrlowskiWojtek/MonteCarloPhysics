set terminal png size 1024,768

set output "plots/".outputName.".png"

set xlabel "Iteracja"
set ylabel "Ilość cząstek"
set title "Ilość cząstek w symulacji w zależności od iteracji"

plot "tmp_filename" u 1 w l lw 2 lt rgb "black" notitle 

set terminal dumb size 60,20 ansi

set colorsequence classic


plot "data/hist.dat" using 1:2 with lines title "mc", "" using 1:3 with lines title "dokladny"
pause 3
reread

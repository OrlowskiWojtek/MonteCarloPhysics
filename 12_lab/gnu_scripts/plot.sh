#!/bin/bash

# Ścieżka do pliku z danymi
DATA_FILE="data/hist.dat"

# Pętla nieskończona do ciągłego rysowania wykresów
while true; do
  # Wyczyść ekran terminala
  clear

  # Użyj gnuplota do rysowania wykresu w terminalu dumb
  gnuplot -e "set terminal dumb ansi; set colorsequence classic ;plot '$DATA_FILE' using 1:2 with lines title 'mc', '' u 1:3 with lines title 'dokładny' "

  # Przerwa 1 sekunda
  sleep 1
done

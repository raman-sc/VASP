#!/usr/bin/bash

gnuplot << EOF
reset
set terminal postscript enhanced "Helvetica" 20
set output "Raman.ps"
set xrange[0:1000]
set style line 1 lt 1 lc -1 lw 4 # black-solid
set style line 2 lt 2 lc 1 lw 4 # red-dashed
set title "Off-resonance Raman: PAW @ PBE level\nSolid - Li1Si64-Td; Dashed - Si64-diamond bulk"
#set nokey
set xlabel "Energy, cm^{-1}"
#set ylabel "Intensity, a.u."
plot "vasp_raman.dat-a-pw91.out" u 1:2 w l ls 1,  "Raman-data-a-pw91.out" u 1:2 w l ls 2 
EOF

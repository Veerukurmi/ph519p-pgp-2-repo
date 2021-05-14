#1 au = a_0 = 0.529177 angstorm        
#1 a_0^2 = 0.280028297 A^2, a_0 = bohr radius  

reset
set terminal postscript eps enhanced "Helvetica" 18 color
set output "l_6.eps"
set size 1,1





set encoding iso_8859_1
#set multiplot layout 2,1
set auto fix
show label

set tics scale 2

#First main plot




set xlabel "Electron Energy (a.u.)" font "Times-Roman,24"
#set xrange[15:105]
#set yrange [-2.0:7.5]

set y2tics
#set logscale y2
set mxtics
set ytics nomirror

set ylabel "Time delay (a.u.)" font "Times-Roman,18"
set y2label "Phase-Shift in radians" font "Times-Roman,18"

set xrange[-0.01:0.55]
set yrange[-12:8]
set y2range[0:1.2]

#set label 1 "Phase-Shift (radians)" at screen 0.05,1 center front rotate font "Times-Roman,24"
#set label 2 "Na_{20}" font "Times-Roman,26" at 7,100 center
#set label 3 "(b)" font "Times-Roman,26" at 7,1e-5 center

#set format y2 "10^{%L}"



set key at graph .75, .25 spacing 1.5 font "Times-Roman, 24"
set key box width -0.5
set key box lt -1 lw 2

plot "l_6_matlab_phase_shift.txt" using 1:2 w lp lw 4.0 lc 'red' title "l=6" axis x1y1,"phase_l_6.txt" using 1:2 w lp lw 4.0 dt 3.0 lc 'black' title "l=6" axis x1y2










#set xrange [0:10]
set term x11 1

p "./dmft.005.dat" u 1:2 w l title "G", "" u 1:4 w l title "G0", "" u 1:6 w l title "Sigma"
set term x11 2
p "./dmft.005.dat" u 1:3 w l title "G", "" u 1:5 w l title "G0", "" u 1:7 w l title "Sigma"

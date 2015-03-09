 set term x11 enhanced font "arial,15"
 set output "plot1.png"        
 set border linewidth 0
 set lmargin screen 0.1
 set rmargin screen 0.9
 set tmargin screen 0.9
 set bmargin screen 0.1
 set palette maxcolors 2
 set palette defined ( -1 "#0066ff", 1 "#ff3300")
 set cbrange [-1:1]
 set cbtics ("+" 1, "-" -1)
 set title "initial state"
 set pm3d map
 count = 0
 load "loop.plt"

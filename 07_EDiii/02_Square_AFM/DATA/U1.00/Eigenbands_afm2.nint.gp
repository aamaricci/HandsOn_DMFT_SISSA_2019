 #set term wxt
 #set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
 #set out 'Eigenbands_afm2.nint.png'
 
 #set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'
 #set out 'Eigenbands_afm2.nint.svg'
 
 #set term postscript eps enhanced color 'Times'
 #set output '|ps2pdf -dEPSCrop - Eigenbands_afm2.nint.pdf'
 unset key
 set xtics ('G'0.000000,'X'3.141593,'M'6.283185,'G'10.726068)
 set grid noytics xtics
set label 'Orb 1' tc rgb 16711680 at graph 0.9,0.950000 font 'Times-Italic,11'
set label 'Orb 2' tc rgb 255 at graph 0.9,0.900000 font 'Times-Italic,11'
set label 'Orb 3' tc rgb 16711680 at graph 0.9,0.850000 font 'Times-Italic,11'
set label 'Orb 4' tc rgb 255 at graph 0.9,0.800000 font 'Times-Italic,11'
 plot 'Eigenbands_afm2.nint' u 1:2:3 w l lw 3 lc rgb variable,\
 'Eigenbands_afm2.nint' u 1:4:5 w l lw 3 lc rgb variable,\
 'Eigenbands_afm2.nint' u 1:6:7 w l lw 3 lc rgb variable,\
 'Eigenbands_afm2.nint' u 1:8:9 w l lw 3 lc rgb variable

#unset key
set terminal postscript eps enhanced color font "Times-Roman,16"
set output "|ps2pdf -sEPSCrop - spH0dw.pdf"
set size ratio -1
set xlabel "<--- J --->"
set ylabel "<--- I --->"
set title "    60 nonzeros for spH0dw"
set timestamp
plot [x=1:20] [y=20:1] "spH0dw_data.dat" w p pt 5 ps 0.4 lc rgb "red"

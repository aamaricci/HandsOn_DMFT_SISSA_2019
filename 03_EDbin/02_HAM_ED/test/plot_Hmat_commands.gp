#unset key
set terminal postscript eps enhanced color font "Times-Roman,16"
set output "|ps2pdf -sEPSCrop - Hmat.pdf"
set size ratio -1
set xlabel "<--- J --->"
set ylabel "<--- I --->"
set title "  2800 nonzeros for Hmat"
set timestamp
plot [x=1:400] [y=400:1] "Hmat_data.dat" w p pt 5 ps 0.4 lc rgb "red"

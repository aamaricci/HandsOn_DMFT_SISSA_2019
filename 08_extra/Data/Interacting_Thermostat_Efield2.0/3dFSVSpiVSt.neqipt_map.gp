 reset
 #set term gif animate
 #set output '3dFSVSpiVSt.neqipt.gif'
 set term wxt
 set pm3d map
 set size square
 set xrange [0.000000:5.983986]
 set yrange [0.000000:5.983986]
 set cbrange [-.018174:0.999836]
 n=100
 do for [i=0:n-1]{
 set title sprintf('3dFSVSpiVSt.neqipt; i=%i',i+1)
 splot '3dFSVSpiVSt.neqipt' every :::(21*i)::(21*i+20) title '' 
 }
 #set output

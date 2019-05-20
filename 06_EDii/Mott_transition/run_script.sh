HERE=$(pwd)
while read u;
do 
    echo $U;
    DIR=U$u;
    mkdir $DIR;
    cp inputED.in *.restart $DIR/;
    cd $DIR;
    ed_hm_bethe uloc=$u > OUT 2>&1;
    cp *.restart $HERE/;
    cd $HERE;
done <list_u

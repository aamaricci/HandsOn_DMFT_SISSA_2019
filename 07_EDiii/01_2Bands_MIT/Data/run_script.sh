HERE=$(pwd)
while read u;
do
    ULOC=$u
    UST=$(echo $u/2|bc -l)
    JH=$(echo $u/4|bc -l)
    echo "Uloc=$ULOC"
    echo "Ust =$UST"
    echo "Jh  =$JH"
    DIR=U$u;
    if [ ! -d $DIR ];then
	mkdir -pv $DIR;
	cp inputED.conf *.restart $DIR/;
	cd $DIR;
	edn_hm_2b_square uloc=$ULOC,$ULOC ust=$UST jh=$JH |tee LOG.out #> OUT 2>&1;
	cp *.restart $HERE/;
	cd $HERE;
    else
	echo "$DIR exists: skip"
	cp $DIR/*.restart $HERE/
    fi
done <list_u

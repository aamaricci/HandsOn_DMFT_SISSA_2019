HERE=$(pwd)
while read u;
do
    ULOC=$u
    echo "Uloc=$ULOC"
    DIR=U$u;
    if [ ! -d $DIR ];then
	mkdir -pv $DIR;
	cp inputED.conf *.restart $DIR/;
	cd $DIR;
	edn_hm_square_afm2 uloc=$ULOC |tee LOG.out
	cp *.restart $HERE/
	cd $HERE;
    else
	echo "$DIR exists: skip"
	cp $DIR/*.restart $HERE/
    fi
done <list_u

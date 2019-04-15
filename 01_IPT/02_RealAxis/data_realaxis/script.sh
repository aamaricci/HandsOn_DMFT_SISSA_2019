while read U;do
    DIR=U$U;
    mkdir -p $DIR;
    cp -v *restart inputIPT.conf $DIR/;
    cd $DIR/
    ipt_hm_real uloc=$U |tee LOG.out
    cp *.restart ../
    cd ../
done<list_u

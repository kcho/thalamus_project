for side in lh rh
do
    for i in ../[CN]*
    do
        echo bash 9_thalamus_to_whole_brain_tractography.sh ${i} ${side}
    done
done|parallel -j 20

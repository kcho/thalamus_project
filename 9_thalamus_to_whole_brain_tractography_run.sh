for i in ../[CN]*
do
    #echo bash 2_segmentation.sh ${i} lh

    echo bash 9_thalamus_to_whole_brain_tractography.sh ${i} lh
    echo bash 9_thalamus_to_whole_brain_tractography.sh ${i} rh

done|parallel -j 20

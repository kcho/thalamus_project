for i in ../[CN]*
do
    #echo bash 2_segmentation.sh ${i} lh
    echo bash 2_segmentation.sh ${i} rh
done|parallel

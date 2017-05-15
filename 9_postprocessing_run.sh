for i in ../[CNF]*
do
    echo bash 9_postprocessing.sh ${i}
done|parallel -j 12


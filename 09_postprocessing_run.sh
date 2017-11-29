for i in ../[CNF]*
do
    echo ${i}
    #echo bash 9_postprocessing.sh ${i}/thalamus_tractography_MNI
    bash 9_postprocessing.sh ${i}
done


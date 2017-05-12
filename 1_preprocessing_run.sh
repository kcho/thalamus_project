for i in ../[CN]*
do 
    echo bash 1_preprocessing.sh ${i}
done|parallel

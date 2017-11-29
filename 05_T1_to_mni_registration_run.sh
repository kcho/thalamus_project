for i in ../[CN]*
do
    echo bash 5_T1_to_mni_registration.sh ${i}
done|parallel

for i in ../[CN]*
do
    echo bash 7_DKI_to_mni_registration.sh ${i}
done|parallel

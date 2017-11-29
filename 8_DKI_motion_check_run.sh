for i in $@
do
    echo bash diff_motion_check.sh `ls ${i}/DKI/*.ecclog`
done|parallel

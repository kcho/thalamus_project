for i in $@
do
    echo diff_motion_check.sh `ls ${i}/DKI/*.ecclog`
done|parallel

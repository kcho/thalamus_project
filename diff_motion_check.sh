#!/bin/sh

# https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;e2a42b99.1101
# by Mark Jenkinson
# edited by kcho - output location into the eddy ecclog location

if [ $# -lt 1 ] ; then 
    echo "Usage: `basename $0` <eddy current ecclog file>"
    exit 1;
fi

logfile=$1;
dir=`dirname $logfile`;
basenm=`basename $logfile .ecclog`;

nums=`grep -n 'Final' $logfile | sed 's/:.*//'`; 

#touch ${dir}/grot_ts.txt
touch ${dir}/grot.mat

firsttime=yes;
m=1;

for tmpFile in disp rot trans
do
    tmpFileLoc=${dir}/ec_${tmpFile}.txt 
    if [ -e ${tmpFileLoc} ]
    then
        rm ${tmpFileLoc}
    fi
done

for n in $nums ; do 
    echo "Timepoint $m"
    n1=`echo $n + 1 | bc` ; 
    n2=`echo $n + 5 | bc` ;
    echo $n1 $n2
    sed -n  "${n1},${n2}p" $logfile > ${dir}/grot.mat ; 

    if [ $firsttime = yes ] 
    then 
        firsttime=no
        cp ${dir}/grot.mat ${dir}/grot.refmat
        cp ${dir}/grot.mat ${dir}/grot.oldmat
    fi
    #absval=`$FSLDIR/bin/rmsdiff grot.mat grot.refmat $basenm`;
    #relval=`$FSLDIR/bin/rmsdiff grot.mat grot.oldmat $basenm`;
    absval=`$FSLDIR/bin/rmsdiff ${dir}/grot.mat ${dir}/grot.refmat ${dir}/${basenm}`;
    relval=`$FSLDIR/bin/rmsdiff ${dir}/grot.mat ${dir}/grot.oldmat ${dir}/${basenm}`;
    cp ${dir}/grot.mat ${dir}/grot.oldmat
    echo $absval $relval >> ${dir}/ec_disp.txt ;

    $FSLDIR/bin/avscale --allparams ${dir}/grot.mat ${dir}/${basenm} | grep 'Rotation Angles' | sed 's/.* = //' >> ${dir}/ec_rot.txt ;

    $FSLDIR/bin/avscale --allparams ${dir}/grot.mat ${dir}/${basenm} | grep 'Translations' | sed 's/.* = //' >> ${dir}/ec_trans.txt ;

    m=`echo $m + 1 | bc`;
done

echo "absolute" > ${dir}/grot_labels.txt
echo "relative" >> ${dir}/grot_labels.txt

$FSLDIR/bin/fsl_tsplot -i ${dir}/ec_disp.txt -t 'Eddy Current estimated mean displacement (mm)' -l ${dir}/grot_labels.txt -o ${dir}/ec_disp.png

echo "x" > ${dir}/grot_labels.txt
echo "y" >> ${dir}/grot_labels.txt
echo "z" >> ${dir}/grot_labels.txt

$FSLDIR/bin/fsl_tsplot -i ${dir}/ec_rot.txt -t 'Eddy Current estimated rotations (radians)' -l ${dir}/grot_labels.txt -o ${dir}/ec_rot.png
$FSLDIR/bin/fsl_tsplot -i ${dir}/ec_trans.txt -t 'Eddy Current estimated translations (mm)' -l ${dir}/grot_labels.txt -o ${dir}/ec_trans.png

# clean up temp files
/bin/rm ${dir}/grot_labels.txt ${dir}/grot.oldmat ${dir}/grot.refmat ${dir}/grot.mat 

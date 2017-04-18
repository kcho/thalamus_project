#!/bin/sh
Usage() {
    echo ""
    echo "Usage: 2_segmentation <subjDir> <lh>"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

# edit here for different folder structure
subj=${1}
fsDir=${subj}/FREESURFER
regDir=${subj}/registration
dtiDir=${subj}/DTI
dkiDir=${subj}/DKI
roiDir=${subj}/ROI
tractDir=${subj}/thalamus_tractography
segDir=${subj}/segmentation
bedpostDir=${subj}/DTI.bedpostX
side_s=$2

if [ $2 == 'lh' ]
then
    side=left
elif [ $2 == 'rh' ]
then
    side=right
else
    Usage
fi

#Segmentation
if [ ! -d ${tractDir} ]
then
    mkdir ${tractDir}
else
    if [ ! -d ${tractDir}/${side} ]
    then
        mkdir ${tractDir}/${side}
    fi
fi

# Freesurfer ROI registration 
fsThalROI=${regDir}/${side_s}_thalamus.nii.gz 
dtiThalROI=${regDir}/${side_s}_thalamus_DTI.nii.gz 

if [ ! -e ${dtiThalROI} ]
then 
    flirt \
        -in ${fsThalROI} \
        -ref ${bedpostDir}/nodif_brain_mask.nii.gz \
        -applyxfm -init ${regDir}/FREESURFERT1toNodif.mat \
        -interp nearestneighbour \
        -out ${dtiThalROI}
fi

if [ ! -e ${tractDir}/${side}/fdt_paths.nii.gz ]
then
    rm -rf ${tractDir}/${side} 
    mkdir -p ${tractDir}/${side}
    probtrackx2 \
        -x ${dtiThalROI} \
        -l \
        --onewaycondition \
        --omatrix2 \
        --target2=${bedpostDir}/nodif_brain_mask.nii.gz \
        -c 0.2 \
        -S 2000 \
        --steplength=0.5 \
        -P 5000 \
        --fibthresh=0.01 \
        --distthresh=0.0 \
        --sampvox=0.0 \
        --forcedir \
        --opd \
        -s ${bedpostDir}/merged \
        -m ${bedpostDir}/nodif_brain_mask \
        --dir=${tractDir}/${side}
    echo ${subj} thalamo-whole brain tractography on the ${side} done
else
    echo ${subj} thalamo-whole brain tractography on the ${side} done
fi

#!/bin/sh

Usage() {
    echo ""
    echo "Usage: 9_postprocessing.sh <subject_space_tractography_dir>"
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

for side in left right
do
    reconImg=${tractDir}/${side}/fdt_matrix2_recontructed.nii.gz
    reconImgMNI=${tractDir}/${side}/fdt_matrix2_recontructed_MNI.nii.gz
    reconImgMNI_ds=${tractDir}/${side}/fdt_matrix2_recontructed_MNI_ds.nii.gz
    #Segmentation
    if [ ! -e ${reconImg} ]
    then 
        python tracktography/postprocessing/probtrackx_postprocessing.py ${tractDir}/${side}
    fi

    if [ ! -e ${reconImgMNI} ]
    then
        flirt \
            -in ${reconImg} \
            -ref ${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz \
            -applyxfm -init ${regDir}/nodifToMNI.mat \
            -out ${reconImgMNI}
    fi

    if [ ! -e ${reconImgMNI_ds} ]
    then
        flirt \
            -in ${reconImgMNI} \
            -ref ${reconImgMNI}
            -applyisoxfm 3 \
            -out ${reconImgMNI_ds}
    fi
done
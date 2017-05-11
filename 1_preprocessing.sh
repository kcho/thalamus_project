#!/bin/bash
for subj in $@
do
    # edit here for different folder structure
    fsDir=${subj}/FREESURFER
    regDir=${subj}/registration
    dtiDir=${subj}/DTI
    dkiDir=${subj}/DKI
    roiDir=${subj}/ROI
    segDir=${subj}/segmentation
    bedpostDir=${subj}/DTI.bedpostX

    if [ ! -d ${bedpostDir} ]
    then
        bedpostx ${dtiDir}
    fi

    # convert the Freesurfer output : brain.mgz --> brain.nii.gz 
    if [ ! -e ${fsDir}/mri/brain.nii.gz ]
    then
        mri_convert --out_orientation RAS \
            ${fsDir}/mri/brain.mgz \
            ${fsDir}/mri/brain.nii.gz
    fi

    # Registration
    if [ ! -d ${regDir} ]
    then
        mkdir ${regDir}
    fi

    # FREESURFER output T1 --> diffusion tensor space
    if [ ! -e ${regDir}/FREESURFERT1toNodif.mat ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref ${dtiDir}/nodif_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toNodif \
            -omat ${regDir}/FREESURFERT1toNodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # FREESURFER output T1 --> diffusion kurtosis space
    if [ ! -e ${regDir}/FREESURFERT1toDKINodif.mat ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref ${dkiDir}/nodif_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toDKINodif \
            -omat ${regDir}/FREESURFERT1toDKINodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # FREESURFER output T1 --> MNI 1mm
    if [ ! -e ${regDir}/FREESURFERT1toMNI.mat ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref /usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toMNI \
            -omat ${regDir}/FREESURFERT1toMNI.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # DTI --> FREESURFER output T1 --> MNI 1mm

    if [ ! -e ${regDir}/nodifToFreesurfer.mat ]
    then
        convert_xfm \
            -omat ${regDir}/nodifToFreesurfer.mat \
            -inverse ${regDir}/FREESURFERT1toNodif.mat
    fi

    if [ ! -e ${regDir}/nodifToMNI.mat ]
    then
        convert_xfm -omat ${regDir}/nodifToMNI.mat \
            -concat ${regDir}/FREESURFERT1toMNI.mat ${regDir}/nodifToFreesurfer.mat
    fi

    # ROI extraction
    if [ ! -d ${roiDir} ] 
    then
        mkdir ${roiDir}
    fi

    if [ ! -e ${roiDir}/rh_OCC.nii.gz ]
    then
        python roiExtraction.py \
            -s ${subj} \
            -r `basename ${roiDir}` \
            -f `basename ${fsDir}`
    else
        echo ${subj} 'ROI extraction done'
    fi
done

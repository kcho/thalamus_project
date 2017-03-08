#!/bin/bash
for subj in $@
do
    # edit here for different folder structure
    fsDir=${i}/FREESURFER
    regDir=${i}/registration
    dtiDir=${i}/DTI
    dkiDir=${i}/DKI
    roiDir=${i}/ROI
    segDir=${i}/segmentation
    bedpostDir=${i}/DTI.bedpostX

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
            -ref ${dtiDir}/b0_brain.nii.gz \
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
            -ref ${dkiDir}/b0_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toDKINodif \
            -omat ${regDir}/FREESURFERT1toDKINodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi


    # ROI extraction
    if [ ! -d ${roiDir} ] 
    then
        mkdir ${roiDir}
    fi

    if [ ! -e ${roiDir}/rh_OCC.nii.gz ]
    then
        python roiExtraction.py -S ${subj}
    else
        echo ${subj} 'ROI extraction done'
    fi

done

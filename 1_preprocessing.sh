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

    # MNI 2mm --> FREESURFER output T1 
    if [ ! -e ${regDir}/MNItoFREESURFER.mat ]
    then
        flirt \
            -in /usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz \
            -ref ${fsDir}/mri/brain.nii.gz \
            -out ${regDir}/MNItoFREESURFER \
            -omat ${regDir}/MNItoFREESURFER.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # Estimation of transformation matrix
    # - MNI to DTI
    # - This is more accurately done by going through DTI --> FREESURFER output T1 --> MNI 1mm
    # Usage: convert_xfm [options] <input-matrix-filename>
    # e.g. convert_xfm -omat <outmat> -inverse <inmat>
    #        convert_xfm -omat <outmat_AtoC> -concat <mat_BtoC> <mat_AtoB>

    # MNI --> FREESURFER --> DTI
    if [ ! -e ${regDir}/MNItoNodif.mat ]
    then
        convert_xfm -omat ${regDir}/MNItoNodif.mat \
            -concat ${regDir}/FREESURFERT1toNodif.mat ${regDir}/MNItoFREESURFER.mat
    fi

    # MNI --> FREESURFER --> DTI
    if [ ! -e ${regDir}/nodifToMNI.mat ]
    then
        convert_xfm -omat ${regDir}/nodifToMNI.mat \
            -inverse ${regDir}/MNItoNodif.mat
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

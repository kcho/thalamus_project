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

    mni2mm=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
    mni2mmMask=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz


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
    fs2nodifMat=${regDir}/FREESURFERT1toNodif.mat
    if [ ! -e ${fs2nodifMat} ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref ${dtiDir}/nodif_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toNodif \
            -omat ${fs2nodifMat} \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # FREESURFER output T1 --> diffusion kurtosis space
    fs2dkinodifMat=${regDir}/FREESURFERT1toDKINodif.mat
    if [ ! -e ${fs2dkinodifMat} ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref ${dkiDir}/nodif_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toDKINodif \
            -omat ${fs2dkinodifMat} \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # MNI 2mm --> FREESURFER output T1 
    mni2fs=${regDir}/MNItoFREESURFER.mat
    if [ ! -e ${mni2fs} ]
    then
        flirt \
            -in ${mni2mm} \
            -ref ${fsDir}/mri/brain.nii.gz \
            -out ${regDir}/MNItoFREESURFER \
            -omat ${mni2fs} \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    #============================================================
    # FNIRT MNI --> Freesurfer
    #============================================================
    mni2fs_fnirt=${regDir}/MNI_to_FREESURFER_fnirt_coeff.nii.gz
    if [ ! -e ${mni2fs_fnirt} ]
    then
        fnirt \
            --ref=${fsDir}/mri/brain.nii.gz \
            --in=${mni2mm} \
            --aff=${mni2fs} \
            --inmask=${mni2mmMask} \
            --refmask=${fsDir}/mri/brainmask.nii.gz \
            --iout=${regDir}/MNI_to_FREESURFER_fnirt \
            --cout=${mni2fs_fnirt}
    fi
    
    # Freesurfer --> MNI
    mni2fs2nodif=${regDir}/MNI_to_FREESURFER_fnirt_to_Nodif_flirt.nii.gz
    if [ ! -e ${mni2fs2nodif} ]
    then
        convertwarp \
            --ref=${dtiDir}/nodif_brain_mask.nii.gz \
            --warp1=${mni2fs_fnirt} \
            --postmat=${fs2nodifMat} \
            --out=${mni2fs2nodif}
    fi


    nodif2fs2mni=${regDir}/Nodif_to_FREESURFER_flirt_to_MNI_fnirt.nii.gz
    if [ ! -e ${nodif2fs2mni} ]
    then
        invwarp \
            --ref=${mni2mm} \
            --warp=${mni2fs2nodif} \
            --out=${nodif2fs2mni}
    fi

    mni2nodif_fnfl_img=${regDir}/MNI_in_Nodif_FNIRT_FLIRT.nii.gz
    if [ ! -e ${mni2nodif_fnfl_img} ]
    then
        applywarp \
            --ref=${dtiDir}/nodif_brain.nii.gz \
            --in=${mni2mm} \
            --warp=${nodif2fs2mni} \
            --out=${mni2nodif_fnfl_img}
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
            -concat ${fs2nodifMat} ${mni2fs}
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

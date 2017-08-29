#!/bin/bash
for subj in $@
do
    # Path variables
    t1Dir=${subj}/T1

    # Freesurfer
    fsDir=${subj}/FREESURFER
    fsDirName=${fsDir##*/}
    fsMriDir=${fsDir}/mri
    fsT1Img=${fsMriDir}/T1.nii.gz
    fsT1BetImg=${fsMriDir}/brain.nii.gz
    fsT1MaskImg=${fsMriDir}/brainmask.nii.gz
    # MNI
    mni2mmImg=${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz
    mni2mmBetImg=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
    mni2mmMaskImg=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz 
    # DTI
    dtiDir=${subj}/DTI
    dtiNodifBetImg=${dtiDir}/nodif_brain.nii.gz
    dtiNodifImg=${dtiDir}/nodif.nii.gz
    dtiNodifMaskImg=${dtiDir}/nodif_brain_mask.nii.gz 
    # DKI
    dkiDir=${subj}/DKI
    dkiNodifBetImg=${dkiDir}/nodif_brain.nii.gz
    dkiNodifImg=${dkiDir}/nodif.nii.gz
    dkiNodifMaskImg=${dkiDir}/nodif_brain_mask.nii.gz
    # other dirs
    regDir=${subj}/registration
    roiDir=${subj}/ROI
    bedpostDir=${subj}/DTI.bedpostX

    # Recon-all
    if [ ! -d ${fsDir} ]
    then
        export SUBJECTS_DIR=${subj}
        recon-all \
            -all \
            -subjid ${fsDirName} \
            -i ${t1Dir}/20*.nii.gz
    fi

    # freesurfer mgz --> nii.gz
    if [ ! -e ${fsT1BetImg} ]
    then
        mri_convert \
            --out_orientation RAS \
            ${fsMriDir}/brain.mgz \
            ${fsT1BetImg}
    fi

    if [ ! -e ${fsT1Img} ]
    then
        mri_convert \
            --out_orientation RAS \
            ${fsMriDir}/T1.mgz \
            ${fsT1Img}
    fi

    if [ ! -e ${fsT1MaskImg} ]
    then
        mri_convert \
            --out_orientation RAS \
            ${fsMriDir}/brainmask.mgz \
            ${fsT1MaskImg}
        fslmaths \
            ${fsT1MaskImg} \
            -bin ${fsT1MaskImg}
    fi

    # Registration
    if [ ! -d ${regDir} ]
    then
        mkdir ${regDir}
    fi

    # flirt : freesurfer --> DTI/nodif_brain.nii.gz
    fs2dti=${regDir}/fs2dti.mat
    if [ ! -e ${fs2dti} ]
    then
        flirt \
            -in ${fsT1BetImg} \
            -ref ${dtiNodifBetImg} \
            -omat ${fs2dti} \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # fnirt : freesurfer nobet T1 --> DTI/nodif.nii.gz
    fs2dti_fnirt=${regDir}/fs2dti_fnirt_coeff.nii.gz
    fs2dti_fnirt_inv=${regDir}/fs2dti_fnirt_coeff_inv.nii.gz
    if [ ! -e ${fs2dti_fnirt} ]
    then
        fnirt \
            --ref=${dtiNodifImg} \
            --in=${fsT1Img} \
            --aff=${fs2dti} \
            --inmask=${fsT1MaskImg} \
            --refmask=${dtiNodifMaskImg} \
            --cout=${fs2dti_fnirt}
    fi

    if [ ! -e ${fs2dti_fnirt_inv} ]
    then
        invwarp \
            -w ${fs2dti_fnirt} \
            -o ${fs2dti_fnirt_inv} \
            -r ${fsT1Img}
    fi

    # flirt : freesurfer --> DKI/nodif_brain.nii.gz
    fs2dki=${regDir}/fs2dki.mat
    if [ ! -e ${fs2dki} ]
    then
        flirt \
            -in ${fsT1BetImg} \
            -ref ${dkiNodifBetImg} \
            -omat ${fs2dki} \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # FREESURFER output t1 --> diffusion kurtosis space fnirt
    fs2dki_fnirt=${regDir}/fs2dki_fnirt_coeff.nii.gz
    fs2dki_fnirt_inv=${regDir}/fs2dki_fnirt_coeff_inv.nii.gz
    if [ ! -e ${fs2dki_fnirt} ]
    then
        fnirt \
            --ref=${dkiNodifImg} \
            --in=${fsT1Img} \
            --aff=${fs2dki} \
            --inmask=${fsT1MaskImg} \
            --refmask=${dkiNodifMaskImg} \
            --cout=${fs2dki_fnirt}
    fi

    if [ ! -e ${fs2dki_fnirt_inv} ]
    then
        invwarp \
            -w ${fs2dki_fnirt} \
            -o ${fs2dki_fnirt_inv} \
            -r ${fsT1Img}
    fi

    # MNI 2mm --> FREESURFER output T1 
    mni2fs=${regDir}/mni2fs.mat
    if [ ! -e ${mni2fs} ]
    then
        flirt \
            -in ${mni2mmImg} \
            -ref ${fsT1BetImg} \
            -omat ${mni2fs} \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 12 -interp trilinear
    fi
    mni2fs_fnirt=${regDir}/mni2fs_fnirt_coeff.nii.gz
    if [ ! -e ${mni2fs_fnirt} ]
    then
        fnirt \
            --ref=${fsT1BetImg} \
            --in=${mni2mmImg} \
            --config=T1_2_MNI152_2mm \
            --aff=${mni2fs} \
            --inmask=${mni2mmMaskImg} \
            --refmask=${fsT1MaskImg} \
            --cout=${mni2fs_fnirt}
    fi
    # MNI --fnirt--> Freesurfer --flirt--> nodif
    mni2fs2nodif=${regDir}/mni2fs2nodif_coeff.nii.gz
    if [ ! -e ${mni2fs2nodif} ]
    then
        convertwarp \
            --ref=${dtiNodifMaskImg} \
            --warp1=${mni2fs_fnirt} \
            --postmat=${fs2nodifMat} \
            --out=${mni2fs2nodif}
    fi

    # nodif --flirt--> Freesurfer --fnirt--> MNI
    nodif2fs2mni=${regDir}/nodif2fs2mni_coeff.nii.gz
    if [ ! -e ${nodif2fs2mni} ]
    then
        invwarp \
            --ref=${mni2mmBetImg} \
            --warp=${mni2fs2nodif} \
            --out=${nodif2fs2mni}
    fi

    # Check ouput image
    mni2nodif_img=${regDir}/mni2nodif_img.nii.gz
    if [ ! -e ${mni2nodif_img} ]
    then
        applywarp \
            --ref=${dtiNodifMaskImg} \
            --in=${mni2mmBetImg} \
            --warp=${mni2fs2nodif} \
            --out=${mni2nodif_img}
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

    if [ ! -d ${bedpostDir} ]
    then
        bedpostx ${dtiDir}
    fi
done

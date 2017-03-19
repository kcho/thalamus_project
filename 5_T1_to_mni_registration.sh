for i in $@
do
    #echo ${i}
    rawT1=${i}/FREESURFER/mri/brain.nii.gz
    rawT1_mask=${i}/FREESURFER/mri/brainmask.nii.gz
    rawT1_mask_mgz=${i}/FREESURFER/mri/brainmask.mgz
    sourceImg=${rawT1}
    targetImg=/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz
    targetImgMask=/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain_mask.nii.gz

    #make T1_brain.nii.gz
    if [ ! -e ${rawT1_mask} ]
    then
        mri_convert --out_orientation RAS \
            ${rawT1_mask_mgz} \
            ${rawT1_mask}
        fslmaths ${rawT1_mask} -bin ${rawT1_mask}
    fi

    flirtMat=${i}/registration/T1_to_MNI.mat
    #flirt & fnirt
    if [ ! -e ${flirtMat} ]
    then
        flirt -in ${sourceImg} \
            -ref ${targetImg} \
            -out ${i}/registration/T1_to_MNI.nii.gz \
            -interp nearestneighbour \
            -omat ${flirtMat}
    fi

    if [ ! -e ${i}/registration/T1_to_MNI_fnirt_coeff.nii.gz ]
    then
        fnirt \
            --in=${sourceImg} \
            --ref=${targetImg} \
            --aff=${flirtMat} \
            --inmask=${rawT1_mask} \
            --refmask=${targetImgMask} \
            --iout=${i}/registration/T1_to_MNI_fnirt.nii.gz \
            --cout=${i}/registration/T1_to_MNI_fnirt_coeff
    fi

    for side in left right
    do
        if [ ! -e ${i}/segmentation/${side}/biggest_mni.nii.gz ]
        then
        applywarp \
            --ref=${targetImg} \
            --in=${i}/segmentation/${side}/biggest.nii.gz \
            --warp=${i}/registration/T1_to_MNI_fnirt_coeff \
            --out=${i}/segmentation/${side}/biggest_mni.nii.gz \
            --interp=nn

        fi
    done
done



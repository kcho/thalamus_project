for i in $@
do
    echo ${i}
    #echo ${i}
    fsDir=${i}/FREESURFER
    rawT1=${fsDir}/mri/brain.nii.gz
    rawT1_mask=${fsDir}/mri/brainmask.nii.gz
    rawT1_mask_mgz=${fsDir}/mri/brainmask.mgz

    dkiDir=${i}/DKI
    dkib0=${dkiDir}/b0.nii.gz
    dkib0_brain=${dkiDir}/b0_brain.nii.gz
    dkib0_brain_mask=${dkiDir}/b0_brain_mask.nii.gz

    dki_mask_from_fs=${dkiDir}/b0_mask_from_fs.nii.gz
    mk=${dkiDir}/kmean.nii
    mk_brain=${dkiDir}/kmean_brain.nii.gz

    mk_in_T1=${dkiDir}/kmean_T1.nii.gz
    mk_in_MNI=${dkiDir}/kmean_MNI.nii.gz
    mk_in_MNI_smooth_10=${dkiDir}/kmean_MNI_10s.nii.gz
    mk_in_MNI_smooth_5=${dkiDir}/kmean_MNI_5s.nii.gz
    
    regDir=${i}/registration
    T1_to_DKI_flirt_mat=${regDir}/FREESURFERT1toDKINodif.mat
    DKI_to_T1_flirt_mat=${regDir}/DKINodiftoFREESURFERT1.mat
    T1_to_MNI_fnirt_coeff=${regDir}/T1_to_MNI_fnirt_coeff


    targetImg=/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz
    targetImgMask=/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain_mask.nii.gz

    # make DKI_to_T1_flirt_mat
    if [ ! -e ${DKI_to_T1_flirt_mat} ]
    then
        convert_xfm -omat ${DKI_to_T1_flirt_mat} -inverse ${T1_to_DKI_flirt_mat}
    fi

    flirtMat=${i}/registration/T1_to_MNI.mat

    # bet MK using T1 mask from FREESURFER
    if [ ! -e ${dki_mask_from_fs} ]
    then
        flirt -in ${rawT1_mask} -ref ${dkib0} -applyxfm -init ${T1_to_DKI_flirt_mat} -out ${dki_mask_from_fs}
    fi

    if [ -e ${mk_brain} ]
    then
        fslmaths ${mk} -mas ${dki_mask_from_fs} ${mk_brain}
    fi

    echo ha
    # DKI to T1 flirt
    if [ -e ${mk_in_T1} ]
    then
        flirt -in ${mk_brain} -ref ${rawT1} -applyxfm -init ${DKI_to_T1_flirt_mat} -out ${mk_in_T1}
    fi

    # MK_T1 to MNI fnirt
    if [ -e ${mk_in_MNI} ]
    then
        applywarp \
            --ref=${targetImg} \
            --in=${mk_in_T1} \
            --warp=${T1_to_MNI_fnirt_coeff} \
            --out=${mk_in_MNI}
    fi


    # Smoothing
    if [ -e ${mk_in_MNI_smooth_5} ]
    then
        fslmaths ${mk_in_MNI} -s 5 ${mk_in_MNI_smooth_5}
    fi
done

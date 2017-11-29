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
regDir=${subj}/registration
dtiDir=${subj}/DTI
dkiDir=${subj}/DKI
roiDir=${subj}/ROI
tractDir=${subj}/thalamus_tractography_MNI
#tractDirMNI=${subj}/thalamus_tractography_MNI
segDir=${subj}/segmentation
bedpostDir=${subj}/DTI.bedpostX

for side in left right
do
    reconImg=${tractDir}/${side}/fdt_matrix2_reconstructed.nii.gz
    reconImgMNI=${tractDir}/${side}/fdt_matrix2_reconstructed_MNI.nii.gz
    reconImgMNI_ds=${tractDir}/${side}/fdt_matrix2_reconstructed_MNI_ds.nii.gz

    # probtracks postprocessing
    #if [ ! -e ${reconImg} ]
    #then 
        echo "Convert ${tractDir}/${side}/fdt_matrix2 --> ${i}"
        python tracktography/postprocessing/probtrackx_postprocessing.py \
            ${tractDir}/${side} \
            ${bedpostDir}/nodif_brain_mask.nii.gz
    #fi

    # do it in melodic
    #if [ ! -e ${reconImgMNI} ]
    #then
        echo 'Registration'
        flirt \
            -in ${reconImg} \
            -ref ${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz \
            -applyxfm -init ${regDir}/nodifToMNI.mat \
            -out ${reconImgMNI}
    #fi

    #if [ ! -e ${reconImgMNI_ds} ]
    #then
        echo 'Downsampling'
        flirt \
            -in ${reconImgMNI} \
            -ref ${reconImgMNI} \
            -applyisoxfm 4 \
            -out ${reconImgMNI_ds}
    #fi
done

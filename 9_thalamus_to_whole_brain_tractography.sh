#!/bin/sh
Usage() {
    echo ""
    echo "Usage: 9_thalamus_to_whole_brain_tractography.sh <subjDir> <lh>"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

# edit here for different folder structure
subj=${1}
side_s=${2}
fsDir=${subj}/FREESURFER
regDir=${subj}/registration
dtiDir=${subj}/DTI
roiDir=${subj}/ROI
tractDir=${subj}/thalamus_tractography
tractDir_MNI=${subj}/thalamus_tractography_MNI
bedpostDir=${subj}/DTI.bedpostX

if [ ${side_s} == 'lh' ]
then
    side=left
elif [ ${side_s} == 'rh' ]
then
    side=right
else
    Usage
fi

echo ${subj} ${side}

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
# fsThalROI=${roiDir}/${side_s}_thalamus.nii.gz 
# dtiThalROI=${roiDir}/${side_s}_thalamus_DTI.nii.gz 
mniThalROI_raw=${side_s}_thalamus_HOSC_60.nii.gz
mniThalROI=${roiDir}/${side_s}_thalamus_DTI_HO.nii.gz 

if [ ! -e ${mniThalROI_raw} ]
then
    fslroi ${FSLDIR}/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz lh_thalamus_HOSC.nii.gz 3 1
    fslmaths lh_thalalmus_HOSC.nii.gz -thr 60 -bin lh_thalamus_HOSC_60.nii.gz

    fslroi ${FSLDIR}/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz rh_thalamus_HOSC.nii.gz 14 1
    fslmaths rh_thalalmus_HOSC.nii.gz -thr 60 -bin rh_thalamus_HOSC_60.nii.gz
fi

if [ ! -e ${mniThalROI} ]
then 
    flirt \
        -in ${mniThalROI_raw} \
        -ref ${bedpostDir}/nodif_brain_mask.nii.gz \
        -applyxfm -init ${regDir}/MNItoNodif.mat \
        -interp nearestneighbour \
        -out ${mniThalROI}
fi

# Tractography using fnirt
mni2fs2nodif=${regDir}/MNI_to_FREESURFER_fnirt_to_Nodif_flirt.nii.gz
nodif2fs2mni=${regDir}/Nodif_to_FREESURFER_flirt_to_MNI_fnirt.nii.gz
if [ ! -e ${tractDir_MNI}/${side}/fdt_paths.nii.gz ]
then
    rm -rf ${tractDir_MNI}/${side} 
    mkdir -p ${tractDir_MNI}/${side}
    probtrackx2 \
        -x ${mniThalROI_raw} \
        -l \
        --onewaycondition \
        --omatrix2 \
        --target2=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz \
        -c 0.2 \
        -S 2000 \
        --steplength=0.5 \
        -P 5000 \
        --fibthresh=0.01 \
        --distthresh=0.0 \
        --sampvox=0.0 \
        --xfm=${mni2fs2nodif} \
        --invxfm=${nodif2fs2mni} \
        --forcedir \
        --opd \
        -s ${bedpostDir}/merged \
        -m ${bedpostDir}/nodif_brain_mask \
        --dir=${tractDir_MNI}/${side}
    echo ${subj} MNI thalamo-whole brain tractography on the ${side} done
else
    echo ${subj} MNI thalamo-whole brain tractography on the ${side} done
fi

reconImg=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed.nii.gz
reconImg4s=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed_4s.nii.gz
reconImg_ds=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed_ds.nii.gz

# probtracks postprocessing
if [ ! -e ${reconImg} ]
then 
    echo "Convert ${tractDir_MNI}/${side}/fdt_matrix2 --> ${i}"
    python tracktography/postprocessing/probtrackx_postprocessing.py \
        -i ${tractDir_MNI}/${side}
        #${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz
fi

## smoothing
fslmaths ${reconImg} -kernel gauss 1.69865806 -fmean ${reconImg4s}

#if [ ! -e ${reconImgMNI_ds} ]
#then
#echo 'Downsampling'
flirt \
    -in ${reconImg} \
    -ref ${reconImg} \
    -applyisoxfm 4 \
    -out ${reconImg_ds}
#fi

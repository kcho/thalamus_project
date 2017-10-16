#!/bin/sh
Usage() {
    echo ""
    echo "Usage: 9_thalamus_to_whole_brain_tractography.sh <subjDir> <lh>"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

################################
# Input from commandline
################################
subj=${1}
side_s=${2}

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

# Edit here for different folder structure

fsDir=${subj}/FREESURFER
regDir=${subj}/registration
dtiDir=${subj}/DTI
roiDir=${subj}/ROI
tractDir=${subj}/thalamus_tractography
tractDir_MNI=${subj}/wb_thalamus_tractography
bedpostDir=${subj}/DTI.bedpostX

rawT1=${fsDir}/mri/brain.nii.gz
rawT1_mask=${fsDir}/mri/brainmask.nii.gz
rawT1_mask_mgz=${fsDir}/mri/brainmask.mgz
nodif_brain=${dtiDir}/nodif_brain.nii.gz

#sourceImg=${rawT1}
mni=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
mniMask=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz

################################################################
# 0. nifti from freesurfer
################################################################
if [ ! -e ${rawT1_mask} ]
then
    mri_convert \
        --out_orientation RAS \
        ${rawT1_mask_mgz} \
        ${rawT1_mask}
    fslmaths ${rawT1_mask} -bin ${rawT1_mask}
fi

################################################################
# 1.1 Registration : fs --> mni fnirt
################################################################
fs2mni_flirt=${regDir}/fs2mni.mat
#flirt & fnirt
if [ ! -e ${fs2mni_flirt} ]
then
    flirt \
        -in ${rawT1} \
        -ref ${mni} \
        -omat ${fs2mni_flirt}
fi

fs2mni_fnirt=${regDir}/fs2mni_fnirt_coeff.nii.gz
fs2mni_fnirt_img=${regDir}/fs2mni_fnirt_img.nii.gz
if [ ! -e ${fs2mni_fnirt} ]
then
    fnirt \
        --in=${rawT1} \
        --ref=${mni} \
        --aff=${fs2mni_flirt} \
        --inmask=${rawT1_mask} \
        --refmask=${mniMask} \
        --cout=${fs2mni_fnirt} \
        --iout=${fs2mni_fnirt_img} 
fi

mni2fs_fnirt=${regDir}/mni2fs_fnirt_coeff.nii.gz
if [ ! -e ${mni2fs_fnirt} ]
then
    invwarp \
        -w ${fs2mni_fnirt} \
        -o ${mni2fs_fnirt} \
        -r ${rawT1}
fi

################################################################
# 1.2 Registration : fs --> DTI flirt
################################################################
fs2nodif=${regDir}/fs2nodif.mat
if [ ! -e ${fs2nodif} ]
then
    flirt \
        -in ${rawT1} \
        -ref ${nodif_brain} \
        -omat ${fs2nodif}
fi

################################################################
# 1.3 Registration : MNI --> fs --> DTI 
################################################################
mni2fs2nodif=${regDir}/mni2fs2nodif_coeff.nii.gz
if [ ! -e ${mni2fs2nodif} ]
then
    convertwarp \
        --ref=${nodif_brain} \
        --warp1=${mni2fs_fnirt} \
        --postmat=${fs2nodif} \
        --out=${mni2fs2nodif}

fi

nodif2fs2mni=${regDir}/nodif2fs2mni_coeff.nii.gz
if [ ! -e ${nodif2fs2mni} ]
then
    invwarp \
        -w ${mni2fs2nodif} \
        -o ${nodif2fs2mni} \
        -r ${mni}
fi


################################################################
# 2. Extract thalamic ROI from the Harvard Oxford template
################################################################
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
    applywarp \
        --ref=${nodif_brain} \
        --in=${mniThalROI_raw} \
        --warp=${mni2fs2nodif} \
        --out=${mniThalROI} \
        --interp=nn
fi

################################################################
# 3. Thalamus seeded whole-brain tractography
################################################################
if [ ! -e ${tractDir_MNI}/${side}/fdt_paths.nii.gz ]
then
    rm -rf ${tractDir_MNI}/${side} 
    mkdir -p ${tractDir_MNI}/${side}
    probtrackx2 \
        -x ${mniThalROI_raw} \
        -l \
        --onewaycondition \
        --omatrix2 \
        --target2=${mniMask} \
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

################################################################
# 4. Downsampling --> requires editing 
################################################################
reconImg=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed.nii.gz
reconImg_ds_3=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed_ds_3.nii.gz
reconImg_ds_4=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed_ds_4.nii.gz
reconImg4s=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed_4s.nii.gz
reconImg_ds_3_4s=${tractDir_MNI}/${side}/fdt_matrix2_reconstructed_ds_3_4s.nii.gz

# probtracks postprocessing
#if [ ! -e ${reconImg} ]
#then 
    #echo "Convert ${tractDir_MNI}/${side}/fdt_matrix2 --> ${i}"
    #python tracktography/postprocessing/probtrackx_postprocessing.py \
        #-i ${tractDir_MNI}/${side}
        ##${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz
#fi

## smoothing
fslmaths ${reconImg} -kernel gauss 1.69865806 -fmean ${reconImg4s}

#if [ ! -e ${reconImgMNI_ds} ]
#then
#echo 'Downsampling'
#flirt \
    #-in ${reconImg} \
    #-ref ${reconImg} \
    #-applyisoxfm 4 \
    #-out ${reconImg_ds_3}
if [ ! -e ${reconImg_ds_4} ]
then
    flirt \
        -in ${reconImg} \
        -ref ${reconImg} \
        -applyisoxfm 3 \
        -out ${reconImg_ds_4}
fi

flirt \
    -in ${reconImg4s} \
    -ref ${reconImg4s} \
    -applyisoxfm 3 \
    -out ${reconImg_ds_3_4s}
#fi

import os
from os import join, basename, dirname, isfile, isdir
from glob import glob
import re
from nipype.interfaces import fsl
from nipype.interfaces.freesurfer import ReconAll, MRIConvert
import sys

def run_reconall(subjDir, t1_location):
    res = ReconAll(subject_id = 'FREESURFER',
                   directive = 'all',
                   subjects_dir = subjDir,
                   T1_files = t1_location).run()
    reconall.run()

def mgz_to_nifti(mgz_location, Mask=False):
    outfile = re.sub('.mgz$', '.nii.gz$', mgz_location)

    res = MRIConvert(in_file = mgz_location,
                           out_file = outfile,
                           out_type = 'nii.gz',
                           out_orientation = 'RAS').run()

    if Mask:
        command = 'fslmaths {imgInput} -bin {imgOutput}'.format(
            imgInput = outfile,
            imgOutput = outfile)
        os.popen(command).read()

def flirt_within_subject(img, ref, omat):
    flt = fsl.FLIRT(bins=256, 
                    cost_func='mutualinfo',
                    searchrx=[-180, 180],
                    searchry=[-180, 180],
                    searchrz=[-180, 180],
                    dof=6,
                    interp='trilinear',
                    in_file=img,
                    reference=ref,
                    out_matrix_file=omat)
    flt.run()
    
def fnirt(img, ref, aff, out, inmask, refmask, mni=False):
    frt = fsl.FNIRT(in_file=img,
                    ref_file=ref,
                    affine_file=aff, 
                    inmask_file=inmask, 
                    refmask_file=refmask,
                    fieldcoeff_file=out)

    if mni:
        frt.inputs.config_file='T1_2_MNI152_2mm'

    frt.run()

def preprocessing(subjDir):
    raw_t1_location = join(subjDir, 'T1', 'T1.nii.gz')
    run_reconall(subjDir, raw_t1_location)
    
    mgz_imgs_to_convert = ['brain', 'T1', 'brainmask']
    mgz_imgs_loc = [join(subjDir, 'FREESURFER', 'mri', x+'.mgz') for x in mgz_imgs_to_convert]
    for mgz_img in mgz_imgs_loc:
        if 'mask' in mgz_img:
            mgz_to_nifti(mgz_img, Mask=True)
        else:
            mgz_to_nifti(mgz_img)


    regDir = join(subjDir, 'registration')
    if not isdir(regDir):
        os.mkdir(regDir)

    # DTI / DKI flirt
    fsDir = join(subjDir, 'FREESURFER', 'mri')
    fs_t1 = join(fsDir, 'T1.nii.gz')
    fs_t1_mask = join(fsDir, 'brainmask.nii.gz')
    dtiDir = join(subjDir, 'DTI')
    dkiDir = join(subjDir, 'DKI')
    dti_nodif_brain = join(dtiDir, 'nodif_brain.nii.gz')
    dki_nodif_brain = join(dkiDir, 'nodif_brain.nii.gz')

    fs2dti=join(regDir, 'fs2dti.mat')
    fs2dki=join(regDir, 'fs2dki.mat')

    flirt_within_subject(dti_nodif_brain, fs_t1, fs2dti)
    flirt_within_subject(dki_nodif_brain, fs_t1, fs2dki)

    fs2dti_fnirt=join(regDir, 'fs2dti_coeff.nii.gz')
    fs2dki_fnirt=join(regDir, 'fs2dki_coeff.nii.gz')
    fnirt(dti_nodif_brain, fs_t1, fs2dti, d


    fs2dki_fnirt_inv=join(regDir, 'fs2dki_coeff_inv.nii.gz')

    

    
def fnirt(img, ref, aff, out, inmask, refmask, mni=False):
    invwarp = InvWarp(warp=fnirt_out, 
                      reference=ref,
                     inverse_warp=out)
    invwarp.run()


    convwarp = ConvertWarp(warp1=inv_out,
                           reference=ref,
                           postmat=postmat,
                           out_file=out)
    convwarp.run()

    applywarp = ApplyWarp(in_file=in_img,
                          ref_file=ref,
                          field_file=warp,
                          out_file=out)
    applywarp.run()

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

if __name__=='__main__':
    preprocessing(sys.argv[1])

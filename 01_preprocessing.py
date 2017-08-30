import os
from os.path import join, basename, dirname, isfile, isdir
from glob import glob
import re
from nipype.interfaces import fsl
from nipype.interfaces.freesurfer import ReconAll, MRIConvert
import sys

class subject:
    def get_img_raw(self, t1Dir):
        file_list = os.listdir(t1Dir)
        co_list = [x for x in file_list if re.search('^20.*\.nii\.gz$', x)]
        if len(co_list)==1:
            return co_list[0]
        else:
            return 'No T1 images'

    def run_recon_all(self):
        if not isdir(self.fsDir):
            ReconAll(subject_id = 'FREESURFER',
                     directive = 'all',
                     subjects_dir = subjDir,
                     T1_files = self.t1raw).run()

    def convert_mgzs_into_niftis(self):
        if not isfile(self.fs_t1):
            MRIConvert(in_file=self.fs_t1_mgz,
                       out_file=self.fs_t1,
                       out_type='nii.gz',
                       out_orientation='RAS').run()

        if not isfile(self.fs_t1_brain):
            MRIConvert(in_file=self.fs_t1_brain_mgz,
                       out_file=self.fs_t1_brain,
                       out_type='nii.gz',
                       out_orientation='RAS').run()

        if not isfile(self.fs_t1_brain_mask):
            MRIConvert(in_file=self.fs_t1_brain_mask_mgz,
                                  out_file=self.fs_t1_brain_mask,
                                  out_type='nii.gz',
                                  out_orientation='RAS').run()
            command = 'fslmaths {imgInput} -bin {imgOutput}'.format(
                imgInput = self.fs_t1_brain_mask,
                imgOutput = self.fs_t1_brain_mask)
            os.popen(command).read()

    def dti_preproc(self):
        if not isfile(self.dti_eddy_out):
            fsl.EddyCorrect(in_file=self.dtiraw, 
                            out_file=self.dti_eddy_out).run()

        if not isfile(self.dti_nodif):
            fsl.ExtractROI(in_file=self.dti_eddy_out,
                           t_min=self.nodif_number,
                           t_size=1,
                           roi_file=self.dti_nodif).run()

        if not isfile(self.nodif_brain):
            fsl.Bet(in_file = self.nodif,
                    frac = 0.3,
                    mask = True,
                    out_file = self.nodif_brain).run()
    def dki_preproc(self):
        if not isfile(self.dki_eddy_out):
            fsl.EddyCorrect(in_file=self.dkiraw, 
                            out_file=self.dki_eddy_out).run()

        if not isfile(self.dki_nodif):
            fsl.ExtractROI(in_file=self.dki_eddy_out,
                           t_min=self.nodif_number,
                           t_size=1,
                           roi_file=self.dki_nodif).run()

        if not isfile(self.nodif_brain):
            fsl.Bet(in_file = self.nodif,
                    frac = 0.3,
                    mask = True,
                    out_file = self.nodif_brain).run()
    def fs_to_dti_flirt(self):
        if not isfile(self.fs2dti_mat):
            fsl.FLIRT(bins=256, 
                      cost_func='mutualinfo',
                      searchrx=[-180, 180],
                      searchry=[-180, 180],
                      searchrz=[-180, 180],
                      dof=6,
                      interp='trilinear',
                      in_file=self.fs_t1_brain,
                      reference=self.dti_nodif_brain,
                      out_matrix_file=self.fs2dti_mat).run()

    def fs_to_dki_flirt(self):
        if not isfile(self.fs2dki_mat):
            fsl.FLIRT(bins=256, 
                      cost_func='mutualinfo',
                      searchrx=[-180, 180],
                      searchry=[-180, 180],
                      searchrz=[-180, 180],
                      dof=6,
                      interp='trilinear',
                      in_file=self.fs_t1_brain,
                      reference=self.dki_nodif_brain,
                      out_matrix_file=self.fs2dki_mat).run()
    def fs_to_dti_fnirt(self):
        if not isfile(self.fs2dti_coeff):
            fsl.FNIRT(in_file=self.fs_t1,
                      ref_file=self.dti_nodif,
                      affine_file=self.fs2dti_mat, 
                      inmask_file=self.fs_t1_brain_mask, 
                      refmask_file=self.dti_nodif_brain_mask,
                      fieldcoeff_file=self.fs2dti_coeff).run()

    def fs_to_dki_fnirt(self):
        if not isfile(self.fs2dki_coeff):
            fsl.FNIRT(in_file=self.fs_t1,
                      ref_file=self.dki_nodif,
                      affine_file=self.fs2dki_mat, 
                      inmask_file=self.fs_t1_brain_mask, 
                      refmask_file=self.dki_nodif_brain_mask,
                      fieldcoeff_file=self.fs2dki_coeff).run()

    def fs_to_dti_fnirt(self):
        if not isfile(self.dti2fs_coeff):
            fsl.ConvertWarp(warp1=self.fs2dti_coeff,
                            reference=self.dti_nodif,
                            postmat=self.fs2dti_mat,
                            out_file=self.dti2fs_coeff).run()

    def fs_to_dki_fnirt(self):
        if not isfile(self.dki2fs_coeff):
            fsl.ConvertWarp(warp1=self.fs2dki_coeff,
                            reference=self.dki_nodif,
                            postmat=self.fs2dki_mat,
                            out_file=self.dki2fs_coeff).run()

    def mni_to_fs_flirt(self):
        if not isfile(self.mni2fs_mat):
            fsl.FLIRT(bins=256, cost_func='mutualinfo',
                      searchrx=[-180, 180],
                      searchry=[-180, 180],
                      searchrz=[-180, 180], dof=12,
                      interp='trilinear',
                      in_file=self.mni_brain,
                      reference=self.fs_t1_brain,
                      out_matrix_file=self.mni2fs_mat).run()
    def mni_to_fs_fnirt(self):
        if not isfile(self.mni2fs_coeff):
            fsl.FNIRT(in_file=self.mni, 
                      ref_file=self.fs_t1,
                      affine_file=self.mni2fs_mat, 
                      inmask_file=self.mni_brain_mask, 
                      config_file='T1_2_MNI152_2mm',
                      refmask_file=self.fs_t1_brain_mask, 
                      fieldcoeff_file=self.mni2fs_coeff).run()

    def mni_to_fs_to_dti_convwarp(self):
        if not isfile(self.mni2fs2dti):
            fsl.ConvertWarp(warp1=self.mni2fs_coeff,
                            reference=self.dti_nodif,
                            postmat=self.fs2dti_mat,
                            out_file=self.mni2fs2dti).run()

    def dti_to_fs_to_mni_from_invwarp(self):
        if not isfile(self.dti2fs2mni):
            fsl.InvWarp(warp=self.mni2fs2dti,
                        reference=self.mni,
                        inverse_warp=self.dti2fs2mni).run()

    def visual_check_nodif_to_mni(self):
        if not isfile(self.nodif2mni_check):
            fsl.ApplyWarp(in_file=self.dti_nodif,
                          ref_file=self.mni,
                          out_file=self.dti_nodif2mni_check,
                          field_file=self.dti2fs2mni).run()

    def __init__(self, subjDir):
        fsldir = os.environ['FSLDIR']
        fslstandard = join(fsldir, 'data', 'standard')
        self.mni = join(fslstandard, 'MNI152_T1_2mm.nii.gz')
        self.mni_brain = join(fslstandard, 'MNI152_T1_2mm_brain.nii.gz')
        self.mni_brain_mask = join(fslstandard, 'MNI152_T1_2mm_brain_mask.nii.gz')

        self.dir = subjDir
        self.t1Dir = join(subjDir, 'T1')
        self.t1raw = self.get_img_raw(self.t1Dir)
        self.fsDir = join(subjDir, 'FREESURFER')

        self.fsMriDir = join(self.fsDir, 'mri')
        self.fs_t1_mgz = join(self.fsMriDir, 'T1.mgz')
        self.fs_t1 = join(self.fsMriDir, 'T1.nii.gz')

        self.fs_t1_brain_mgz = join(self.fsMriDir, 'brain.mgz')
        self.fs_t1_brain = join(self.fsMriDir, 'brain.nii.gz')

        self.fs_t1_brain_mask_mgz = join(self.fsMriDir, 'brainmask.mgz')
        self.fs_t1_brain_mask = join(self.fsMriDir, 'brainmask.nii.gz')
        self.dtiDir = join(subjDir, 'DTI')
        self.nodif_number = 0
        self.dtiraw = self.get_img_raw(self.dtiDir)
        self.dti_eddy_out = join(self.dtiDir, 'data.nii.gz')
        self.dti_nodif = join(self.dtiDir, 'nodif.nii.gz')
        self.dti_nodif_brain = join(self.dtiDir, 'nodif_brain.nii.gz')
        self.dti_nodif_brain_mask = join(self.dtiDir, 'nodif_brain_mask.nii.gz')


        self.dkiDir = join(subjDir, 'DKI')
        self.nodif_number = 0
        self.dkiraw = self.get_img_raw(self.dkiDir)
        self.dki_eddy_out = join(self.dkiDir, 'data.nii.gz')
        self.dki_nodif = join(self.dkiDir, 'nodif.nii.gz')
        self.dki_nodif_brain = join(self.dkiDir, 'nodif_brain.nii.gz')
        self.dki_nodif_brain_mask = join(self.dkiDir, 'nodif_brain_mask.nii.gz')


        self.regDir = join(subjDir, 'registration')
        if not isdir(self.regDir):
            os.mkdir(self.regDir)

        self.fs2dti_mat = join(self.regDir, 'fs2dti.mat')
        self.fs2dki_mat = join(self.regDir, 'fs2dki.mat')

        self.fs2dti_coeff = join(self.regDir, 'fs2dti_coeff.nii.gz')
        self.fs2dki_coeff = join(self.regDir, 'fs2dki_coeff.nii.gz')
        self.dti2fs_coeff = join(self.regDir, 'dti2fs_coeff.nii.gz')
        self.dki2fs_coeff = join(self.regDir, 'dki2fs_coeff.nii.gz')

        #self.mni2dti_mat = join(self.regDir, 'mni2dti.mat')
        #self.mni2dki_mat = join(self.regDir, 'mni2dki.mat')
        #if not isfile(self.mni2dti_mat):
            #fsl.FLIRT(bins=256, 
                      #cost_func='mutualinfo', 
                      #searchrx=[-180, 180],
                      #searchry=[-180, 180],
                      #searchrz=[-180, 180],
                      #dof=9,
                      #interp='trilinear',
                      #in_file=self.mni_brain,
                      #reference=self.dti_nodif_brain,
                      #out_matrix_file=self.mni2dti_mat).run()
        #if not isfile(self.mni2dki_mat):
            #fsl.FLIRT(bins=256, 
                      #cost_func='mutualinfo', 
                      #searchrx=[-180, 180],
                      #searchry=[-180, 180],
                      #searchrz=[-180, 180],
                      #dof=9,
                      #interp='trilinear',
                      #in_file=self.mni_brain,
                      #reference=self.dki_nodif_brain,
                      #out_matrix_file=self.mni2dki_mat).run()

        #self.mni2dti_coeff = join(self.regDir, 'mni2dti_coeff.nii.gz')
        #if not isfile(self.mni2dti_coeff):
            #fsl.FNIRT(in_file=self.mni, 
                      #ref_file=self.dti_nodif,
                      #affine_file=self.mni2dti_mat, 
                      #inmask_file=self.mni_brain_mask, 
                      #config_file='T1_2_MNI152_2mm',
                      #refmask_file=self.dti_nodif_brain_mask, 
                      #fieldcoeff_file=self.mni2dti_coeff).run()
        #self.mni2dki_coeff = join(self.regDir, 'mni2dki_coeff.nii.gz')
        #if not isfile(self.mni2dki_coeff):
            #fsl.FNIRT(in_file=self.mni, 
                      #ref_file=self.dki_nodif,
                      #affine_file=self.mni2dki_mat, 
                      #inmask_file=self.mni_brain_mask, 
                      #config_file='T1_2_MNI152_2mm',
                      #refmask_file=self.dki_nodif_brain_mask, 
                      #fieldcoeff_file=self.mni2dki_coeff).run()

        #self.dti2mni_coeff = join(self.regDir, 'dti2mni_coeff.nii.gz')
        #if not isfile(self.dti2mni_coeff):
            #fsl.ConvertWarp(warp1=self.mni2dti_coeff,
                            #reference=self.dti_nodif,
                            #postmat=self.mni2dti_mat,
                            #out_file=self.dti2mni_coeff)

        #self.dki2mni_coeff = join(self.regDir, 'dki2mni_coeff.nii.gz')
        #if not isfile(self.dki2mni_coeff):
            #fsl.ConvertWarp(warp1=self.mni2dki_coeff,
                            #reference=self.dki_nodif,
                            #postmat=self.mni2dki_mat,
                            #out_file=self.dki2mni_coeff)

        self.mni2fs_mat = join(self.regDir, 'mni2fs.mat')
        self.mni2fs_coeff = join(self.regDir, 'mni2fs_coeff.nii.gz')
        self.mni2fs2dti = join(self.regDir, 'mni2fs2dti_coeff.nii.gz')
        self.dti2fs2mni = join(self.regDir, 'dti2fs2mni_coeff.nii.gz')

        self.dti_nodif2mni_check = join(self.regDir,
                                    'nodif2mni_check.nii.gz')

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

#for subj in $@
#do
    ## Path variables
    #t1Dir=${subj}/T1

    ## Freesurfer
    #fsDir=${subj}/FREESURFER
    #fsDirName=${fsDir##*/}
    #fsMriDir=${fsDir}/mri
    #fsT1Img=${fsMriDir}/T1.nii.gz
    #fsT1BetImg=${fsMriDir}/brain.nii.gz
    #fsT1MaskImg=${fsMriDir}/brainmask.nii.gz
    ## MNI
    #mni2mmImg=${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz
    #mni2mmBetImg=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
    #mni2mmMaskImg=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz 
    ## DTI
    #dtiDir=${subj}/DTI
    #dtiNodifBetImg=${dtiDir}/nodif_brain.nii.gz
    #dtiNodifImg=${dtiDir}/nodif.nii.gz
    #dtiNodifMaskImg=${dtiDir}/nodif_brain_mask.nii.gz 
    ## DKI
    #dkiDir=${subj}/DKI
    #dkiNodifBetImg=${dkiDir}/nodif_brain.nii.gz
    #dkiNodifImg=${dkiDir}/nodif.nii.gz
    #dkiNodifMaskImg=${dkiDir}/nodif_brain_mask.nii.gz
    ## other dirs
    #regDir=${subj}/registration
    #roiDir=${subj}/ROI
    #bedpostDir=${subj}/DTI.bedpostX


    ## freesurfer mgz --> nii.gz
    #if [ ! -e ${fsT1BetImg} ]
    #then
        #mri_convert \
            #--out_orientation RAS \
            #${fsMriDir}/brain.mgz \
            #${fsT1BetImg}
    #fi

    #if [ ! -e ${fsT1Img} ]
    #then
        #mri_convert \
            #--out_orientation RAS \
            #${fsMriDir}/T1.mgz \
            #${fsT1Img}
    #fi

    #if [ ! -e ${fsT1MaskImg} ]
    #then
        #mri_convert \
            #--out_orientation RAS \
            #${fsMriDir}/brainmask.mgz \
            #${fsT1MaskImg}
        #fslmaths \
            #${fsT1MaskImg} \
            #-bin ${fsT1MaskImg}
    #fi

    ## Registration
    #if [ ! -d ${regDir} ]
    #then
        #mkdir ${regDir}
    #fi

    ## flirt : freesurfer --> DTI/nodif_brain.nii.gz
    #fs2dti=${regDir}/fs2dti.mat
    #if [ ! -e ${fs2dti} ]
    #then
        #flirt \
            #-in ${fsT1BetImg} \
            #-ref ${dtiNodifBetImg} \
            #-omat ${fs2dti} \
            #-bins 256 \
            #-cost mutualinfo \
            #-searchrx -180 180 \
            #-searchry -180 180 \
            #-searchrz -180 180 \
            #-dof 6  -interp trilinear
    #fi

    ## fnirt : freesurfer nobet T1 --> DTI/nodif.nii.gz
    #fs2dti_fnirt=${regDir}/fs2dti_fnirt_coeff.nii.gz
    #fs2dti_fnirt_inv=${regDir}/fs2dti_fnirt_coeff_inv.nii.gz
    #if [ ! -e ${fs2dti_fnirt} ]
    #then
        #fnirt \
            #--ref=${dtiNodifImg} \
            #--in=${fsT1Img} \
            #--aff=${fs2dti} \
            #--inmask=${fsT1MaskImg} \
            #--refmask=${dtiNodifMaskImg} \
            #--cout=${fs2dti_fnirt}
    #fi

    #if [ ! -e ${fs2dti_fnirt_inv} ]
    #then
        #invwarp \
            #-w ${fs2dti_fnirt} \
            #-o ${fs2dti_fnirt_inv} \
            #-r ${fsT1Img}
    #fi

    ## flirt : freesurfer --> DKI/nodif_brain.nii.gz
    #fs2dki=${regDir}/fs2dki.mat
    #if [ ! -e ${fs2dki} ]
    #then
        #flirt \
            #-in ${fsT1BetImg} \
            #-ref ${dkiNodifBetImg} \
            #-omat ${fs2dki} \
            #-bins 256 \
            #-cost mutualinfo \
            #-searchrx -180 180 \
            #-searchry -180 180 \
            #-searchrz -180 180 \
            #-dof 6  -interp trilinear
    #fi

    ## FREESURFER output t1 --> diffusion kurtosis space fnirt
    #fs2dki_fnirt=${regDir}/fs2dki_fnirt_coeff.nii.gz
    #fs2dki_fnirt_inv=${regDir}/fs2dki_fnirt_coeff_inv.nii.gz
    #if [ ! -e ${fs2dki_fnirt} ]
    #then
        #fnirt \
            #--ref=${dkiNodifImg} \
            #--in=${fsT1Img} \
            #--aff=${fs2dki} \
            #--inmask=${fsT1MaskImg} \
            #--refmask=${dkiNodifMaskImg} \
            #--cout=${fs2dki_fnirt}
    #fi

    #if [ ! -e ${fs2dki_fnirt_inv} ]
    #then
        #invwarp \
            #-w ${fs2dki_fnirt} \
            #-o ${fs2dki_fnirt_inv} \
            #-r ${fsT1Img}
    #fi

    ## MNI 2mm --> FREESURFER output T1 
    #mni2fs=${regDir}/mni2fs.mat
    #if [ ! -e ${mni2fs} ]
    #then
        #flirt \
            #-in ${mni2mmImg} \
            #-ref ${fsT1BetImg} \
            #-omat ${mni2fs} \
            #-bins 256 \
            #-cost mutualinfo \
            #-searchrx -180 180 \
            #-searchry -180 180 \
            #-searchrz -180 180 \
            #-dof 12 -interp trilinear
    #fi
    #mni2fs_fnirt=${regDir}/mni2fs_fnirt_coeff.nii.gz
    #if [ ! -e ${mni2fs_fnirt} ]
    #then
        #fnirt \
            #--ref=${fsT1BetImg} \
            #--in=${mni2mmImg} \
            #--config=T1_2_MNI152_2mm \
            #--aff=${mni2fs} \
            #--inmask=${mni2mmMaskImg} \
            #--refmask=${fsT1MaskImg} \
            #--cout=${mni2fs_fnirt}
    #fi
    ## MNI --fnirt--> Freesurfer --flirt--> nodif
    #mni2fs2nodif=${regDir}/mni2fs2nodif_coeff.nii.gz
    #if [ ! -e ${mni2fs2nodif} ]
    #then
        #convertwarp \
            #--ref=${dtiNodifMaskImg} \
            #--warp1=${mni2fs_fnirt} \
            #--postmat=${fs2nodifMat} \
            #--out=${mni2fs2nodif}
    #fi

    ## nodif --flirt--> Freesurfer --fnirt--> MNI
    #nodif2fs2mni=${regDir}/nodif2fs2mni_coeff.nii.gz
    #if [ ! -e ${nodif2fs2mni} ]
    #then
        #invwarp \
            #--ref=${mni2mmBetImg} \
            #--warp=${mni2fs2nodif} \
            #--out=${nodif2fs2mni}
    #fi

    ## Check ouput image
    #mni2nodif_img=${regDir}/mni2nodif_img.nii.gz
    #if [ ! -e ${mni2nodif_img} ]
    #then
        #applywarp \
            #--ref=${dtiNodifMaskImg} \
            #--in=${mni2mmBetImg} \
            #--warp=${mni2fs2nodif} \
            #--out=${mni2nodif_img}
    #fi

    ## ROI extraction
    #if [ ! -d ${roiDir} ] 
    #then
        #mkdir ${roiDir}
    #fi

    #if [ ! -e ${roiDir}/rh_OCC.nii.gz ]
    #then
        #python roiExtraction.py \
            #-s ${subj} \
            #-r `basename ${roiDir}` \
            #-f `basename ${fsDir}`
    #else
        #echo ${subj} 'ROI extraction done'
    #fi

    #if [ ! -d ${bedpostDir} ]
    #then
        #bedpostx ${dtiDir}
    #fi
#done

#if __name__=='__main__':
    #preprocessing(sys.argv[1])

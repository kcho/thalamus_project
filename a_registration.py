# -*- coding: utf-8 -*-
import sys
from thalamus_subject import *

class thalamus_subject(thalamus_subject):
    '''
    based on thalamus_subject class in thalamus_subject.py
    eg)
    f = subject(sys.argv[1])
    f.run_recon_all()
    f.convert_mgzs_into_niftis()
    f.dti_preproc()
    f.fs_to_dti_flirt()

    f.mni_to_fs_flirt()
    f.mni_to_fs_fnirt()

    f.mni_to_fs_to_dti_convwarp()
    f.dti_to_fs_to_mni_from_invwarp()
    f.visual_check_nodif_to_mni()
    
    
    Example output)
    subjectdir
      ├─ DTI
      │   ├─ 20*.nii.gz : raw dwi
      │   ├─ 20*.bvec : raw bvectors
      │   └─ 20*.bval : raw bvalues
      │ 
      ├─ T1
      │   └─ 20*.nii.gz : raw T1
      │
      ├─ DKI (optional)
      │   ├─ 20*.nii.gz : raw dwi
      │   ├─ 20*.bvec : raw bvectors
      │   └─ 20*.bval : raw bvalues
      │
      ├─ FREESURFER (optional)
      │   ├─ ...
      │   └─ mri
      │       ├─ T1.nii.gz : T1 image
      │       ├─ brain.nii.gz : T1_brain image
      │       └─ brainmask.nii.gz : T1_brain mask image
      │
      └─ registration
          ├─ fs2dti.mat : freesurfer to dti flirt
          ├─ mni2fs.mat : MNI to freesurfer flirt
          ├─ mni2fs_coeff.nii.gz : MNI to freesurfer fnirt
          ├─ mni2fs2dti_coeff.nii.gz : MNI to freesurfer to DTI fnirt
          ├─ dti2fs2mni_coeff.nii.gz : DTI to freesurfer to MNI frnit
          └─ nodif2mni_check.nii.gz : DTI image registered on MNI

    '''
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
        for img_mgz, img_nifti in zip([self.fs_t1_mgz,
                                       self.fs_t1_brain_mgz,
                                       self.fs_t1_brain_mask_mgz],
                                      [self.fs_t1,
                                       self.fs_t1_brain,
                                       self.fs_t1_brain_mask]):
            if not isfile(img_nifti):
                MRIConvert(in_file=img_mgz,
                           out_file=img_nifti,
                           out_type='niigz',
                           out_orientation='RAS').run()
            if 'mask' in img_nifti:
                command = 'fslmaths {imgInput} -bin {imgOutput}'.format(
                    imgInput = img_nifti,
                    imgOutput = img_nifti)
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

        if not isfile(self.dti_nodif_brain):
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
                      searchr_x=[-180, 180],
                      searchr_y=[-180, 180],
                      searchr_z=[-180, 180],
                      dof=6,
                      interp='trilinear',
                      in_file=self.fs_t1_brain,
                      reference=self.dti_nodif_brain,
                      out_matrix_file=self.fs2dti_mat).run()

    def fs_to_dki_flirt(self):
        if not isfile(self.fs2dki_mat):
            fsl.FLIRT(bins=256, 
                      cost_func='mutualinfo',
                      searchr_x=[-180, 180],
                      searchr_y=[-180, 180],
                      searchr_z=[-180, 180],
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
                      searchr_x=[-180, 180],
                      searchr_y=[-180, 180],
                      searchr_z=[-180, 180], dof=12,
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
        if not isfile(self.dti_nodif2mni_check):
            fsl.ApplyWarp(in_file=self.dti_nodif,
                          ref_file=self.mni,
                          out_file=self.dti_nodif2mni_check,
                          field_file=self.dti2fs2mni).run()


if __name__ == '__main__':
    print(sys.argv[1])
    f = thalamus_subject(sys.argv[1])
    f.run_recon_all()
    f.convert_mgzs_into_niftis()
    f.dti_preproc()
    f.fs_to_dti_flirt()
    #f.fs_to_dti_fnirt()

    f.mni_to_fs_flirt()
    f.mni_to_fs_fnirt()

    f.mni_to_fs_to_dti_convwarp()
    f.dti_to_fs_to_mni_from_invwarp()
    f.visual_check_nodif_to_mni()

    roiExtraction(f.dir, 'ROI', f.fsDir)





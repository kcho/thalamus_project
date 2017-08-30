# -*- coding: utf-8 -*-

import os
from os.path import join, isfile, isdir
from glob import glob
import re
from nipype.interfaces import fsl
from nipype.interfaces.freesurfer import ReconAll, MRIConvert

class thalamus_subject:
    '''
    f = subject('/subjectdir')

    subjectdir requires below structures

    subjectdir
      ├─ DTI
      │   ├─ 20*.nii.gz : raw dwi
      │   ├─ 20*.bvec : raw bvectors
      │   └─ 20*.bval : raw bvalues
      │ 
      ├─ T1
      │   └─ 20*.nii.gz : raw T1
      │
      └─ DKI (optional)
          ├─ 20*.nii.gz : raw dwi
          ├─ 20*.bvec : raw bvectors
          └─ 20*.bval : raw bvalues
    
    '''

    def __init__(self, subjDir):
        self.fsldir = os.environ['FSLDIR']
        self.atlasdir = join(self.fsldir, 'data', 'atlases')
        self.fslstandard = join(self.fsldir, 'data', 'standard')
        self.mni_thalamus_left = 'lh_thal.nii.gz'
        self.mni_thalamus_right = 'rh_thal.nii.gz'
        self.HO_template = join(self.atlasdir,
                        'HarvardOxford',
                        'HarvardOxford-sub-prob-2mm.nii.gz')

        self.mni = join(self.fslstandard, 'MNI152_T1_2mm.nii.gz')
        self.mni_brain = join(self.fslstandard, 'MNI152_T1_2mm_brain.nii.gz')
        self.mni_brain_mask = join(self.fslstandard, 'MNI152_T1_2mm_brain_mask.nii.gz')

        self.dir = subjDir
        self.t1Dir = join(subjDir, 'T1')
        self.t1raw = join(self.t1Dir, self.get_img_raw(self.t1Dir))
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

        self.mni2fs_mat = join(self.regDir, 'mni2fs.mat')
        self.mni2fs_coeff = join(self.regDir, 'mni2fs_coeff.nii.gz')
        self.mni2fs2dti = join(self.regDir, 'mni2fs2dti_coeff.nii.gz')
        self.dti2fs2mni = join(self.regDir, 'dti2fs2mni_coeff.nii.gz')

        self.dti_nodif2mni_check = join(self.regDir, 
                                        'dti_nodif2mni_check.nii.gz')

        self.bedpostxDir = join(subjDir, 
                                'DTI.bedpostX')


        self.tractDir = join(subjDir, 'tractography')
        self.tractDir_left = join(self.tractDir, 'left')
        self.tractDir_right = join(self.tractDir, 'right')


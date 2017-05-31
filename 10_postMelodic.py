#!/home/kangik/anaconda2/bin/python

import nibabel as nb
import numpy as np
import re
import sys, os
from os.path import join, basename
from scipy.sparse import csr_matrix, coo_matrix
from scipy import stats
import argparse
import textwrap

def grepImgInputs(melodicDir): 
    log = join(melodicDir, 'log.txt')
    with open(log, 'r') as f:
        lines = f.readlines()

    imgInputs = []
    for line in lines:
        try:
            dataLoc = re.search('Reading data file (\S+)', line).group(1)
            imgInputs.append(dataLoc+'.nii.gz')
        except:
            pass
    return imgInputs

def postMelodic(melodicDir):
    '''
    Process melodic outputs
    '''
    # Loading thalamic ROI in MNI space 
    print('\tLoading thalamic ROI') 
    thalamus = nb.load('/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60.nii.gz') 
    thalData = thalamus.get_data() 
    thalNZ = np.nonzero(thalData) 
    thalVoxelNum = len(thalNZ[0]) 

    print('\tLoading track space coordinate')
    probtrackx_dir = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/CHR04_PJH/thalamus_tractography_MNI/left'
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    coord = np.loadtxt(coords_file, dtype='int')

    # Return the index number of the thalamic voxels in ravel format
    thalInd = np.ravel_multi_index((thalNZ[0],thalNZ[1],thalNZ[2]),
                               dims=thalData.shape, 
                               order='F')

    imgInputs = grepImgInputs(melodicDir)
    icMap = join(melodicDir, 'melodic_IC.nii.gz')
    
    # fsl_glm
    for num, imgInput in enumerate(imgInputs[0:2]):
        print('\tRunning fsl_glm processes on {0}'.format(
            re.search('\w{3}\d{2}_\w{3,4}', imgInput).group(0)))
        fsl_reg_out = 'fsl_glm_output_{0}'.format(num)
        
        # fsl_glm 1
        command = 'fsl_glm -i {subjectMap} -d {melodicIC} -o {fsl_reg_out}'.format(
            subjectMap = imgInput,
            melodicIC = icMap,
            fsl_reg_out = fsl_reg_out)
        if not os.path.isfile(fsl_reg_out):
            print('\t\tRunning fsl_glm')
            os.popen(command).read()

        # read output from the fsl_glm
        fsl_glm_mat = np.loadtxt(fsl_reg_out)

        # z-score conversion
        fsl_glm_mat_z = stats.zscore(fsl_glm_mat, axis=1)

        # threshold 
        fsl_glm_threshold = 2

        # Make empty array
        componentNum = fsl_glm_mat.shape[1]

        mask_ravel = np.zeros_like(thalData).ravel()
        mask_ravel_rep = np.tile(mask_ravel[:, np.newaxis], componentNum)
        thalInd_rep = np.tile(thalInd[:, np.newaxis], componentNum)

        # Thalamic coord greater than the threshold
        thalInd_one = thalInd_rep * (fsl_glm_mat_z > fsl_glm_threshold)

        for i in range(componentNum):
            mask_ravel_rep[thalInd_one[:,i], i] = 1

        mask4D = mask_ravel_rep.reshape([thalData.shape[0], thalData.shape[1], thalData.shape[2], 
                                         componentNum], order='F')

        img = nb.Nifti1Image(mask4D, thalamus.affine)
        img.to_filename('{0}_recon.nii.gz'.format(num))

        #print('\tWriting Data')
        #for i in range(melodicMat.shape[0]):
            #mask_ravel[thalInd] = melodicMat[i,:,:].argmax(axis=1)
            ## edit each volume of the 4D array 
            ## with the maps amended as M[i,:] at the coordinates
            ## Edit here
            #mask4D[:,:,:,i] = mask_ravel.reshape(thalData.shape, order='F')

        # Save the results
        #img = nb.Nifti1Image(mask4D, thalamus.affine)
        #img.to_filename(join(melodicDir, 'recon.nii.gz'))

def postMelodic_pre(melodicDir):
    # Loading thalamic ROI in MNI space 
    print('\tLoading thalamic ROI') 
    thalamus = nb.load('/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60.nii.gz') 
    thalData = thalamus.get_data() 
    thalNZ = np.nonzero(thalData) 
    thalVoxelNum = len(thalNZ[0]) 

    print('\tLoading melodic_mix')
    #melodicMix = join(melodicDir, 'melodic_mix')
    melodicMat = np.loadtxt(melodicMix)

    #z-score conversion
    melodicMat = stats.zscore(melodicMat)

    # This statistic was thresholded for each component at 3.1,
    # representing a p-value of less than 0.001 that a thalamic
    # voxel contributes towards an individual component
    melodicMat = melodicMat
    threshold = 3.1
    melodicMat[melodicMat < threshold] = 0

    #melodicMat = melodicMat.reshape(thalVoxelNum, melodicMat.shape[1], melodicMat.shape[0]/thalVoxelNum)
    melodicMat = melodicMat.reshape(melodicMat.shape[0]/thalVoxelNum, thalVoxelNum, melodicMat.shape[1])
    #out = np.concatenate([x[...,np.newaxis] for x in np.split(a, 5)], axis=2)

    # Load coordinate information to save results
    print('\tLoading thalamic coordinate')
    probtrackx_dir = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/CHR04_PJH/thalamus_tractography_MNI/left'
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    coord = np.loadtxt(coords_file, dtype='int')
    thalInd = np.ravel_multi_index((thalNZ[0],thalNZ[1],thalNZ[2]),
                               dims=thalData.shape, 
                               order='F')

    # ravel mask map
    mask_ravel = np.zeros_like(thalData).ravel()
    mask4D = np.tile(np.zeros_like(thalData)[:,:,:,np.newaxis], 
                     melodicMat.shape[0])

    # M.shape[0] : number of voxels used in the tractography
    print('\tWriting Data')
    for i in range(melodicMat.shape[0]):
        mask_ravel[thalInd] = melodicMat[i,:,:].argmax(axis=1)
        # edit each volume of the 4D array 
        # with the maps amended as M[i,:] at the coordinates

        # Edit here
        mask4D[:,:,:,i] = mask_ravel.reshape(thalData.shape, order='F')

    # Save the results
    img = nb.Nifti1Image(mask4D, thalamus.affine)
    img.to_filename(join(melodicDir, 'recon.nii.gz'))


if __name__=='__main__':
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=textwrap.dedent('''\
                {codeName} : Post-process ICA
                1. Load fdt_matrix2.dot
                2. Convert to full matrix
                3. Save it as 4D matrix according to the 
                   x,y,z coordinate in tract_coord.txt
                ========================================
                eg) {codeName} -i /Users/kevin/NOR04_CKI/tractography
                eg) {codeName} -i tractographyDir
                eg) {codeName} -i tractographyDir -t b0_brain.nii.gz
                '''.format(codeName=os.path.basename(__file__))))
    parser.add_argument(
                '-i', '--dir',
                help='ICA output directory, default=pwd',
                default=os.getcwd())
    #parser.add_argument(
                #'-t', '--template',
                #help='Template nii image, with the same dimension as the tractography target space',
                #action='store_true')

    args = parser.parse_args()

    if not args.dir:
        parser.error('No directory given, Please read help')

    #if args.template:
        #voxel_probtrackx(args.dir, args.template)
    #else:
    postMelodic(args.dir)
    #grepImgInputs(args.dir)

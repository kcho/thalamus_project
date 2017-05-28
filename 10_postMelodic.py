#!/home/kangik/anaconda2/bin/python

import nibabel as nb
import numpy as np
import sys, os
from os.path import join, basename
from scipy.sparse import csr_matrix, coo_matrix
import argparse
import textwrap

def postMelodic(melodicDir):
    '''
    Converts ROI to whole brain probabilistic tractography output
    to voxelwise connectivity profile 4D nifti map
    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d3f65c4e.1405
    '''

    thalamus = nb.load('/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60.nii.gz')
    thalData = thalamus.get_data()
    thalNZ = np.nonzero(thalData)
    thalVoxelNum = len(thalNZ)

    probtrackx_dir = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/CHR04_PJH/thalamus_tractography_MNI/left'
    melodicMix = join(melodicDir, 'melodic_mix')
    melodicMat = np.loadtxt(melodicMix)
    melodicMat = np.reshape(thalVoxelNum, melodicMat.shape[0]/thalVoxelNum, melodicMat.shape[1])
    #out = np.concatenate([x[...,np.newaxis] for x in np.split(a, 5)], axis=2)

    # Load coordinate information to save results
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    coord = np.loadtxt(coords_file, dtype='int')
    thalInd = np.ravel_multi_index((thalNZ[0],thalNZ[1],thalNZ[2]),
                               dims=data.shape, 
                               order='F')

    # ravel mask map
    mask_ravel = np.zeros_like(data).ravel()
    mask4D = np.tile(np.zeros_like(data)[:,:,:,np.newaxis], 
                     M.shape[0])

    # M.shape[0] : number of voxels used in the tractography
    print('\tWriting Data')
    for i in range(melodicMat.shape[0]):
        mask_ravel[thalInd] = melodicMat.argmax(axis=1)[i,:,:]
        # edit each volume of the 4D array 
        # with the maps amended as M[i,:] at the coordinates

        # Edit here
        mask4D[:,:,:,i] = mask_ravel.reshape(data.shape, order='F')

    # Save the results
    img = nb.Nifti1Image(mask4D, f.affine)
    img.to_filename(join(probtrackx_dir, 'fdt_matrix2_reconstructed.nii.gz'))


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

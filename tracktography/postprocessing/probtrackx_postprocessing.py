#!/home/kangik/anaconda2/bin/python
# Friday, May 12, 2017

import nibabel as nb
import numpy as np
import sys, os
from os.path import join, basename
from scipy.sparse import csr_matrix, coo_matrix
import argparse
import textwrap


def printTree(directory):
    tabNum=1
    for dirName in directory.split('/')[1:]:
        print('\t'*tabNum + '/{0}'.format(dirName))
        tabNum+=1

def voxel_probtrackx(probtrackx_dir, fdt_paths):
    '''
    Converts ROI to whole brain probabilistic tractography output
    to voxelwise connectivity profile 4D nifti map
    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d3f65c4e.1405
    '''

    fdt_mat2 = join(probtrackx_dir, 'fdt_matrix2.dot')
    #fdt_paths = join(probtrackx_dir, 'fdt_paths.nii.gz')
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')

    # Load Matrix2
    print('\tLoading fdt_mat2.dot')
    try:
        x = np.loadtxt(fdt_mat2, dtype='int')
    except IOError as e:
        printTree(fdt_mat2)
        sys.exit("\tis missing. EXITING.".format(fdt_mat2))

    M = full(spconvert(x))

    # Load 'fdt_path.nii.gz'
    f = nb.load(fdt_paths)
    data = f.get_data()

    # Load coordinate information to save results
    coord = np.loadtxt(coords_file, dtype='int')
    newCoord = coord+1

    ind = np.ravel_multi_index((coord[:,0], coord[:,1], coord[:,2]), 
                               dims=data.shape, 
                               order='F')

    # ravel mask map
    mask_ravel = np.zeros_like(data).ravel()
    mask4D = np.tile(np.zeros_like(data)[:,:,:,np.newaxis], 
                     M.shape[0])

    # M.shape[0] : number of voxels used in the tractography
    print('\tWriting Data')
    for i in range(M.shape[0]):
        mask_ravel[ind] = M[i,:]
        # edit each volume of the 4D array 
        # with the maps amended as M[i,:] at the coordinates

        # Edit here
        mask4D[:,:,:,i] = mask_ravel.reshape(data.shape, order='F')

    # Save the results
    img = nb.Nifti1Image(mask4D, f.affine)
    img.to_filename(join(probtrackx_dir, 'fdt_matrix2_reconstructed.nii.gz'))


def full(DATA):
    '''
    Convert Data to dense matrix
    http://stackoverflow.com/questions/16505670/generating-a-dense-matrix-from-a-sparse-matrix-in-numpy-python
    '''
    return csr_matrix(DATA).todense()

def spconvert(DATA):
    '''
    Convert Data to sparse matrix
    https://gist.github.com/kevinavery/9613505
    '''
    dims = DATA.shape[1] - 1
    shape = [np.max(DATA[:,i]) for i in range(dims)]
    M = np.zeros(shape=shape)
    for row in DATA:
        index = tuple(row[:-1] - 1)
        M.itemset(index, row[-1])
    return M

if __name__=='__main__':
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=textwrap.dedent('''\
                {codeName} : Post-process tractography ran with matrix2 option
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
                help='Tractography directory, default=pwd',
                default=os.getcwd())
    parser.add_argument(
                '-t', '--template',
                help='Template nii image, with the same dimension as the tractography target space',
                action='store_true')

    args = parser.parse_args()

    if not args.dir:
        parser.error('No directory given, Please read help')

    if args.template:
        voxel_probtrackx(args.dir, args.template)
    else:
        voxel_probtrackx(args.dir, join(args.dir, 'fdt_paths.nii.gz'))

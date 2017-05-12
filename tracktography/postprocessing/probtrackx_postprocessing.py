import nibabel as nb
import numpy as np
import sys
from os.path import join, basename
from scipy.sparse import csr_matrix, coo_matrix

def voxel_probtrackx(probtrackx_dir):
    '''
    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d3f65c4e.1405
    '''
    fdt_mat2 = join(probtrackx_dir, 'fdt_matrix2.dot')
    fdt_paths = join(probtrackx_dir, 'fdt_paths.nii.gz')
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')

    # Load Matrix2
    print('Loading fdt_matrix2.dot')
    x = np.loadtxt(fdt_mat2, dtype='int')
    print('\t{0}'.format(x.shape))

    print('Converting fdt_matrix2 to full matrix')
    M = full(spconvert(x))
    print('\t{0}'.format(M.shape))

    # Load coordinate information to save results
    print('Loading fdt_paths.nii.gz')
    f = nb.load(fdt_paths)
    data = f.get_data()
    print('\t{0}'.format(data.shape))

    print('Loading tract_space_coords')
    coord = np.loadtxt(coords_file, dtype='int')
    print('\t{0}'.format(coord.shape))

    print('Ravel x,y,z coordinates into index'.format(coords_file))
    ind = np.ravel_multi_index((coord[:,0], coord[:,1], coord[:,2]), 
                               dims=(data.shape), 
                               order='F')
    print('\t{0}'.format(ind.shape))

    
    mask_ravel = np.zeros_like(data).ravel()
    mask4D = np.tile(np.zeros_like(data)[:,:,:,np.newaxis], 
                     M.shape[0])

    for i in range(M.shape[0]):
        print(i)
        mask_ravel[ind] = M[i,:]
        mask4D[:,:,:,i] = mask_ravel.reshape(data.shape)

    img = nb.Nifti1Image(mask4D, f.affine)
    img.to_filename(join(probtrackx_dir, 'fdt_matrix2_recontructed.nii.gz'))


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
    voxel_probtrackx(sys.argv[1])

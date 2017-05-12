import nibabel as nb
import numpy as np
import sys
from os.path import join
from scipy.sparse import csr_matrix, coo_matrix

def voxel_probtrackx(probtrackx_dir):
    '''
    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d3f65c4e.1405
    '''

    # load Matrix2
    fdt_mat2 = join(probtrackx_dir, 'fdt_matrix2.dot')
    fdt_paths = join(probtrackx_dir, 'fdt_paths.nii.gz')
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    x = np.loadtxt(fdt_mat2, dtype='int')

    M = full(spconvert(x))

    #print(M.shape)

    # Load coordinate information to save results
    f = nb.load(fdt_paths)
    data = f.get_data()

    # it starts from zero
    # need to add one
    # change below
    # Thursday, May 11, 2017
    coord = np.loadtxt(coords_file, dtype='int')
    print(coord.shape)
    #mask4D = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2], M.shape[0]-1))
    #print(mask4D.shape)

    # need to vectorize
    #np.ravel_multi_index([2, 2, 2], dims=(d.shape), order='F')
    ind = np.ravel_multi_index((coord[:,0], coord[:,1], coord[:,2]), 
                               dims=(data.shape), 
                               order='F')
    print(ind.shape)

    mask = np.zeros_like(data)
    mask4D = np.tile(np.zeros_like(data)[:,:,:,np.newaxis], 
                     M.shape[0])

    for i in M.shape[0]:
        #print(i)
        #mask[ind] = M[i,:]
        mask[ind] = 
        mask4D[:,:,:,i] = M[i,:]

    img = nb.Nifti1Image(mask4D, f.affine())
    img.to_filename(join(probtrackx_dir, i+'.nii.gz'))


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

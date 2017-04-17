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
    x = np.loadtxt(fdt_mat2)

    #print(x.shape[0],x.shape[1])
    #sparse_matrix = csr_matrix(x, shape(2, x.shape)#.todense()
    #sparse_matrix = coo_matrix((x, (x.shape[0],x.shape[1])))
    #print(x[:,0].toarray().shape)
    x_cord = x[:,0].astype('int')
    y_cord = x[:,1].astype('int')

    
    #sparse_matrix = csr_matrix((x[:,2], (x_cord, y_cord)))#, shape=(x[:,0].max(), x[:,1].max())).toarray()
    #M = sparse_matrix.todense()
    #M = full(spconvert(x))

    #print(M.shape)

    # Load coordinate information to save results
    f = nb.load(fdt_paths)
    data = f.get_data()

    coord = np.loadtxt(coords_file)
    print(coord.shape)
    #mask4D = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2], M.shape[0]-1))
    #print(mask4D.shape)

    #ind = sub2ind(size(mask), coord(:,1), coord(:,2), coord(:,3));
    ind = np.unravel_index(data.shape, dims=(coord[:,0].astype('int'), coord[:,1].astype('int'), coord[:,2].astype('int')))
    print(ind.shape)

    #for i in M.shape[0]:
        #mask = np.zeros_like(data)
        #mask[ind] = M[i,:]
        #print(i)
        #img = nb.Nifti1Image(mask, f.affine())
        #img.to_filename(join(probtrackx_dir, i+'.nii.gz'))


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

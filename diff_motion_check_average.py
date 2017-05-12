import os
from os.path import join
import numpy as np
import sys

def main(motionDir):
    transF = np.loadtxt(join(motionDir, 'ec_trans.txt'))
    rotF = np.loadtxt(join(motionDir, 'ec_rot.txt'))

    motion_mat = np.concatenate((transF, rotF), axis=1)
    FD_mat = np.sum(np.absolute(motion_mat), axis=1)

    FD_average = FD_mat.mean()
    FD_sum = FD_mat.sum()

    print(FD_average, FD_sum)

    np.savetxt(join(motionDir, 'FD_average.txt'), FD_average.reshape(1,), fmt='%1.3f')
    np.savetxt(join(motionDir, 'FD_sum.txt'), FD_sum.reshape(1,), fmt='%1.3f')

if __name__ == '__main__':
    main(sys.argv[1])

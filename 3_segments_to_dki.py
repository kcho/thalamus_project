import os
from os.path import join, dirname, basename, isfile
import re
import argparse
import sys
from multiprocessing import Pool


def segments_to_dki(args):
    segLoc = 'segmentation'
    kmeanLoc = 'DKI/kmean.nii'
    regMatLoc = 'registration/FREESURFERT1toDKINodif.mat'

    raw_seed_imgs = []
    for subject in args.subjects:
        segDir = join(subject, segLoc)
        for root, dirs, files in os.walk(segDir):
            for f in files:
                if re.search('^seeds_to', f):
                    fLoc = join(root, f)
                    raw_seed_imgs.append(fLoc)

    out_imgs = [join(dirname(x), 'dki_' + basename(x)) for x in raw_seed_imgs]
    refs = [join(x.split(segLoc)[0], kmeanLoc) for x in raw_seed_imgs]
    mats = [join(x.split(segLoc)[0], regMatLoc) for x in raw_seed_imgs]

    pool = Pool()
    out = pool.map(flirt, zip(raw_seed_imgs, out_imgs, refs, mats))


def flirt(parameters):
    imgIn, imgOut, ref, mat = parameters
    if not isfile(imgOut):
        command = 'flirt \
                -in {imgIn} \
                -ref {ref} \
                -applyxfm -init {mat} \
                -out {imgOut}'.format(imgIn = imgIn,
                                      imgOut = imgOut,
                                      ref = ref,
                                      mat = mat)
        out = os.popen(command).read()
        print out

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--subjects', '-s',
                        help = 'list of subjects',
                        nargs='+')
    parser.add_argument('--j', '-j',
                        help = 'number of cores to use',
                        default = 1)
    args = parser.parse_args()
    segments_to_dki(args)

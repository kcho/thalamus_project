import os
from os.path import join, dirname, basename, isfile
import re
import argparse
import sys
from multiprocessing import Pool


def segments_to_dki(args):
    segLoc = 'segmentation'
    kmeanLoc = 'DTI/kmean.nii'
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

#FLIRT version 6.0

#Usage: flirt [options] -in <inputvol> -ref <refvol> -out <outputvol>
       #flirt [options] -in <inputvol> -ref <refvol> -omat <outputmatrix>
       #flirt [options] -in <inputvol> -ref <refvol> -applyxfm -init <matrix> -out <outputvol>

  #Available options are:
        #-in  <inputvol>                    (no default)
        #-ref <refvol>                      (no default)
        #-init <matrix-filname>             (input 4x4 affine matrix)
        #-omat <matrix-filename>            (output in 4x4 ascii format)
        #-out, -o <outputvol>               (default is none)
        #-datatype {char,short,int,float,double}                    (force output data type)
        #-cost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}        (default is corratio)
        #-searchcost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}  (default is corratio)
        #-usesqform                         (initialise using appropriate sform or qform)
        #-displayinit                       (display initial matrix)
        #-anglerep {quaternion,euler}       (default is euler)
        #-interp {trilinear,nearestneighbour,sinc,spline}  (final interpolation: def - trilinear)
        #-sincwidth <full-width in voxels>  (default is 7)
        #-sincwindow {rectangular,hanning,blackman}
        #-bins <number of histogram bins>   (default is 256)
        #-dof  <number of transform dofs>   (default is 12)
        #-noresample                        (do not change input sampling)
        #-forcescaling                      (force rescaling even for low-res images)
        #-minsampling <vox_dim>             (set minimum voxel dimension for sampling (in mm))
        #-applyxfm                          (applies transform (no optimisation) - requires -init)
        #-applyisoxfm <scale>               (as applyxfm but forces isotropic resampling)
        #-paddingsize <number of voxels>    (for applyxfm: interpolates outside image by size)
        #-searchrx <min_angle> <max_angle>  (angles in degrees: default is -90 90)
        #-searchry <min_angle> <max_angle>  (angles in degrees: default is -90 90)
        #-searchrz <min_angle> <max_angle>  (angles in degrees: default is -90 90)
        #-nosearch                          (sets all angular search ranges to 0 0)
        #-coarsesearch <delta_angle>        (angle in degrees: default is 60)
        #-finesearch <delta_angle>          (angle in degrees: default is 18)
        #-schedule <schedule-file>          (replaces default schedule)
        #-refweight <volume>                (use weights for reference volume)
        #-inweight <volume>                 (use weights for input volume)
        #-wmseg <volume>                    (white matter segmentation volume needed by BBR cost function)
        #-wmcoords <text matrix>            (white matter boundary coordinates for BBR cost function)
        #-wmnorms <text matrix>             (white matter boundary normals for BBR cost function)
        #-fieldmap <volume>                 (fieldmap image in rads/s - must be already registered to the reference image)
        #-fieldmapmask <volume>             (mask for fieldmap image)
        #-pedir <index>                     (phase encode direction of EPI - 1/2/3=x/y/z & -1/-2/-3=-x/-y/-z)
        #-echospacing <value>               (value of EPI echo spacing - units of seconds)
        #-bbrtype <value>                   (type of bbr cost function: signed [default], global_abs, local_abs)
        #-bbrslope <value>                  (value of bbr slope)
        #-setbackground <value>             (use specified background value for points outside FOV)
        #-noclamp                           (do not use intensity clamping)
        #-noresampblur                      (do not use blurring on downsampling)
        #-2D                                (use 2D rigid body mode - ignores dof)
        #-verbose <num>                     (0 is least and default)
        #-v                                 (same as -verbose 1)
        #-i                                 (pauses at each stage: default is off)
        #-version                           (prints version number)
        #-help

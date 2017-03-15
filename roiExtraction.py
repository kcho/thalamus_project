#!/usr/bin/python

import os
from os.path import join
import re
import argparse

def roiExtraction(subject, roiName, fsName):
    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = '{0}'.format(subject)

    roiDir = join(subject, roiName)
    fsDir = join(subject, fsName)

    # Change orientation to RAS
    aaLoc = join(fsDir, 'mri/aparc+aseg.mgz')
    oaaLoc = join(fsDir, 'mri/oaparc+aseg.mgz')
    command = 'mri_convert --out_orientation RAS \
            {0} {1}'.format(aaLoc, oaaLoc)
    os.system(command)

    #to extract thalamus
    thalamusNumDict={'lh':10, 'rh':49}
    for side in 'lh','rh':
        command = 'mri_binarize --i {oaaLoc} \
                --match {thalamusNum} \
                --o {roiDir}/{side}_thalamus.nii.gz'.format(
                    oaaLoc = oaaLoc,
                    fsDir=fsDir,
                    thalamusNum=thalamusNumDict[side],
                    roiDir=roiDir,
                    side=side)
        os.system(command)

    # to extract cortex
    cortexNumDict = {'OFC':[1019, 1014, 1012],
                     'MPFC':[1002, 1026, 1028],
                     'LPFC':[1020, 1027, 1032, 1018],
                     'SMC':[1024, 1003, 1022, 1017],
                     'PC':[1008, 1031, 1025, 1023, 1010, 1029],
                     'MTC':[1006, 1016, 1007],
                     'LTC':[1034, 1030, 1001, 1009, 1015, 1033],
                     'OCC':[1021, 1013, 1011, 1005]}
    for roiName, roiNums in cortexNumDict.iteritems():
        for side in 'lh', 'rh':
            if side == 'lh':
                newRoiNum = ' '.join([str(x) for x in roiNums])
            else:
                newRoiNum = ' '.join(['2'+str(x)[1:] for x in roiNums])

            command = 'mri_binarize --i {oaaLoc} \
                            --match {numbers} \
                            --o {roiDir}/{side}_{cortex}.nii.gz'.format(
                                    fsDir=fsDir,
                                    numbers=newRoiNum,
                                    roiDir=roiDir,
                                    side=side,
                                    cortex=roiDir)
            os.system(command)

if __name__=='__main__':
    argparser = argparse.ArgumentParser(prog='roiExtraction',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description='''\
changing values of fillthresh made no change
changing values of proj frac made very small changes of voxel assignment
to one or the other subgroup(layer)
''',epilog="Kevin Cho 2013_07_06")

    argparser.add_argument("--subject","-s",
                           nargs=1,
                           type=str,
                            help='''
                            specify location of the subject folder
                            ''')
    argparser.add_argument("--roiName","-r",
                           nargs=1,
                           type=str,
                            help='''
                            specify ROI output folder name
                            ''')
    argparser.add_argument("--fsName","-f",
                           nargs=1,
                           type=str,
                           help='''
                            specify Freesurfer folder name
                            ''')
    args = argparser.parse_args()

    print args
    roiExtraction(args.subject[0], args.roiName[0], args.fsName[0])

#OFC
#1019    ctx-lh-parsorbitalis                20  100 50  0
#1014    ctx-lh-medialorbitofrontal          200 35  75  0
#1012    ctx-lh-lateralorbitofrontal         35  75  50  0

#MPFC
#1002    ctx-lh-caudalanteriorcingulate      125 100 160 0
#1026    ctx-lh-rostralanteriorcingulate     80  20  140 0
#1028    ctx-lh-superiorfrontal              20  220 160 0

#LPFC
#1020    ctx-lh-parstriangularis             220 60  20  0
#1027    ctx-lh-rostralmiddlefrontal         75  50  125 0
#1032    ctx-lh-frontalpole                  100 0   100 0
#1018    ctx-lh-parsopercularis              220 180 140 0

#SMC
#1024    ctx-lh-precentral                   60  20  220 0
#1003    ctx-lh-caudalmiddlefrontal          100 25  0   0
#1022    ctx-lh-postcentral                  220 20  20  0
#1017    ctx-lh-paracentral                  60  220 60  0

#PC
#1008    ctx-lh-inferiorparietal             220 60  220 0
#1031    ctx-lh-supramarginal                80  160 20  0
#1025    ctx-lh-precuneus                    160 140 180 0
#1023    ctx-lh-posteriorcingulate           220 180 220 0
#1010    ctx-lh-isthmuscingulate             140 20  140 0
#1029    ctx-lh-superiorparietal             20  180 140 0

#MTC
#1006    ctx-lh-entorhinal                   220 20  10  0
#1016    ctx-lh-parahippocampal              20  220 60  0
#1007    ctx-lh-fusiform                     180 220 140 0

#LTC
#1034    ctx-lh-transversetemporal           150 150 200 0
#1030    ctx-lh-superiortemporal             140 220 220 0
#1001    ctx-lh-bankssts                     25  100 40  0
#1009    ctx-lh-inferiortemporal             180 40  120 0
#1015    ctx-lh-middletemporal               160 100 50  0
#1033    ctx-lh-temporalpole                 70  70  70  0

#OCC
#1021    ctx-lh-pericalcarine                120 100 60  0
#1013    ctx-lh-lingual                      225 140 140 0
#1011    ctx-lh-lateraloccipital             20  30  140 0
#1005    ctx-lh-cuneus                       220 20  100 0

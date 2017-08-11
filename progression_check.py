import re
import os
from os.path import basename, join, dirname, abspath, isfile
import numpy as np
import pandas as pd

def freesurfer_check(subjects):
    print '='*80
    print 'Freesurfer check'
    print '='*80

    df = pd.DataFrame()
    for subject in subjects:
        logLoc = join(subject, 'FREESURFER/scripts/recon-all.log')
        with open(logLoc, 'r') as f:
            lastLine = f.readlines()[-1]
        if lastLine.startswith('recon-all -s FREESURFER finished without error'):
            df_subject = pd.DataFrame.from_dict({'subject':[basename(subject)], 'fs':['completed']})
        else:
            df_subject = pd.DataFrame.from_dict({'subject':[basename(subject)], 'fs':['not_completed']})
        df = pd.concat([df, df_subject])
    gb = df.groupby('fs')
    print gb.count()
    print gb.get_group('not_completed').subject.unique()

def tractography_check(subjects):
    print '='*80
    print 'Tractography check'
    print '='*80
    df = pd.DataFrame()
    for subject in subjects:
        for side, sside in zip(['left', 'right'], ['lh', 'rh']):
            segDir = join(subject, 'segmentation_fnirt', side)
            fdt_paths = join(segDir, 'fdt_paths.nii.gz')
            if isfile(fdt_paths):
                df_subject = pd.DataFrame.from_dict({'subject':[basename(subject)], 
                                                     'side':side,
                                                     'segmentation':['completed']})
            else:
                df_subject = pd.DataFrame.from_dict({'subject':[basename(subject)], 
                                                     'side':side,
                                                     'segmentation':['not_completed']})
            
            df = pd.concat([df, df_subject])
    gb = df.groupby(['side', 'segmentation'])
    print gb.count().reset_index().pivot_table(index='segmentation', columns='side')
    print df.groupby('segmentation').get_group('not_completed').subject.unique()

def roi_check(subjects):
    print '='*80
    print 'ROI extraction check'
    print '='*80
    cortices = ['OFC', 'MPFC', 'LPFC', 'SMC', 'PC', 'MTC', 'LTC', 'OCC', 'thalamus']
    df = pd.DataFrame()
    for subject in subjects:
        for side, sside in zip(['left', 'right'], ['lh', 'rh']):
            roiDir = join(subject, 'ROI')
            for cortex in cortices:
                roiLoc = join(roiDir, '{0}_{1}.nii.gz'.format(sside, cortex))
                if isfile(roiLoc):
                    df_subject = pd.DataFrame.from_dict({'subject':[basename(subject)], 
                                                         'side':side,
                                                         'cortex':cortex,
                                                         'roi':['completed']})
                else:
                    df_subject = pd.DataFrame.from_dict({'subject':[basename(subject)], 
                                                         'side':side,
                                                         'cortex':cortex,
                                                         'roi':['not_completed']})
                df = pd.concat([df, df_subject])

    gb = df.groupby(['side', 'cortex', 'roi'])
    print gb.count().reset_index().pivot_table(index='roi', columns='cortex')
    print df.groupby('roi').get_group('not_completed').subject.unique()

if __name__=='__main__':
    subjects = [abspath(x) for x in os.listdir(os.getcwd()) \
                if x.startswith('FEP') or x.startswith('CHR') or x.startswith('NOR')]
    freesurfer_check(subjects)
    tractography_check(subjects)
    roi_check(subjects)

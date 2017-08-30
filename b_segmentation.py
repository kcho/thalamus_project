# -*- coding: utf-8 -*-

from roiExtraction import *
from a_registration import *
import sys

if __name__ == '__main__':
    f = thalamus_subject(sys.argv[1])
    roiExtraction(f.dir, 'ROI', f.fsDir)


##!/bin/sh
#Usage() {
    #echo ""
    #echo "Usage: 2_segmentation <subjDir> <lh>"
    #echo ""
    #exit 1
#}

#[ "$1" = "" ] && Usage

## edit here for different folder structure
#subj=${1}
#fsDir=${subj}/FREESURFER
#regDir=${subj}/registration
#dtiDir=${subj}/DTI
#dkiDir=${subj}/DKI
#roiDir=${subj}/ROI
#segDir=${subj}/segmentation
#bedpostDir=${subj}/DTI.bedpostX
#side_s=$2

#if [ $2 == 'lh' ]
#then
    #side=left
#elif [ $2 == 'rh' ]
#then
    #side=right
#else
    #Usage
#fi

##Segmentation
#if [ ! -d ${segDir} ]
#then
    #mkdir ${segDir}
#else
    #if [ ! -d ${segDir}/${side} ]
    #then
        #mkdir ${segDir}/${side}
    #fi
#fi

#if [ ! -e ${segDir}/${side}/fdt_paths.nii.gz ]
#then
    #rm -rf ${segDir}/${side} 
    #mkdir -p ${segDir}/${side}
    #echo ${roiDir}/${side_s}_LPFC.nii.gz > ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_LTC.nii.gz >> ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_MPFC.nii.gz >> ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_MTC.nii.gz >> ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_OCC.nii.gz >> ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_OFC.nii.gz >> ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_PC.nii.gz >> ${segDir}/${side}/targets.txt
    #echo ${roiDir}/${side_s}_SMC.nii.gz >> ${segDir}/${side}/targets.txt

    #probtrackx \
        #--mode=seedmask \
        #-x ${roiDir}/${side_s}_thalamus.nii.gz \
        #-l \
        #-c 0.2 \
        #-s 2000 \
        #--steplength=0.5 \
        #-P 5000 \
        #--xfm=${regDir}/mniToNodif.mat \
        #--forcedir \
        #--opd \
        #-s ${bedpostDir}/merged \
        #-m ${bedpostDir}/nodif_brain_mask \
        #--dir=${segDir}/${side} \
        #--targetmasks=${segDir}/${side}/targets.txt \
        #--os2t
#else
    #echo ${subj} segmentation on the ${side} done
#fi

## Output thresholding
#if [ ! -e ${segDir}/${side}/biggest.nii.gz ]
#then
    #echo 'find the biggest in ${subj} ${side}'
    #find_the_biggest ${segDir}/${side}/seed* \
        #${segDir}/${side}/biggest.nii.gz
#fi

#for thr in 5 10 20 90 95
#do
    #if [ ! -d ${segDir}/${side}/${thr}thrP ]
    #then
        #mkdir ${segDir}/${side}/${thr}thrP

    #fi

    #if ls ${segDir}/${side}/${thr}thrP/*SMC* 1> /dev/null 2>&1
    #then
        #echo ${subj} segmentation thresholding done
    #else
        
        #for imgs in ${segDir}/${side}/seeds*
        #do
            #filename=`basename "${imgs}"`
            #echo ${subj} ${side} ${filename} ${thr}% thresholding
            #echo ${filename}

            #fslmaths ${imgs} -thrP ${thr} \
                #${segDir}/${side}/${thr}thrP/${thr}_${filename}
        #done
    
        #find_the_biggest ${segDir}/${side}/${thr}thrP/*seeds* \
            #${segDir}/${side}/${thr}thrP/${thr}_biggest.nii.gz
    #fi
#done

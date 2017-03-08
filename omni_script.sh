#!/bin/bash
for subj in $@
do
    # convert the Freesurfer output : brain.mgz --> brain.nii.gz 
    if [ ! -e ${subj}/FREESURFER/mri/brain.nii.gz ]
    then
        mri_convert --out_orientation RAS \
            ${subj}/FREESURFER/mri/brain.mgz \
            ${subj}/FREESURFER/mri/brain.nii.gz
    fi


    # Registration
    if [ ! -d ${subj}/registration ]
    then
        mkdir ${subj}/registration
    fi

    # FREESURFER output T1 --> diffusion tensor space
    if [ ! -e ${subj}/registration/FREESURFERT1toNodif.mat ]
    then
        flirt \
            -in ${subj}/FREESURFER/mri/brain.nii.gz \
            -ref ${subj}/DTI/b0_brain.nii.gz \
            -out ${subj}/registration/FREESURFERT1toNodif \
            -omat ${subj}/registration/FREESURFERT1toNodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # FREESURFER output T1 --> diffusion kurtosis space
    if [ ! -e ${subj}/registration/FREESURFERT1toDKINodif.mat ]
    then
        flirt \
            -in ${subj}/FREESURFER/mri/brain.nii.gz \
            -ref ${subj}/DKI/b0_brain.nii.gz \
            -out ${subj}/registration/FREESURFERT1toDKINodif \
            -omat ${subj}/registration/FREESURFERT1toDKINodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi


    # ROI extraction
    if [ ! -d ${subj}/ROI ] 
    then
        mkdir ${subj}/ROI
    fi

    if [ ! -e ${subj}/ROI/rh_OCC.nii.gz ]
    then
        python roiExtraction.py -S ${subj}
    else
        echo ${subj} 'ROI extraction done'
    fi


    #Segmentation
    if [ ! -d ${subj}/segmentation ]
    then
        mkdir ${subj}/segmentation
    else
        if [ ! -d ${subj}/segmentation/left ]
        then
            mkdir ${subj}/segmentation/left
            mkdir ${subj}/segmentation/right
        fi
    fi

    if [ ! -e ${subj}/segmentation/left/fdt_paths.nii.gz ]
    then
        rm -rf ${subj}/segmentation/left 
        mkdir -p ${subj}/segmentation/left
        echo ${subj}/ROI/lh_LPFC.nii.gz > ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_LTC.nii.gz >> ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_MPFC.nii.gz >> ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_MTC.nii.gz >> ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_OCC.nii.gz >> ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_OFC.nii.gz >> ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_PC.nii.gz >> ${subj}/segmentation/left/targets.txt
        echo ${subj}/ROI/lh_SMC.nii.gz >> ${subj}/segmentation/left/targets.txt
        probtrackx --mode=seedmask -x ${subj}/ROI/lh_thalamus.nii.gz -l -c 0.2 -s 2000 --steplength=0.5 -P 5000 --xfm=${subj}/registration/FREESURFERT1toNodif.mat --forcedir --opd -s ${subj}/DTI.bedpostX/merged -m ${subj}/DTI.bedpostX/nodif_brain_mask --dir=${subj}/segmentation/left --targetmasks=${subj}/segmentation/left/targets.txt --os2t
    else
        echo ${subj} segmentation on the left done
    fi


    if [ ! -e ${subj}/segmentation/right/fdt_paths.nii.gz ]
    then
        rm -rf ${subj}/segmentation/right
        mkdir -p ${subj}/segmentation/right
        echo ${subj}/ROI/rh_LPFC.nii.gz > ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_LTC.nii.gz >> ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_MPFC.nii.gz >> ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_MTC.nii.gz >> ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_OCC.nii.gz >> ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_OFC.nii.gz >> ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_PC.nii.gz >> ${subj}/segmentation/right/targets.txt
        echo ${subj}/ROI/rh_SMC.nii.gz >> ${subj}/segmentation/right/targets.txt
        probtrackx --mode=seedmask -x ${subj}/ROI/rh_thalamus.nii.gz -l -c 0.2 -s 2000 --steplength=0.5 -P 5000 --xfm=${subj}/registration/FREESURFERT1toNodif.mat --forcedir --opd -s ${subj}/DTI.bedpostX/merged -m ${subj}/DTI.bedpostX/nodif_brain_mask --dir=${subj}/segmentation/right --targetmasks=${subj}/segmentation/right/targets.txt --os2t 
    else
        echo ${subj} segmentation on the right done
    fi


    # Output thresholding
    for side in left right
    do
        echo ${side}
        if [ ! -e ${subj}/segmentation/${side}/biggest.nii.gz ]
        then
            echo 'find the biggest in ${subj} ${side}'
            find_the_biggest ${subj}/segmentation/${side}/seed* \
                ${subj}/segmentation/${side}/biggest.nii.gz
        fi

        for thr in 5 10 20 90 95
        do
            if [ ! -d ${subj}/segmentation/${side}/${thr}thrP ]
            then
                mkdir ${subj}/segmentation/${side}/${thr}thrP

            fi

            if [ ! -e ${subj}/segmentation/${side}/${thr}thrP/*SMC* ]
            then
                for imgs in ${subj}/segmentation/${side}/seeds*
                do
                    filename=`basename "${imgs}"`
                    echo ${subj} ${side} ${filename} ${thr}% thresholding
                    echo ${filename}

                    fslmaths ${imgs} -thrP ${thr} \
                        ${subj}/segmentation/${side}/${thr}thrP/${thr}_${filename}
                done
            
                find_the_biggest ${subj}/segmentation/${side}/${thr}thrP/*seeds* \
                    ${subj}/segmentation/${side}/${thr}thrP/${thr}_biggest.nii.gz
            else
                echo ${subj} segmentation thresholding done
            fi
        done
    done
done

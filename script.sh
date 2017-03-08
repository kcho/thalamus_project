#!/bin/bash
for subj in $@
do
    # edit here for different folder structure
    fsDir=${i}/FREESURFER
    regDir=${i}/registration
    dtiDir=${i}/DTI
    dkiDir=${i}/DKI
    roiDir=${i}/ROI
    segDir=${i}/segmentation
    bedpostDir=${i}/DTI.bedpostX

    # convert the Freesurfer output : brain.mgz --> brain.nii.gz 
    if [ ! -e ${fsDir}/mri/brain.nii.gz ]
    then
        mri_convert --out_orientation RAS \
            ${fsDir}/mri/brain.mgz \
            ${fsDir}/mri/brain.nii.gz
    fi


    # Registration
    if [ ! -d ${regDir} ]
    then
        mkdir ${regDir}
    fi

    # FREESURFER output T1 --> diffusion tensor space
    if [ ! -e ${regDir}/FREESURFERT1toNodif.mat ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref ${dtiDir}/b0_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toNodif \
            -omat ${regDir}/FREESURFERT1toNodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi

    # FREESURFER output T1 --> diffusion kurtosis space
    if [ ! -e ${regDir}/FREESURFERT1toDKINodif.mat ]
    then
        flirt \
            -in ${fsDir}/mri/brain.nii.gz \
            -ref ${dkiDir}/b0_brain.nii.gz \
            -out ${regDir}/FREESURFERT1toDKINodif \
            -omat ${regDir}/FREESURFERT1toDKINodif.mat \
            -bins 256 \
            -cost mutualinfo \
            -searchrx -180 180 \
            -searchry -180 180 \
            -searchrz -180 180 \
            -dof 6  -interp trilinear
    fi


    # ROI extraction
    if [ ! -d ${roiDir} ] 
    then
        mkdir ${roiDir}
    fi

    if [ ! -e ${roiDir}/rh_OCC.nii.gz ]
    then
        python roiExtraction.py -S ${subj}
    else
        echo ${subj} 'ROI extraction done'
    fi


    #Segmentation
    if [ ! -d ${segDir} ]
    then
        mkdir ${segDir}
    else
        if [ ! -d ${segDir}/left ]
        then
            mkdir ${segDir}/left
            mkdir ${segDir}/right
        fi
    fi

    if [ ! -e ${segDir}/left/fdt_paths.nii.gz ]
    then
        rm -rf ${segDir}/left 
        mkdir -p ${segDir}/left
        echo ${roiDir}/lh_LPFC.nii.gz > ${segDir}/left/targets.txt
        echo ${roiDir}/lh_LTC.nii.gz >> ${segDir}/left/targets.txt
        echo ${roiDir}/lh_MPFC.nii.gz >> ${segDir}/left/targets.txt
        echo ${roiDir}/lh_MTC.nii.gz >> ${segDir}/left/targets.txt
        echo ${roiDir}/lh_OCC.nii.gz >> ${segDir}/left/targets.txt
        echo ${roiDir}/lh_OFC.nii.gz >> ${segDir}/left/targets.txt
        echo ${roiDir}/lh_PC.nii.gz >> ${segDir}/left/targets.txt
        echo ${roiDir}/lh_SMC.nii.gz >> ${segDir}/left/targets.txt

        probtrackx \
            --mode=seedmask \
            -x ${roiDir}/lh_thalamus.nii.gz \
            -l \
            -c 0.2 \
            -s 2000 \
            --steplength=0.5 \
            -P 5000 \
            --xfm=${regDir}/FREESURFERT1toNodif.mat \
            --forcedir \
            --opd \
            -s ${bedpostDir}/merged \
            -m ${bedpostDir}/nodif_brain_mask \
            --dir=${segDir}/left \
            --targetmasks=${segDir}/left/targets.txt \
            --os2t
    else
        echo ${subj} segmentation on the left done
    fi


    if [ ! -e ${segDir}/right/fdt_paths.nii.gz ]
    then
        rm -rf ${subj}/segmentation/right
        mkdir -p ${subj}/segmentation/right
        echo ${roiDir}/rh_LPFC.nii.gz > ${segDir}/right/targets.txt
        echo ${roiDir}/rh_LTC.nii.gz >> ${segDir}/right/targets.txt
        echo ${roiDir}/rh_MPFC.nii.gz >> ${segDir}/right/targets.txt
        echo ${roiDir}/rh_MTC.nii.gz >> ${segDir}/right/targets.txt
        echo ${roiDir}/rh_OCC.nii.gz >> ${segDir}/right/targets.txt
        echo ${roiDir}/rh_OFC.nii.gz >> ${segDir}/right/targets.txt
        echo ${roiDir}/rh_PC.nii.gz >> ${segDir}/right/targets.txt
        echo ${roiDir}/rh_SMC.nii.gz >> ${segDir}/right/targets.txt

        probtrackx \
            --mode=seedmask \
            -x ${roiDir}/rh_thalamus.nii.gz \
            -l \
            -c 0.2 \
            -s 2000 \
            --steplength=0.5 \
            -P 5000 \
            --xfm=${regDir}/FREESURFERT1toNodif.mat \
            --forcedir \
            --opd \
            -s ${bedpostDir}/merged \
            -m ${bedpostDir}/nodif_brain_mask \
            --dir=${segDir}/right \
            --targetmasks=${segDir}/right/targets.txt \
            --os2t 
    else
        echo ${subj} segmentation on the right done
    fi


    # Output thresholding
    for side in left right
    do
        echo ${side}
        if [ ! -e ${segDir}/${side}/biggest.nii.gz ]
        then
            echo 'find the biggest in ${subj} ${side}'
            find_the_biggest ${segDir}/${side}/seed* \
                ${segDir}/${side}/biggest.nii.gz
        fi

        for thr in 5 10 20 90 95
        do
            if [ ! -d ${segDir}/${side}/${thr}thrP ]
            then
                mkdir ${segDir}/${side}/${thr}thrP

            fi

            if [ ! -e ${segDir}/${side}/${thr}thrP/*SMC* ]
            then
                for imgs in ${segDir}/${side}/seeds*
                do
                    filename=`basename "${imgs}"`
                    echo ${subj} ${side} ${filename} ${thr}% thresholding
                    echo ${filename}

                    fslmaths ${imgs} -thrP ${thr} \
                        ${segDir}/${side}/${thr}thrP/${thr}_${filename}
                done
            
                find_the_biggest ${segDir}/${side}/${thr}thrP/*seeds* \
                    ${segDir}/${side}/${thr}thrP/${thr}_biggest.nii.gz
            else
                echo ${subj} segmentation thresholding done
            fi
        done
    done
done

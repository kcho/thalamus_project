#!/ccnc/anaconda2/bin/python

import nibabel as nb
import multiprocessing
import numpy as np
import re
import sys, os
from os.path import join, basename
from scipy.sparse import csr_matrix, coo_matrix
from scipy import stats
import argparse
import textwrap
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from nilearn import plotting
from nilearn.image import iter_img
from copy import copy


#def graph(melodic_ICf, meanMapD):
def graph(icMap_nb, sumMap, affine):
    mni = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
    mni_data = nb.load(mni).get_data()

    my_cmap = copy(matplotlib.cm.get_cmap('jet')) # make a copy so we don't mess up system copy
    my_cmap.set_under('r', alpha=0) # make locations over vmax translucent red
    #my_cmap.set_bad('white', 0)

    #fig, axes = plt.subplots(ncols = 2, figsize=(10,10))
    #print melodic_ICf.get_data().shape
    #axes[0].imshow(meanMapD[20,:,:,0])
    #axes[1].imshow(melodic_ICf.get_data()[20,:,:,0])
    ##axes[2].imshow([20,20,20,0])
    #plt.savefig('prac.png')
    sumMap_img = nb.Nifti1Image(sumMap, affine)
            #img = nb.Nifti1Image(mask4D, thalamus.affine)

    for (threeD, meanImg), ICout in zip(enumerate(iter_img(sumMap_img)), iter_img(icMap_nb)):
        print threeD+1
        fig = plt.figure(figsize=(20,7))
        gs = gridspec.GridSpec(2, 5)
        ax1 = plt.subplot(gs[0, :])
        ax2 = plt.subplot(gs[1, 0])
        ax3 = plt.subplot(gs[1, 1])
        ax4 = plt.subplot(gs[1, 2])
        ax5 = plt.subplot(gs[1, 3])
        ax6 = plt.subplot(gs[1, 4])
        axes= [ax2, ax3, ax4, ax5, ax6]

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        #     plotting.plot_glass_brain(meanImg, threshold=2, axes=axes[0])
        #plotting.plot_glass_brain(ICout, threshold=1, title='Component %d' %(threeD+1), axes=ax1)
        plotting.plot_glass_brain(ICout, title='Component %d' %(threeD+1), axes=ax1)

        startSlice = 35
        for axNum in range(5):
            axes[axNum].imshow(np.flipud(mni_data[15:75,40:75,startSlice].T), cmap=matplotlib.cm.gray)

            s = axes[axNum].imshow(np.flipud(sumMap[15:75,40:75,startSlice, threeD].T), 
                                   #cmap=my_cmap)
                                   cmap=plt.get_cmap('hot'), clim=(-10, 10))
                                   #vmax=2, 
                                   #vmin=1)
            #    axes[axNum].imshow(averageMap[:,:, startSlice, threeD], cmap=my_cmap, alpha=0.4)
            startSlice+=3
            #     cbar=plt.colorbar(s, ax=axes[-1])
            #     fig.suptitle('component '+ str(threeD+1))
        fig.colorbar(s, cax=cbar_ax)
        plt.savefig('/home/kangik/IC_component_4resample_10component_{0:02d}.png'.format(threeD+1))
        #     cbar=plt.colorbar(s, ax=axes[-1])
        #plotting.show()

def fslglm(inputList):
    num, imgInput, icMap = inputList
    fsl_reg_out = 'fsl_glm_output_{0}'.format(num)

    # fsl_glm
    command = 'fsl_glm -i {subjectMap} -d {melodicIC} -o {fsl_reg_out}'.format(
        subjectMap = imgInput,
        melodicIC = icMap,
        fsl_reg_out = fsl_reg_out)

    if not os.path.isfile(fsl_reg_out):
        os.popen(command).read()

def grepImgInputs(melodicDir): 
    log = join(melodicDir, 'log.txt')
    with open(log, 'r') as f:
        lines = f.readlines()

    imgInputs = []
    for line in lines:
        try:
            dataLoc = re.search('Reading data file (\S+)', line).group(1)
            imgInputs.append(dataLoc+'.nii.gz')
        except:
            pass
    return imgInputs

def postMelodic(melodicDir):
    '''
    Process melodic outputs
    '''
    # Loading thalamic ROI in MNI space 
    print('\tLoading thalamic ROI') 
    thalamus = nb.load('/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60.nii.gz') 
    thalData = thalamus.get_data() 
    thalNZ = np.nonzero(thalData) 
    thalVoxelNum = len(thalNZ[0]) 

    print('\tLoading track space coordinate')
    probtrackx_dir = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/CHR04_PJH/thalamus_tractography_MNI/left'
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    coord = np.loadtxt(coords_file, dtype='int')

    # Return the index number of the thalamic voxels in ravel format
    thalInd = np.ravel_multi_index((thalNZ[0],thalNZ[1],thalNZ[2]),
                               dims=thalData.shape, 
                               order='F')

    imgInputs = grepImgInputs(melodicDir)
    icMap = join(melodicDir, 'melodic_IC.nii.gz')
    icMap_nb = nb.load(icMap)
    

    # Mixture-modelling outputs
    # -------------------------
    # Each spatial independent component found by melodic was thresholded 
    # using mixture modelling (Beckmann et al., 2003), 
    #    - with a Gaussian distribution used to model background noise and 
    #    - a single gamma distribution used to model signal of interest, 
    #      in this case modelling the positive tail of the distribution. 
    #      The gamma distribution was thresholded at 
    #     a posterior probability of p < 0.5 
    mmthreshDir = join(melodicDir, 'stats')
    mmthreshImgs = [join(mmthreshDir, x) for x in os.listdir(mmthreshDir)]

    # The contribution of each tractogram to each independent component was converted into a normalised Z statistic, 
    #     (output from fsl_glm for every subjects)
    #     (and simultaneously the voxel from which the tractogram is calculated) 
    # calculated as the distance of a particular tractogram's representation 
    # from the mean representation of all tractograms in units of standard deviation. 
    #    z-score conversion
    #    z-score conversion (axis=1) --> does not survive the threshold
    # This allows the localisation of the seed voxels that contribute towards an independent component, 
    # in this case a presumed white matter pathway, 
    # as well as the quantification of the size of the effect of their contribution. 
    # This method assumes that thalamic clusters are distributed similarly across the healthy controls tested (which is tested in study 2). 
    # This statistic was thresholded 
    #     for each component at 3.1, 
    # representing a p-value of less than 0.001 that a thalamic voxel contributes towards an individual component.
    # 3.1

    # run fsl_glm
    inputList = []
    for num, imgInput in enumerate(imgInputs):
        inputList.append((num, imgInput, icMap))
    pool = multiprocessing.Pool(5)
    print('\t\tRunning fsl_glm in parallel')
    for i,_ in enumerate(pool.imap_unordered(fslglm, inputList), 1):
        sys.stderr.write('\rProgress {0:%}'.format(i/len(inputList)))
    print()

    # read outputs from the fsl_glm
    componentNum = icMap_nb.shape[3]

    fsl_glm_mat_subj = np.zeros(thalVoxelNum * componentNum * len(imgInputs))
    fsl_glm_mat_subj = fsl_glm_mat_subj.reshape(thalVoxelNum, componentNum, len(imgInputs))

    print fsl_glm_mat_subj.shape
    print len(imgInputs)

    if not os.path.isfile('sumMap.npy'):
        # for each subject fsl_glm outputs
        for num, imgInput in enumerate(imgInputs):
            print num
            # read output from the fsl_glm
            fsl_reg_out = 'fsl_glm_output_{0}'.format(num)
            fsl_glm_mat = np.loadtxt(fsl_reg_out)

            # z-score conversion
            # axis=1 option states to estimate z-scores 
            # from the samples in the same component
            # fsl_glm_mat_z = stats.zscore(fsl_glm_mat, axis=1) --> does not survive threshold 3.1
            fsl_glm_mat_z = stats.zscore(fsl_glm_mat)
            fsl_glm_threshold = 3.1
            #fsl_glm_mat_z[fsl_glm_mat_z < fsl_glm_threshold] = 0

            # Make empty array
            # mask_ravel_rep : whole brain ravel x component number
            # thalInd_rep : thalamic indices x component number 
            #mask_ravel_rep = np.tile(np.zeros_like(thalData).ravel()[:, np.newaxis], componentNum)
            #thalInd_rep = np.tile(thalInd[:, np.newaxis], componentNum)

            # Threshold the z-score map
            # Thalamic coord greater than the threshold

            # save to fsl_glm_mat_subj matrix
            fsl_glm_mat_subj[:,:, num] = fsl_glm_mat_z

            #mask_ravel_rep[thalInd, :] = fsl_glm_mat_z

            #print('\t\tWriting IC image for the subject')
            #mask4D = mask_ravel_rep.reshape([thalData.shape[0], thalData.shape[1], thalData.shape[2], 
                                             #componentNum], order='F')
            #Add to mean matrix
            #sumMap += mask4D
            #img = nb.Nifti1Image(mask4D, thalamus.affine)
            #img.to_filename('{0}_IC.nii.gz'.format(num))
        #np.save('sumMap', sumMap)
        np.save('fsl_glm_mat_subj', fsl_glm_mat_subj)


    fsl_glm_mat_subj = np.load('fsl_glm_mat_subj.npy')
    # Make empty array
    # mask_ravel_rep : whole brain ravel x component number
    # mask_ravel_rep_sub : whole brain ravel x component number x subject number
    # thalInd_rep : thalamic indices x component number 
    # thalInd_rep_sub : thalamic indices x component number x subject number
    mask_ravel_rep = np.tile(np.zeros_like(thalData).ravel()[:, np.newaxis], componentNum)
    mask_ravel_rep_sub = np.tile(mask_ravel_rep[:,:, np.newaxis], len(imgInputs))

    #mask_ravel_rep_sub[thalInd_rep_sub, :] = fsl_glm_mat_subj
    mask_ravel_rep_sub[thalInd, :,:] = fsl_glm_mat_subj
    print mask_ravel_rep_sub.shape
    mask5D = mask_ravel_rep_sub.reshape([thalData.shape[0], thalData.shape[1], thalData.shape[2], 
                                      componentNum, len(imgInputs)], order='F')
    
    fig, axes = plt.subplots(ncols=4, nrows=2, figsize=(20,5))
    [a,b,c,d,e,f,g,h,i,j] = axes[0][0].plot(fsl_glm_mat_subj[:,:,0], label='first subject')
    axes[0][1].plot(fsl_glm_mat_subj[:,:,1], label='second subject')
    axes[0][2].plot(fsl_glm_mat_subj[:,:,2], label='third subject')
    axes[0][3].plot(fsl_glm_mat_subj[:,:,3], label='fourth subject')

    plt.legend([a,b,c,d,e,f,g,h,i,j], ['1','2','3','4','5','6','7','8','9','10'], loc=1)

    axes[1][0].imshow(mask5D[15:75, 40:75, 35, 4, 0], label='first subject, first component')
    axes[1][1].imshow(mask5D[15:75, 40:75, 40, 4, 0], label='first subject, second component')
    axes[1][2].imshow(mask5D[15:75, 40:75, 45, 4, 0], label='second subject, first component')
    axes[1][3].imshow(mask5D[15:75, 40:75, 50, 4, 0], label='second subject, second component')


    #plt.savefig('prac_thal.png')

    #sumMap = np.load('sumMap.npy')
    # z-score conversion
    # axis=1 option states to estimate z-scores 
    # from the samples in the same component
    #sumMap = stats.zscore(sumMap, axis=1)


    # Make graph
    #graph(icMap_nb, sumMap, thalamus.affine)

    #img = nb.Nifti1Image(sumMap, thalamus.affine)
    #img.to_filename('sum_IC.nii.gz'.format(num))

def postMelodic_pre(melodicDir):
    # Loading thalamic ROI in MNI space 
    print('\tLoading thalamic ROI') 
    thalamus = nb.load('/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60.nii.gz') 
    thalData = thalamus.get_data() 
    thalNZ = np.nonzero(thalData) 
    thalVoxelNum = len(thalNZ[0]) 

    print('\tLoading melodic_mix')
    #melodicMix = join(melodicDir, 'melodic_mix')
    melodicMat = np.loadtxt(melodicMix)

    #z-score conversion
    melodicMat = stats.zscore(melodicMat)

    # This statistic was thresholded for each component at 3.1,
    # representing a p-value of less than 0.001 that a thalamic
    # voxel contributes towards an individual component
    melodicMat = melodicMat
    threshold = 3.1
    melodicMat[melodicMat < threshold] = 0

    #melodicMat = melodicMat.reshape(thalVoxelNum, melodicMat.shape[1], melodicMat.shape[0]/thalVoxelNum)
    melodicMat = melodicMat.reshape(melodicMat.shape[0]/thalVoxelNum, thalVoxelNum, melodicMat.shape[1])
    #out = np.concatenate([x[...,np.newaxis] for x in np.split(a, 5)], axis=2)

    # Load coordinate information to save results
    print('\tLoading thalamic coordinate')
    probtrackx_dir = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/CHR04_PJH/thalamus_tractography_MNI/left'
    coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    coord = np.loadtxt(coords_file, dtype='int')
    thalInd = np.ravel_multi_index((thalNZ[0],thalNZ[1],thalNZ[2]),
                               dims=thalData.shape, 
                               order='F')

    # ravel mask map
    mask_ravel = np.zeros_like(thalData).ravel()
    mask4D = np.tile(np.zeros_like(thalData)[:,:,:,np.newaxis], 
                     melodicMat.shape[0])

    # M.shape[0] : number of voxels used in the tractography
    print('\tWriting Data')
    for i in range(melodicMat.shape[0]):
        mask_ravel[thalInd] = melodicMat[i,:,:].argmax(axis=1)
        # edit each volume of the 4D array 
        # with the maps amended as M[i,:] at the coordinates

        # Edit here
        mask4D[:,:,:,i] = mask_ravel.reshape(thalData.shape, order='F')

    # Save the results
    img = nb.Nifti1Image(mask4D, thalamus.affine)
    img.to_filename(join(melodicDir, 'recon.nii.gz'))



def makeGraph(melodicDir):
    kkkkkkkkkkkkkkkkkkkk
if __name__=='__main__':
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=textwrap.dedent('''\
                {codeName} : Post-process ICA
                1. Load fdt_matrix2.dot
                2. Convert to full matrix
                3. Save it as 4D matrix according to the 
                   x,y,z coordinate in tract_coord.txt
                ========================================
                eg) {codeName} -i /Users/kevin/NOR04_CKI/tractography
                eg) {codeName} -i tractographyDir
                eg) {codeName} -i tractographyDir -t b0_brain.nii.gz
                '''.format(codeName=os.path.basename(__file__))))
    parser.add_argument(
                '-i', '--dir',
                help='ICA output directory, default=pwd',
                default=os.getcwd())
    #parser.add_argument(
                #'-t', '--template',
                #help='Template nii image, with the same dimension as the tractography target space',
                #action='store_true')

    args = parser.parse_args()

    if not args.dir:
        parser.error('No directory given, Please read help')

    #if args.template:
        #voxel_probtrackx(args.dir, args.template)
    #else:
    postMelodic(args.dir)
    #grepImgInputs(args.dir)

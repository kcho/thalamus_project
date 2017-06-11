#!/ccnc/anaconda2/bin/python
from __future__ import division
import nibabel as nb
import multiprocessing
import numpy as np
import re
import sys, os
from os.path import join, basename, isfile, isdir
from scipy.sparse import csr_matrix, coo_matrix
from scipy import stats
import argparse
import textwrap
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.cm as cmx
import matplotlib.colors as colors


from nilearn import plotting
from nilearn.image import iter_img
from copy import copy


#def graph(melodic_ICf, meanMapD):
def graph(icMap_nb, sumMap, affine, mniLoc):
    mni_data = nb.load(mniLoc).get_data()

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
    num, imgInput, icMap, melodicDir = inputList
    outDir = join(melodicDir, 'glm_out_dir')

    try:
        os.mkdir(outDir)
    except:
        pass

    fsl_reg_out = join(outDir, 'fsl_glm_output_{0}'.format(num))

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
    return set(imgInputs)

def postMelodic(melodicDir):
    '''
    Process melodic outputs
    '''
    imgInputs = grepImgInputs(melodicDir)
    subjectNum = len(imgInputs)

    icMap_loc = join(melodicDir, 'melodic_IC.nii.gz')
    icMap_nb = nb.load(icMap_loc)
    componentNum = icMap_nb.shape[3]

    voxSize = icMap_nb.header.get_zooms()[0]
    #if voxSize == 4:
        #mniLoc = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/MNI152_T1_2mm_ds_mask.nii.gz'
        #thalamusLoc = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60_ds.nii.gz'
    #elif voxSize == 3:
        #mniLoc = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/MNI152_T1_2mm_ds_3_mask.nii.gz'
        #thalamusLoc = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60_ds_3.nii.gz'
    #else:
    mniLoc = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm.nii.gz'
    thalamusLoc = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/scripts/lh_thalamus_HOSC_60.nii.gz'

    # Loading thalamic ROI in MNI space 
    print('\tLoading thalamic ROI') 

    thalamus = nb.load(thalamusLoc) 
    thalData = thalamus.get_data() 
    thalNZ = np.nonzero(thalData) 
    thalVoxelNum = len(thalNZ[0]) 

    # Return the index number of the thalamic voxels in ravel format
    thalInd = np.ravel_multi_index((thalNZ[0],thalNZ[1],thalNZ[2]),
                               dims=thalData.shape, 
                               order='F')


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
    mmthreshImgs = [join(mmthreshDir, x) for x in os.listdir(mmthreshDir) if x.startswith('thresh_zstat')]

    # mmthresh to $FSLDIR/data/standard/MNI152_T1_0.5mm
    for mmthreshImgLoc in mmthreshImgs:
        ref='$FSLDIR/data/standard/MNI152_T1_0.5mm'
        mmthreshImg = basename(mmthreshImgLoc)
        mmthreshNum = re.search('\d', mmthreshImg).group(0)
        hresOut = join(mmthreshDir, 'hres_{0}.nii.gz'.format(mmthreshNum))
        overlayOut = join(mmthreshDir, 'overlay_{0}.nii.gz'.format(mmthreshNum))
        slicerOut = join(mmthreshDir, 'slicer_{0}.ppm'.format(mmthreshNum))
        imgOut = join(mmthreshDir, '{0}_sliced.png'.format(mmthreshNum))
        if not isfile(imgOut):
            command = 'flirt -in {mmthresh} \
                    -ref {ref} \
                    -applyxfm -out {hresOut}'.format(
                        mmthresh = mmthreshLoc, ref = ref, hresOut = hresOut)
            os.popen(command).read()
            command = 'overlay 1 0 {ref} \
                    -A {hresOut} \
                    2.5 10 {overlayOut}'.format(
                        ref = ref, hresOut = hresOut, overlayOut = overlayOut)
            os.popen(command).read()
            command = 'slicer {overlayOut} \
                    -s 24 1200 \
                    {slicerOut}'.format(
                        overlayOut = overlayOut, slicerOut = slicerOut)
            os.popen(command).read()
            command = 'convert {slicerOut} {imgOut}'.format(
                slicerOut = slicerOut, imgOut = imgOut)
            os.popen(command).read()

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

    # Running fsl_glm
    # Make inputList for parallel processing
    inputList = []
    for num, imgInput in enumerate(imgInputs):
        inputList.append((num, imgInput, icMap_loc, melodicDir))
    pool = multiprocessing.Pool(5)
    print('\tRunning fsl_glm in parallel')

    # Run fsl_glm in parallel
    for i,_ in enumerate(pool.imap_unordered(fslglm, inputList), 1):
        sys.stderr.write('\r\tProgress {0:%}'.format(i/len(inputList)))
    print()

    # Reading outputs from the fsl_glm
    # Make a containner
    fsl_glm_mat_sub = np.zeros(thalVoxelNum * componentNum * subjectNum)
    fsl_glm_mat_sub = fsl_glm_mat_sub.reshape(thalVoxelNum, componentNum, subjectNum)

    # for each subject, read fsl_glm outputs
    fsl_glm_mat_loc = join(melodicDir, 'fsl_glm_mat_sub.npy')
    if not isfile(fsl_glm_mat_loc):
        for num, imgInput in enumerate(imgInputs):
            fsl_glm_mat_z = read_glm_outputs(melodicDir, num)
            print(fsl_glm_mat_z.shape)
            print(fsl_glm_mat_sub.shape)
            fsl_glm_mat_sub[:,:, num] = fsl_glm_mat_z
        np.save(join(melodicDir, 'fsl_glm_mat_sub'), fsl_glm_mat_sub)
    else:
        fsl_glm_mat_sub = np.load(join(melodicDir, 'fsl_glm_mat_sub.npy'))

    # Make empty array
    # mask_ravel_rep : whole brain ravel x component number
    # mask_ravel_rep_sub : whole brain ravel x component number x subject number
    mask_ravel_rep = np.tile(np.zeros_like(thalData).ravel()[:, np.newaxis], componentNum)
    mask_ravel_rep_sub = np.tile(mask_ravel_rep[:,:, np.newaxis], subjectNum)
    mask_ravel_rep_sub[thalInd, :,:] = fsl_glm_mat_sub

    mask5D = mask_ravel_rep_sub.reshape([thalData.shape[0], 
                                         thalData.shape[1], 
                                         thalData.shape[2], 
                                         componentNum, 
                                         subjectNum], order='F')
    
    draw_max_segmentation(fsl_glm_mat_sub, thalData, thalInd)

    #plt.savefig('prac_thal.png')

    #sumMap = np.load('sumMap.npy')
    # z-score conversion
    # axis=1 option states to estimate z-scores 
    # from the samples in the same component
    #sumMap = stats.zscore(sumMap, axis=1)

    # Make graph
    #graph(icMap_nb, sumMap, thalamus.affine, mniLoc)

    #img = nb.Nifti1Image(sumMap, thalamus.affine)
    #img.to_filename('sum_IC.nii.gz'.format(num))

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def draw_max_segmentation(mat, thalData, thalInd):
    icNum = mat.shape[1]
    compMax = mat.mean(axis=2).argmax(axis=1)
    threeDmap = np.zeros_like(thalData).ravel()
    threeDmap[thalInd] = compMax[:]
    threeDmap = threeDmap.reshape([thalData.shape[0], 
                                   thalData.shape[1], 
                                   thalData.shape[2]], order='F')

    fig = plt.figure(figsize=(20,15))
    gs = gridspec.GridSpec(4,5)
    ax = plt.subplot(gs[0,:])
    ax1 = plt.subplot(gs[1,0])
    ax2 = plt.subplot(gs[1,1])
    ax3 = plt.subplot(gs[1,2])
    ax4 = plt.subplot(gs[1,3])
    ax5 = plt.subplot(gs[1,4])
    ax6 = plt.subplot(gs[2,0])
    ax7 = plt.subplot(gs[2,1])
    ax8 = plt.subplot(gs[2,2])
    ax9 = plt.subplot(gs[2,3])
    ax10 = plt.subplot(gs[2,4])

    baxes = [ax1, ax2, ax3, ax4, ax5, 
             ax6, ax7, ax8, ax9, ax10]

    cmaps=get_cmap(icNum+1)
    cmap_dict = {}
    cmap_dict[0]= 'white'
    for compNum in range(1, icNum+1):
        col = cmaps(compNum)
        cmap_dict[compNum] = col

    comp_count = {}
    for voxelNum, value in enumerate(compMax):
        col = cmap_dict[value+1]
        rect = plt.Rectangle((voxelNum, -0.5), 1, 1, facecolor=col)
        ax.add_artist(rect)

        try:
            comp_count[value+1]+=1
        except:
            comp_count[value+1]=1

    # Brain
    vmax=3
    cmaps = colors.LinearSegmentedColormap.from_list('hsv', 
                                                     [(compNum/(icNum), col) for compNum,col in cmap_dict.iteritems()]
                                                    )
                                                     #N=icNum)
    startSlice=35
    for axNum in range(10):
        s = baxes[axNum].imshow(np.flipud(threeDmap[15:75, 40:75, startSlice].T),
                                cmap = cmaps)
        baxes[axNum].set_title(startSlice)
        baxes[axNum].axes.yaxis.set_ticklabels([])
        baxes[axNum].axes.xaxis.set_ticklabels([])
        startSlice+=1

    #legend
    patches = [mpatches.Patch(color=cmap_dict[x], 
    label='Component {0} : {1} voxels'.format(x,comp_count[x])) for x in range(1, icNum+1)]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=3, borderaxespad=0.)

    ax.set_title('Assigning each thalamic voxel with the strongest ICA component', 
    fontsize=20)
    ax.set_xlim(0,len(compMax))
    ax.set_xlabel('Thalamic voxels in order')
    ax.set_ylim(0,0.5)
    ax.axes.yaxis.set_ticklabels([])
    plt.savefig('max_segmentation.png')
    #fig.show()

def read_glm_outputs(melodicDir, num):
    # read output from the fsl_glm

    fsl_reg_out = join(melodicDir, 'glm_out_dir', 'fsl_glm_output_{0}'.format(num))
    fsl_glm_mat = np.loadtxt(fsl_reg_out)

    # z-score conversion
    # axis=1 option states to estimate z-scores 
    # from the samples in the same component
    # fsl_glm_mat_z = stats.zscore(fsl_glm_mat, axis=1) --> does not survive threshold 3.1
    fsl_glm_mat_z = stats.zscore(fsl_glm_mat)
    return fsl_glm_mat_z

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
    thalInd = np.ravel_multi_index((thalNZ[0],
                                    thalNZ[1],
                                    thalNZ[2]), dims=thalData.shape, order='F')
    #probtrackx_dir = '/Volume/CCNC_W1_2T/2017_CHR_thalamus_microstructure_kcho/allData/CHR04_PJH/thalamus_tractography_MNI/left'
    #coords_file = join(probtrackx_dir, 'tract_space_coords_for_fdt_matrix2')
    #coord = np.loadtxt(coords_file, dtype='int')


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

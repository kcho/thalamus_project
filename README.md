# thalamus_project

2017_03_08

## 1_preprocessing
- Freesurfer space T1 to nifti
- Registration
  - fs T1 --> DTI
  - fs T1 --> DKI
  - DTI --> DKI [to be done]
- ROI extraction using `roiExtraction.py`

```sh
bash 1_preprocessing.sh /home/kangik/subj01
```


## 2_segmentation
- Connectivity based thalamus segmentation
- Finding the subregions (find_the_biggest function from fsl)
- Thresholding connectivity maps (fslmaths -thrP)

```sh
bash 2_segmentation.sh /home/kangik/subj01
```

## 3_segments_to_dki.py*
- registers Freesurfer spaced connectivity maps (seeds_in_*) to DKI space

```sh
python 3_segments_to_dki.py -s ../*
```

 

## For parallel processing

1. Edit name variables in each shell files
2. Use GNU parallel as below

```sh
for i in /home/kangik/subj*
do
  echo bash 1_preprocessing.sh ${i}
done|parallel
```

#######################################################################################################################################
#######################################################################################################################################
## This code is used  OASIS segmentations
## It requires a brain mask (DuraMask.nii.gz), a FLAIR ('FLAIR.nii.gz'), a proton density ('PD.nii.gz'), a T2-w volume ('T2.nii.gz')
## and a T1-w volume (VolumetricT1.nii.gz')
## The code calls fsl and R
##
## Author: Elizabeth Sweeney
## Last updated: June 16, 2013
#######################################################################################################################################
#######################################################################################################################################

#!/bin/bash

#######################################################################################################################################
## Dilate first segmentation and apply csf mask
#######################################################################################################################################

fslmaths '_segment_run1.nii' -kernel boxv 5  -dilM '_segment_dilate'
fslmaths 'csfmask.nii.gz' -sub '_segment_dilate' -mas 'csfmask.nii.gz' 'lesion_segment_dialte.nii.gz'

#######################################################################################################################################
## Create second smoothed volumes (kernel window size of 10 and 20)
#######################################################################################################################################

fslmaths  'FLAIRnorm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 10  'FLAIRnorm_blur2_10.nii.gz'
fslmaths  'FLAIRnorm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 20 'FLAIRnorm_blur2_20.nii.gz'
fslmaths  'PDnorm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 10  'PDnorm_blur2_10.nii.gz'
fslmaths  'PDnorm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 20 'PDnorm_blur2_20.nii.gz'
fslmaths  'T2norm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 10  'T2norm_blur2_10.nii.gz'
fslmaths  'T2norm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 20 'T2norm_blur2_20.nii.gz'
fslmaths  'T1norm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 10  'T1norm_blur2_10.nii.gz'
fslmaths  'T1norm.nii.gz' -mas 'lesion_segment_dialte.nii.gz' -s 20 'T1norm_blur2_20.nii.gz'
fslmaths  'lesion_segment_dialte.nii.gz' -s 10  'SegmentDilate_blur2_10.nii.gz' 
fslmaths  'lesion_segment_dialte.nii.gz' -s 20  'SegmentDilate_blur2_20.nii.gz'  
fslmaths  'FLAIRnorm_blur2_10.nii.gz' -div 'SegmentDilate_blur2_10.nii.gz'   -mas 'csfmask.nii.gz' 'FLAIR_blur2_10_div.nii.gz'
fslmaths  'FLAIRnorm_blur2_20.nii.gz' -div 'SegmentDilate_blur2_20.nii.gz' -mas 'csfmask.nii.gz' 'FLAIR_blur2_20_div.nii.gz'
fslmaths  'PDnorm_blur2_10.nii.gz' -div 'SegmentDilate_blur2_10.nii.gz'   -mas 'csfmask.nii.gz' 'PD_blur2_10_div.nii.gz'
fslmaths  'PDnorm_blur2_20.nii.gz' -div 'SegmentDilate_blur2_20.nii.gz'  -mas 'csfmask.nii.gz' 'PD_blur2_20_div.nii.gz'
fslmaths  'T2norm_blur2_10.nii.gz' -div 'SegmentDilate_blur2_10.nii.gz'   -mas 'csfmask.nii.gz' 'T2_blur2_10_div.nii.gz'
fslmaths  'T2norm_blur2_20.nii.gz' -div 'SegmentDilate_blur2_20.nii.gz'  -mas 'csfmask.nii.gz' 'T2_blur2_20_div.nii.gz'
fslmaths  'T1norm_blur2_10.nii.gz' -div 'SegmentDilate_blur2_10.nii.gz'   -mas 'csfmask.nii.gz' 'T1_blur2_10_div.nii.gz'
fslmaths  'T1norm_blur2_20.nii.gz' -div 'SegmentDilate_blur2_20.nii.gz'  -mas 'csfmask.nii.gz' 'T1_blur2_20_div.nii.gz'

#######################################################################################################################################
## Do second segmentation (call R)
#######################################################################################################################################

Rscript OASIS_run2.R  'FLAIRnorm.nii.gz' 'PDnorm.nii.gz' 'T2norm.nii.gz' 'T1norm.nii.gz' 'csfmask.nii.gz' 'FLAIR_blur2_10_div.nii.gz' 
'FLAIR_blur2_20_div.nii.gz' 'PD_blur2_10_div.nii.gz' 'PD_blur2_20_div.nii.gz' 'T2_blur2_10_div.nii.gz' 'T2_blur2_20_div.nii.gz' 
'T1_blur2_10_div.nii.gz' 'T1_blur2_20_div.nii.gz' 
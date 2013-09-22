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
## Erodes the brain mask 
#######################################################################################################################################

fslmaths 'DuraMask.nii.gz' -kernel box 5x5x5 -ero 'mask_mask'

#######################################################################################################################################
## Mask the FLAIR 
#######################################################################################################################################

fslmaths 'FLAIR.nii.gz' -mas 'mask_mask' 'flair_strip.nii.gz'

#######################################################################################################################################
## Exclude the bottom 10 percent of the FLAIR (to remove CSF)
#######################################################################################################################################

threshold_flair=`fslstats 'flair_strip.nii.gz' -P 15`
fslmaths  'flair_strip.nii.gz' -thr $threshold_flair -bin 'csfmask.nii.gz' -odt int

#######################################################################################################################################
## Full brain normalization of the images
#######################################################################################################################################

fslmaths 'FLAIR.nii.gz' -mas 'csfmask.nii.gz' 'FLAIRstrip_masked.nii.gz'
fslmaths 'PD.nii.gz' -mas 'csfmask.nii.gz' 'PDstrip_masked.nii.gz'
fslmaths 'T2.nii.gz' -mas 'csfmask.nii.gz' 'T2strip_masked.nii.gz'
fslmaths 'VolumetricT1.nii.gz' -mas 'csfmask.nii.gz' 'T1strip_masked.nii.gz'
mean_flair=`fslstats 'FLAIRstrip_masked.nii.gz' -M`
sd_flair=`fslstats 'FLAIRstrip_masked.nii.gz' -S`
mean_pd=`fslstats 'PDstrip_masked.nii.gz' -M`
sd_pd=`fslstats 'PDstrip_masked.nii.gz' -S`
mean_t2=`fslstats 'T2strip_masked.nii.gz' -M`
sd_t2=`fslstats 'T2strip_masked.nii.gz' -S`
mean_t1=`fslstats 'T1strip_masked.nii.gz' -M`
sd_t1=`fslstats 'T1strip_masked.nii.gz' -S`
fslmaths 'FLAIRstrip_masked.nii.gz' -sub $mean_flair  -div $sd_flair -mas 'csfmask.nii.gz' 'FLAIRnorm.nii.gz'
fslmaths 'PDstrip_masked.nii.gz' -sub $mean_pd  -div $sd_pd -mas 'csfmask.nii.gz' 'PDnorm.nii.gz'
fslmaths 'T2strip_masked.nii.gz' -sub $mean_t2  -div $sd_t2 -mas 'csfmask.nii.gz' 'T2norm.nii.gz'
fslmaths 'T1strip_masked.nii.gz' -sub $mean_t1  -div $sd_t1 -mas 'csfmask.nii.gz' 'T1norm.nii.gz' 

#######################################################################################################################################
## Create first smoothed volumes (kernel window size of 10 and 20)
#######################################################################################################################################

fslmaths  'FLAIRnorm.nii.gz' -mas 'csfmask.nii.gz'  -s 10  'FLAIRnorm_blur1_10.nii.gz'
fslmaths  'FLAIRnorm.nii.gz' -mas 'csfmask.nii.gz' -s 20 'FLAIRnorm_blur1_20.nii.gz'
fslmaths  'PDnorm.nii.gz' -mas 'csfmask.nii.gz' -s 10 'PDnorm_blur1_10.nii.gz'
fslmaths  'PDnorm.nii.gz' -mas 'csfmask.nii.gz' -s 20 'PDnorm_blur1_20.nii.gz'
fslmaths  'T2norm.nii.gz' -mas 'csfmask.nii.gz' -s 10 'T2norm_blur1_10.nii.gz'
fslmaths  'T2norm.nii.gz' -mas 'csfmask.nii.gz' -s 20 'T2norm_blur1_20.nii.gz'
fslmaths  'T1norm.nii.gz' -mas 'csfmask.nii.gz' -s 10 'T1norm_blur1_10.nii.gz'
fslmaths  'T1norm.nii.gz' -mas 'csfmask.nii.gz' -s 20 'T1norm_blur1_20.nii.gz'
fslmaths  'csfmask.nii.gz' -s 10 'csfmask_blur1_10.nii.gz' 
fslmaths  'csfmask.nii.gz' -s 20 'csfmask_blur1_20.nii.gz'  
fslmaths  'FLAIRnorm_blur1_10.nii.gz' -div 'csfmask_blur1_10.nii.gz'   -mas 'csfmask.nii.gz'  'FLAIR_blur1_10_div.nii.gz'
fslmaths  'FLAIRnorm_blur1_20.nii.gz' -div 'csfmask_blur1_20.nii.gz'    -mas 'csfmask.nii.gz'  'FLAIR_blur1_20_div.nii.gz'
fslmaths  'PDnorm_blur1_10.nii.gz' -div 'csfmask_blur1_10.nii.gz'  -mas 'csfmask.nii.gz'  'PD_blur1_10_div.nii.gz'
fslmaths  'PDnorm_blur1_20.nii.gz' -div 'csfmask_blur1_20.nii.gz'    -mas 'csfmask.nii.gz'  'PD_blur1_20_div.nii.gz'
fslmaths  'T2norm_blur1_10.nii.gz' -div 'csfmask_blur1_10.nii.gz'   -mas 'csfmask.nii.gz'  'T2_blur1_10_div.nii.gz'
fslmaths  'T2norm_blur1_20.nii.gz' -div 'csfmask_blur1_20.nii.gz'    -mas 'csfmask.nii.gz'  'T2_blur1_20_div.nii.gz'
fslmaths  'T1norm_blur1_10.nii.gz' -div 'csfmask_blur1_10.nii.gz'   -mas 'csfmask.nii.gz'  'T1_blur1_10_div.nii.gz'
fslmaths  'T1norm_blur1_20.nii.gz' -div 'csfmask_blur1_20.nii.gz'    -mas 'csfmask.nii.gz'  'T1_blur1_20_div.nii.gz'


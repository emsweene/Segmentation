#===============================================================================
# This code creates OASIS segmentations
# Citation: 
# Input:	Name		- File
#		Brain mask	- DuraMask.nii.gz
#		FLAIR		- FLAIRStrip.nii.gz
#		PDw		- PDStrip.nii.gz
#		T2w		- T2Strip.nii.gz
#		T1w		- VolumetricT1Strip.nii.gz
#
# Output:	Name		- File
#
#
# Requirements:
# fsl 
# R (v3.0) with ff, AnalyzeFMRI, ROCR, R.utils
#
# Author: Elizabeth Sweeney, Colin Shea
# Last updated: June 24, 2013
# Change log
# 20130624	- Reformatted comments
#		- Modified paths to use variables
#		- Modified file names to conform to TNU pipeline standard
#===============================================================================
# TODO:
#===============================================================================

#!/bin/bash
CASE="04881278/20120822-05845"
CODE="/Users/deweybe/Desktop/OASIS"
DATA="$CODE/Input/$CASE"
OUTPUT="$CODE/Output/$CASE"


mkdir -p $OUTPUT

## Erodes the brain mask 
fslmaths "$DATA/DuraMask.nii.gz" -kernel box 5x5x5 -ero "$OUTPUT/OASIS_DuraMask"

## Mask the FLAIR 
fslmaths "$DATA/FLAIRStrip.nii.gz" -mas "$OUTPUT/OASIS_DuraMask" "$OUTPUT/OASIS_FLAIRStrip.nii.gz"

## Exclude the bottom 10 percent of the FLAIR (to remove CSF)
threshold_flair=`fslstats "$OUTPUT/OASIS_FLAIRStrip.nii.gz" -P 15`
fslmaths  "$OUTPUT/OASIS_FLAIRStrip.nii.gz" -thr $threshold_flair -bin "$OUTPUT/OASIS_CSFMask.nii.gz" -odt int

## Full brain normalization of the images
fslmaths "$DATA/FLAIRStrip.nii.gz"		-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_FLAIRStrip_masked.nii.gz"
fslmaths "$DATA/PDStrip.nii.gz"			-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_PDStrip_masked.nii.gz"
fslmaths "$DATA/T2Strip.nii.gz"			-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T2Strip_masked.nii.gz"
fslmaths "$DATA/VolumetricT1Strip.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T1Strip_masked.nii.gz"
mean_flair=`fslstats "$OUTPUT/OASIS_FLAIRStrip_masked.nii.gz" -M`
sd_flair=`fslstats "$OUTPUT/OASIS_FLAIRStrip_masked.nii.gz" -S`
mean_pd=`fslstats "$OUTPUT/OASIS_PDStrip_masked.nii.gz" -M`
sd_pd=`fslstats "$OUTPUT/OASIS_PDStrip_masked.nii.gz" -S`
mean_t2=`fslstats "$OUTPUT/OASIS_T2Strip_masked.nii.gz" -M`
sd_t2=`fslstats "$OUTPUT/OASIS_T2Strip_masked.nii.gz" -S`
mean_t1=`fslstats "$OUTPUT/OASIS_T1Strip_masked.nii.gz" -M`
sd_t1=`fslstats "$OUTPUT/OASIS_T1Strip_masked.nii.gz" -S`
fslmaths "$OUTPUT/OASIS_FLAIRStrip_masked.nii.gz"	-sub $mean_flair	-div $sd_flair	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_FLAIRNorm.nii.gz"
fslmaths "$OUTPUT/OASIS_PDStrip_masked.nii.gz"		-sub $mean_pd		-div $sd_pd	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_PDNorm.nii.gz"
fslmaths "$OUTPUT/OASIS_T2Strip_masked.nii.gz"		-sub $mean_t2		-div $sd_t2	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T2Norm.nii.gz"
fslmaths "$OUTPUT/OASIS_T1Strip_masked.nii.gz"		-sub $mean_t1		-div $sd_t1	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T1Norm.nii.gz" 

## Create first smoothed volumes (kernel window size of 10 and 20)
fslmaths  "$OUTPUT/OASIS_FLAIRNorm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 10 "$OUTPUT/OASIS_FLAIRNorm_blur1_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_FLAIRNorm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 20 "$OUTPUT/OASIS_FLAIRNorm_blur1_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 10 "$OUTPUT/OASIS_PDNorm_blur1_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 20 "$OUTPUT/OASIS_PDNorm_blur1_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 10 "$OUTPUT/OASIS_T2Norm_blur1_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 20 "$OUTPUT/OASIS_T2Norm_blur1_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 10 "$OUTPUT/OASIS_T1Norm_blur1_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" -s 20 "$OUTPUT/OASIS_T1Norm_blur1_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_CSFMask.nii.gz" -s 10 "$OUTPUT/OASIS_CSFMask_blur1_10.nii.gz" 
fslmaths  "$OUTPUT/OASIS_CSFMask.nii.gz" -s 20 "$OUTPUT/OASIS_CSFMask_blur1_20.nii.gz"  
fslmaths  "$OUTPUT/OASIS_FLAIRNorm_blur1_10.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_FLAIR_blur1_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_FLAIRNorm_blur1_20.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_FLAIR_blur1_20_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm_blur1_10.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_PD_blur1_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm_blur1_20.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_PD_blur1_20_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm_blur1_10.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_T2_blur1_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm_blur1_20.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_T2_blur1_20_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm_blur1_10.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_T1_blur1_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm_blur1_20.nii.gz" -div "$OUTPUT/OASIS_CSFMask_blur1_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz"  "$OUTPUT/OASIS_T1_blur1_20_div.nii.gz"

# gzip all nii files in Output
gzip $OUTPUT/*.nii

## Do first segmentation (call R)
Rscript "$CODE/OASIS_run1.R"  "$OUTPUT/OASIS_FLAIRNorm.nii.gz" "$OUTPUT/OASIS_PDNorm.nii.gz" "$OUTPUT/OASIS_T2Norm.nii.gz" "$OUTPUT/OASIS_T1Norm.nii.gz" "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_FLAIR_blur1_10_div.nii.gz" "$OUTPUT/OASIS_FLAIR_blur1_20_div.nii.gz" "$OUTPUT/OASIS_PD_blur1_10_div.nii.gz" "$OUTPUT/OASIS_PD_blur1_20_div.nii.gz" "$OUTPUT/OASIS_T2_blur1_10_div.nii.gz" "$OUTPUT/OASIS_T2_blur1_20_div.nii.gz" "$OUTPUT/OASIS_T1_blur1_10_div.nii.gz" "$OUTPUT/OASIS_T1_blur1_20_div.nii.gz"  "$OUTPUT"

gzip $OUTPUT/*.nii

## Dilate first segmentation and apply csf mask
fslmaths "$OUTPUT/OASIS_segment_run1.nii.gz" -kernel boxv 5  -dilM "$OUTPUT/OASIS_segment_dilate.nii.gz"
fslmaths "$OUTPUT/OASIS_CSFMask.nii.gz" -sub "$OUTPUT/OASIS_segment_dilate" -mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz"

gzip $OUTPUT/*.nii

## Create second smoothed volumes (kernel window size of 10 and 20)
fslmaths  "$OUTPUT/OASIS_FLAIRNorm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 10  "$OUTPUT/OASIS_FLAIRNorm_blur2_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_FLAIRNorm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 20 "$OUTPUT/OASIS_FLAIRNorm_blur2_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 10  "$OUTPUT/OASIS_PDNorm_blur2_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 20 "$OUTPUT/OASIS_PDNorm_blur2_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 10  "$OUTPUT/OASIS_T2Norm_blur2_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 20 "$OUTPUT/OASIS_T2Norm_blur2_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 10  "$OUTPUT/OASIS_T1Norm_blur2_10.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm.nii.gz" -mas "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 20 "$OUTPUT/OASIS_T1Norm_blur2_20.nii.gz"
fslmaths  "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 10  "$OUTPUT/OASIS_SegmentDilate_blur2_10.nii.gz" 
fslmaths  "$OUTPUT/OASIS_lesion_segment_dialte.nii.gz" -s 20  "$OUTPUT/OASIS_SegmentDilate_blur2_20.nii.gz"  
fslmaths  "$OUTPUT/OASIS_FLAIRNorm_blur2_10.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_FLAIR_blur2_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_FLAIRNorm_blur2_20.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_FLAIR_blur2_20_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm_blur2_10.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_PD_blur2_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_PDNorm_blur2_20.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_PD_blur2_20_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm_blur2_10.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T2_blur2_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T2Norm_blur2_20.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T2_blur2_20_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm_blur2_10.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_10.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T1_blur2_10_div.nii.gz"
fslmaths  "$OUTPUT/OASIS_T1Norm_blur2_20.nii.gz" -div "$OUTPUT/OASIS_SegmentDilate_blur2_20.nii.gz"	-mas "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_T1_blur2_20_div.nii.gz"

gzip $OUTPUT/*.nii

## Do second segmentation (call R)
Rscript "$CODE/OASIS_run2.R"  "$OUTPUT/OASIS_FLAIRNorm.nii.gz" "$OUTPUT/OASIS_PDNorm.nii.gz" "$OUTPUT/OASIS_T2Norm.nii.gz" "$OUTPUT/OASIS_T1Norm.nii.gz" "$OUTPUT/OASIS_CSFMask.nii.gz" "$OUTPUT/OASIS_FLAIR_blur2_10_div.nii.gz" "$OUTPUT/OASIS_FLAIR_blur2_20_div.nii.gz" "$OUTPUT/OASIS_PD_blur2_10_div.nii.gz" "$OUTPUT/OASIS_PD_blur2_20_div.nii.gz" "$OUTPUT/OASIS_T2_blur2_10_div.nii.gz" "$OUTPUT/OASIS_T2_blur2_20_div.nii.gz" "$OUTPUT/OASIS_T1_blur2_10_div.nii.gz" "$OUTPUT/OASIS_T1_blur2_20_div.nii.gz" "$OUTPUT"

gzip $OUTPUT/*.nii

rm -r $OUTPUT/temp

##########################################################################################################################
##########################################################################################################################
# This file creates probability maps of brain lesions using LIME. The file performs logsitc 
# regression after voxel selection masks  have been made and applied to the brain.  Corrected 
# lesion TOADs masks are used as the Gold Standard. It estimates and saves ROC curves.
#
# Elizabeth Sweeney
# February 14, 2010
##########################################################################################################################
##########################################################################################################################
##CLUSTER commands
##qrsh -l  mem_free=39G,h_vmem=40G
##########################################################################################################################

##This function unzips and reads a .nii.gz file FROM A NODE##

readniigz <- function(file){
	name = gunzip(file, 
		temporary = TRUE, 
		overwrite = TRUE,
		remove = FALSE)
	image <- f.read.nifti.volume(name)
	file.remove(name)
	return(image)
}



candidate_voxels_one_modal <- function(kernel, cutoff, brain_mask, modal) {
	
	dura_mask <-1*(brain_mask>0)[,,,1]	
	sigma.smooth<-diag(3,3)
	
	k.size<-kernel
	blur_modal<-GaussSmoothArray(modal ,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
	brain_modal <- blur_modal[brain_mask == 1]
	
	cutpoint<-quantile(brain_modal, cutoff)
	modal_voxels<-1*(blur_modal>cutpoint)
	
	return(modal_voxels)
}

##load libraries## 
library(AnalyzeFMRI)
library(ROCR)
library(R.utils)
library(pbapply)
library(dplyr)

##SET SEED##
##set.seed(55)


##set data directorys##
hopkins_mri = paste0("/dcl01/smart/data/structural", 
	"/msmri/Hopkins_MRI/")
navid_data = paste0(
	hopkins_mri,
	"Navid_Data")

e_data = file.path(navid_data, "Elizabeth_Processed_MS_Data")
re_data = file.path(e_data, "OASIS_RePro")
##set data directorys##
lesion_dir <- file.path(navid_data, 
	'corrected_lesions/')
lesion_files <- dir(lesion_dir,
	full.names = TRUE)
data_dir <- file.path(e_data, 
	"Register_to_Navid", "NORMALIZED/")

blur_10_dir <- file.path(re_data, 'BLUR/divided/10/')
blur_20_dir <- file.path(re_data, 'BLUR/divided/20/')

#############################
# This is for the doubly-smoothed data
#############################
blur2_10_dir <- file.path(e_data, 'OASIS_final', 
	'BLUR/divided/10/')
blur2_20_dir <- file.path(e_data, 'OASIS_final', 
	'BLUR/divided/20/')

mask_dir <- file.path(re_data, 'Masks/')


h_dir = file.path(hopkins_mri, 
	"Healthy_Brain_Data_process")
flair_dir <- file.path(h_dir, "flair/")
flair_files <- dir(flair_dir,full.names = TRUE)

healthy_dir <- file.path(h_dir, "N3_CoReg_Healthy_with_Spectre/")

dir.exists(lesion_dir)
dir.exists(data_dir)
dir.exists(blur_10_dir)
dir.exists(blur_20_dir)
dir.exists(mask_dir)

dir.exists(blur2_10_dir)
dir.exists(blur2_20_dir)




#####################################################################################################################
#####################################################################################################################
# Train the model
#####################################################################################################################
#####################################################################################################################


#####################################################################################################################
# Reads in the data for the MS subjects
#####################################################################################################################


i = 1

ms_df = vector(mode = "list", length = 98)

for (i in 1:98) {
	
	lesions <- lesion_files[i]
	patient_id <- basename(lesions)
	patient_id <- sub("(.*)_mod.*", "\\1", patient_id)

	#############################################
	##load in files
	#############################################
	gold_lesions <- readniigz(paste0(lesion_dir, 
		patient_id, '_mod.nii.gz'))

	flair <- readniigz(paste0(data_dir, 
		patient_id, '_normalized_flairmni.nii.gz'))

	dura_mask <- readniigz(paste0(mask_dir, 
		patient_id, 'mask.nii.gz'))

	topvoxels <- candidate_voxels_one_modal(1, 
		.85 , dura_mask, flair)

	dura_mask = drop(dura_mask)
	topvoxels = drop(topvoxels)
	dimg = dim(dura_mask)

	df = data_frame(
		topvoxels = c(topvoxels) == 1,
		dura_mask = c(dura_mask) == 1,
		gold_lesions = c(gold_lesions)
		)
	dimg = dim(dura_mask)
	attr(df, "dimg") = dimg
	ms_df[[i]] = df
	print(i)
}

names(ms_df) = paste0("MS_", seq_along(ms_df))


#####################################################################################################################
# Reads in the data for the HEALTHY subjects
#####################################################################################################################
i = 1

healthy_df = vector(mode = "list", length = 33)

for (i in 1:33) {
	
	lesions <- flair_files[i]
	patient_id <- basename(lesions)
	patient_id <- sub("(.*)flair.*", "\\1", patient_id)
	
###################################
##load in files
###################################

	flair = readniigz(
		paste0(healthy_dir, patient_id,
		"FLAIRnorm.nii.gz"))

	csf_fname = paste0(healthy_dir, 
		patient_id, 'csfmask.nii.gz')
	dura_mask <- readniigz(csf_fname)


	topvoxels <- candidate_voxels_one_modal(
		1, .85, 
		dura_mask, flair)

	dura_mask = drop(dura_mask)
	topvoxels = drop(topvoxels)
	dimg = dim(dura_mask)
	df = data_frame(
		topvoxels = c(topvoxels) == 1,
		dura_mask = c(dura_mask) == 1)
	df$gold_lesions = 0
	attr(df, "dimg") = dimg

	healthy_df[[i]] = df
	print(i)
}


names(healthy_df) = paste0("Control_", 
	seq_along(healthy_df))

all_df = c(ms_df, healthy_df)

outfile = file.path(
    hopkins_mri, 
    "full_voxsel_data.rds")

readr::write_rds(all_df, path = outfile)


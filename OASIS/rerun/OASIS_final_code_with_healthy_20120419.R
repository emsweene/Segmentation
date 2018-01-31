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


hopkins_mri = "/dcl01/smart/data/structural/msmri/Hopkins_MRI/"
navid_data = paste0(
	hopkins_mri,
	"Navid_Data")
# navid_data = paste0("/project/taki/msmri/Hopkins_MRI/Navid_Data")

e_data = file.path(navid_data, "Elizabeth_Processed_MS_Data")
re_data = file.path(e_data, "OASIS_RePro")
##set working directory##

##set data directorys##
lesion_dir <- file.path(navid_data, 'corrected_lesions/')
lesion_files <- dir(lesion_dir,full.names = TRUE)
data_dir <- file.path(e_data, "Register_to_Navid", "NORMALIZED/")

blur_10_dir <- file.path(re_data, 'BLUR/divided/10/')
blur_20_dir <- file.path(re_data, 'BLUR/divided/20/')

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




#####################################################################################################################
#####################################################################################################################
# Train the model
#####################################################################################################################
#####################################################################################################################


#####################################################################################################################
# Reads in the data for the MS subjects
#####################################################################################################################

all_df = vector(mode = "list", length = 98)
i = 1

for (i in 1:98) {
	
	print(i)

	lesions <- lesion_files[i]
	patient_id <- basename(lesions)
	patient_id <- sub("(.*)_mod.*", "\\1", patient_id)

	#####################################################################################################################
	##load in files
	#####################################################################################################################
	modes = c("FLAIR", "PD", "T1", "T2")
	suffs = c("mni", "_div_10",
		"_div_20")
	all_imgs = outer(modes, suffs, paste0)
	all_imgs = c(all_imgs)
	n_imgs = all_imgs
	n_imgs = sub("mni", "", n_imgs)
	n_imgs = sub("_div_(\\d\\d)", "_\\1", n_imgs) 
	all_imgs = tolower(all_imgs)
	df = data_frame(patient_id = patient_id,
		img = all_imgs, mod = n_imgs)
	df = df %>% mutate(
		dir = ifelse(grepl("10", mod),
			blur_10_dir, 
			ifelse(grepl("20", mod), 
				blur_20_dir, data_dir)),
		img = ifelse(!grepl("_(10|20)", mod),
			paste0("_normalized_", img), img)
		)
	df = df %>% mutate(
		img = paste0(dir, patient_id, img, ".nii.gz"))
	all_imgs = df$img
	names(all_imgs) = df$mod

	gold_lesions <- readniigz(paste0(lesion_dir, patient_id, '_mod.nii.gz'))
	# flair <- readniigz(paste0(data_dir, patient_id, '_normalized_flairmni.nii.gz'))
	dura_mask <- readniigz(paste0(mask_dir, patient_id, 'mask.nii.gz'))

	all_imgs = pblapply(all_imgs, readniigz)
	flair = all_imgs[["FLAIR"]]

	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	index = which(topvoxels==1)
		
	mat = pbsapply(all_imgs, c)
	mat = mat[index, ]
	mat = as_data_frame(mat)
	mat$GOLD_Lesions = gold_lesions[index]

	all_df[[i]] = mat
	rm(list= c("mat", "topvoxels", "index", 
		"all_imgs", "flair", "dura_mask", "df"))
}


#####################################################################################################################
# Reads in the data for the HEALTHY subjects
#####################################################################################################################
i = 1

healthy_df = vector(mode = "list", length = 33)

for (i in 1:33) {
	
	lesions <- flair_files[i]
	patient_id <- basename(lesions)
	patient_id <- sub("(.*)flair.*", "\\1", patient_id)
	
#####################################################################################################################
##load in files
#####################################################################################################################
	modes = c("FLAIR", "PD", "T1", "T2")
	suffs = c("norm", "_blur1_10_div",
		"_blur1_20_div")
	all_imgs = outer(modes, suffs, paste0)
	all_imgs = c(all_imgs)
	n_imgs = all_imgs
	n_imgs = sub("norm", "", n_imgs)
	n_imgs = sub("_blur1_(\\d\\d)_div", "_\\1", n_imgs) 

	all_imgs = paste0(healthy_dir, patient_id,
		all_imgs, ".nii.gz")
	names(all_imgs) = n_imgs	

	csf_fname = paste0(healthy_dir, patient_id, 'csfmask.nii.gz')
	dura_mask <- readniigz(csf_fname)

	all_imgs = pblapply(all_imgs, readniigz)
	flair = all_imgs[["FLAIR"]]

	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	index = which(topvoxels==1)

	mat = pbsapply(all_imgs, c)
	mat = mat[index, ]
	mat = as_data_frame(mat)
	mat$GOLD_Lesions = 0

	healthy_df[[i]] = mat
	print(i)	
	rm(list= c("mat", "topvoxels", "index", 
		"all_imgs", "flair", "dura_mask"))	
}

all_data = bind_rows(all_df, healthy_df)
print(dim(all_data))

model<- glm(GOLD_Lesions ~ 
	FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
	PD_10 * PD + PD_20 * PD + 
	T2_10 * T2 + T2_20 * T2 + 
	T1_10 * T1 + T1_20 * T1 , 
	data = all_data,
	family = binomial(),
	control= list(trace = TRUE)) 

print(summary(model))

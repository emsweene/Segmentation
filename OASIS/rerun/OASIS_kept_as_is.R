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


library(AnalyzeFMRI)
library(ROCR)
library(R.utils)
library(pbapply)
library(dplyr)

##SET SEED##
##set.seed(55)


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


GOLD_Lesions <- c()
FLAIR <- c()
PD <- c()
T2 <- c()
T1 <- c()

FLAIR_10 <- c()
PD_10 <- c()
T2_10 <- c()
T1_10 <- c()


FLAIR_20 <- c()
PD_20 <- c()
T2_20 <- c()
T1_20 <- c()



#####################################################################################################################
#####################################################################################################################
# Train the model
#####################################################################################################################
#####################################################################################################################


#####################################################################################################################
# Reads in the data for the MS subjects
#####################################################################################################################


for (i in 1:98) {
	print(i)
	lesions <- lesion_files[i]
	patient_id <- unlist(strsplit(lesions, "/"))
	h <- length(patient_id)
	patient_id <- unlist(strsplit(patient_id[h], "_")) 
	patient_id <- patient_id[1]
	
	#####################################################################################################################
	##load in files
	#####################################################################################################################
	

	gold_lesions <- readniigz(paste(lesion_dir, patient_id, '_mod.nii.gz', sep =''))
	flair <- readniigz(paste(data_dir, patient_id, '_normalized_flairmni.nii.gz', sep = ''))
	pd <- readniigz(paste(data_dir, patient_id, '_normalized_pdmni.nii.gz', sep = ''))
	################################### ERROR - switched! ###################################
	t1 <- readniigz(paste(data_dir, patient_id, '_normalized_t2mni.nii.gz', sep = ''))
	t2 <- readniigz(paste(data_dir, patient_id, '_normalized_t1mni.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(mask_dir, patient_id, 'mask.nii.gz', sep = ''))
	
	flair_10 <- readniigz(paste(blur_10_dir, patient_id, 'flair_div_10.nii.gz', sep =''))
	pd_10 <- readniigz(paste(blur_10_dir, patient_id, 'pd_div_10.nii.gz', sep =''))
	t1_10 <- readniigz(paste(blur_10_dir, patient_id, 't1_div_10.nii.gz', sep =''))
	t2_10 <- readniigz(paste(blur_10_dir, patient_id, 't2_div_10.nii.gz', sep =''))
	
	
	flair_20 <- readniigz(paste(blur_20_dir, patient_id, 'flair_div_20.nii.gz', sep =''))
	pd_20 <- readniigz(paste(blur_20_dir, patient_id, 'pd_div_20.nii.gz', sep =''))
	t1_20 <- readniigz(paste(blur_20_dir, patient_id, 't1_div_20.nii.gz', sep =''))
	t2_20 <- readniigz(paste(blur_20_dir, patient_id, 't2_div_20.nii.gz', sep =''))	
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	index = which(topvoxels == 1)

	flair_10 <- flair_10[index]
	pd_10 <- pd_10[index]
	t2_10 <- t2_10[index]
	t1_10 <- t1_10[index]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	PD_10 <- c(PD_10, pd_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(pd_10)
	rm(t2_10)
	rm(t1_10)
	
	
	flair_20 <- flair_20[index]
	pd_20 <- pd_20[index]
	t2_20 <- t2_20[index]
	t1_20 <- t1_20[index]	

	FLAIR_20 <- c(FLAIR_20, flair_20)
	PD_20 <- c(PD_20, pd_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(pd_20)
	rm(t2_20)
	rm(t1_20)

	gold_lesions <- gold_lesions[index]
	flair <- flair[index]
	pd <- pd[index]
	t2 <- t2[index]
	t1 <- t1[index]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	PD <- c(PD, pd)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	
	print(i)
	}


#####################################################################################################################
# Reads in the data for the HEALTHY subjects
#####################################################################################################################


for (i in 1:33) {
	lesions <- flair_files[i]
	patient_id <- unlist(strsplit(lesions, "/"))
	h <- length(patient_id)
	patient_id <- unlist(strsplit(patient_id[h], "f")) 
	patient_id <- patient_id[1]
	
#####################################################################################################################
##load in files
#####################################################################################################################
	

	flair <- readniigz(paste(healthy_dir, patient_id, 'FLAIRnorm.nii.gz', sep = ''))
	pd <- readniigz(paste(healthy_dir, patient_id, 'PDnorm.nii.gz', sep = ''))
	t1 <- readniigz(paste(healthy_dir, patient_id, 'T1norm.nii.gz', sep = ''))
	t2 <- readniigz(paste(healthy_dir, patient_id, 'T2norm.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(healthy_dir, patient_id, 'csfmask.nii.gz', sep = ''))
	
	flair_10 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur1_10_div.nii.gz', sep =''))
	pd_10 <- readniigz(paste(healthy_dir, patient_id, 'PD_blur1_10_div.nii.gz', sep =''))
	t1_10 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur1_10_div.nii.gz', sep =''))
	t2_10 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur1_10_div.nii.gz', sep =''))
	
	
	flair_20 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur1_20_div.nii.gz', sep =''))
	pd_20 <- readniigz(paste(healthy_dir, patient_id, 'PD_blur1_20_div.nii.gz', sep =''))
	t1_20 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur1_20_div.nii.gz', sep =''))
	t2_20 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur1_20_div.nii.gz', sep =''))	
	
	
	gold_lesions <- array(0,dim=dim(flair))
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	index = which(topvoxels == 1)
	
	
	flair_10 <- flair_10[index]
	pd_10 <- pd_10[index]
	t2_10 <- t2_10[index]
	t1_10 <- t1_10[index]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	PD_10 <- c(PD_10, pd_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(pd_10)
	rm(t2_10)
	rm(t1_10)
	
	
	flair_20 <- flair_20[index]
	pd_20 <- pd_20[index]
	t2_20 <- t2_20[index]
	t1_20 <- t1_20[index]	
	
	FLAIR_20 <- c(FLAIR_20, flair_20)
	PD_20 <- c(PD_20, pd_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(pd_20)
	rm(t2_20)
	rm(t1_20)
	
	gold_lesions <- gold_lesions[index]
	flair <- flair[index]
	pd <- pd[index]
	t2 <- t2[index]
	t1 <- t1[index]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	PD <- c(PD, pd)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	print(i)
}

model<-glm(GOLD_Lesions ~ 
	FLAIR_10 *FLAIR +  
	FLAIR_20*FLAIR + 
	PD_10 *PD  + PD_20 *PD + 
	T2_10 * T2 + T2_20 *T2 + 
	T1_10 *T1  + T1_20 *T1 , 
	family=binomial(),
	control= list(trace = TRUE)) 

print(summary(model))


keep_mod = function(model){
  model$y = c()
  model$model = c()
  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$qr$qr = c()  
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()
  attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()
  model
}



original_oasis_model = keep_mod(model)
save(original_oasis_model, 
	file = file.path(hopkins_mri, 
		"original_oasis_model.rda"),
	compress = "xz")

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
	blur2_df = df %>% 
		filter(grepl("_(10|20)", mod))
	blur2_df$dir = sub("OASIS_RePro", "OASIS_final",
		blur2_df$dir)

	df = df %>% mutate(
		img = paste0(dir, patient_id, img, ".nii.gz"))
	blur2_df = blur2_df %>% mutate(
		img = paste0(dir, patient_id, img, ".nii.gz"))
	blur2_df$mod = sub("_", "_2_", blur2_df$mod)

	all_imgs = c(df$img, blur2_df$img)
	names(all_imgs) = c(df$mod, blur2_df$mod)

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
		"_blur1_20_div",
		"_blur2_10_div",
		"_blur2_20_div")
	all_imgs = outer(modes, suffs, paste0)
	all_imgs = c(all_imgs)
	n_imgs = all_imgs
	n_imgs = sub("norm", "", n_imgs)
	n_imgs = sub("_blur1_(\\d\\d)_div", "_\\1", 
		n_imgs) 
	n_imgs = sub("_blur2_(\\d\\d)_div", "_2_\\1", 
		n_imgs) 

	all_imgs = paste0(healthy_dir, patient_id,
		all_imgs, ".nii.gz")
	names(all_imgs) = n_imgs	

	csf_fname = paste0(healthy_dir, 
		patient_id, 'csfmask.nii.gz')
	dura_mask <- readniigz(csf_fname)

	all_imgs = pblapply(all_imgs, readniigz)
	flair = all_imgs[["FLAIR"]]

	topvoxels <- candidate_voxels_one_modal(1, .85 , 
		dura_mask, flair)
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

names(all_df) = paste0("MS_", 
	seq_along(all_df))
names(healthy_df) = paste0("Control_", 
	seq_along(healthy_df))

all_data = bind_rows(all_df, 
	healthy_df, .id = "id")
print(dim(all_data))

###################################
# Final Data
###################################
test_data = all_data[ 
	all_data$id %in% c("MS_1", "Control_1"),]
saveRDS(test_data, 
	file = file.path(hopkins_mri, 
		"test_oasis_data.rds"),
	compress = "xz")



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


###################################
# Correct Full Model
###################################
model<- glm(GOLD_Lesions ~ 
	FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
	PD_10 * PD + PD_20 * PD + 
	T2_10 * T2 + T2_20 * T2 + 
	T1_10 * T1 + T1_20 * T1 , 
	data = all_data,
	family = binomial(),
	control= list(trace = TRUE)) 

print(summary(model))

oasis_model = keep_mod(model)
save(oasis_model, 
	file = file.path(hopkins_mri, "oasis_model.rda"),
	compress = "xz")

####################################
# refinement model from the paper
####################################
model_blur_2 <- glm(GOLD_Lesions ~ 
	FLAIR_2_10 * FLAIR + FLAIR_2_20 * FLAIR + 
	PD_2_10 * PD + PD_2_20 * PD + 
	T2_2_10 * T2 + T2_2_20 * T2 + 
	T1_2_10 * T1 + T1_2_20 * T1, 
	data = all_data, 
	family=binomial(),
	control= list(trace = TRUE)
	) 

print(summary(model_blur_2))
oasis_model_blur2 = keep_mod(model_blur_2)
save(oasis_model_blur2, 
	file = file.path(hopkins_mri, 
    "oasis_model_blur2.rda"),
	compress = "xz")

###################################
# Dropping PD modality
###################################
nopd_oasis_model<- glm(GOLD_Lesions ~ 
	FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
	# PD_10 * PD + PD_20 * PD + 
	T2_10 * T2 + T2_20 * T2 + 
	T1_10 * T1 + T1_20 * T1 , 
	data = all_data,
	family = binomial(),
	control= list(trace = TRUE)) 

print(summary(nopd_oasis_model))

nopd_oasis_model = keep_mod(nopd_oasis_model)
save(nopd_oasis_model, 
	file = file.path(hopkins_mri, "nopd_oasis_model.rda"),
	compress = "xz")


###################################
# Making incorrect data to compare 
# to previous results
###################################
bad_data = all_data
bad_data = bad_data %>% 
	mutate(
		T1_orig = T1,
		T2_orig = T2,
		ms = grepl("MS", id),
		T1 = ifelse(ms, T2_orig, T1_orig),
		T2 = ifelse(ms, T1_orig, T2_orig)
		) %>% 
	select(-T1_orig, -T2_orig, -ms)



original_oasis_model<- glm(GOLD_Lesions ~ 
	FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
	PD_10 * PD + PD_20 * PD + 
	T2_10 * T2 + T2_20 * T2 + 
	T1_10 * T1 + T1_20 * T1 , 
	data = bad_data,
	family = binomial(),
	control= list(trace = TRUE)) 
print(summary(original_oasis_model))

original_oasis_model = keep_mod(original_oasis_model)
save(original_oasis_model, 
	file = file.path(hopkins_mri, 
		"original_oasis_model.rda"),
	compress = "xz")


####################################
# refinement model from the paper
####################################
original_model_blur_2 <- glm(GOLD_Lesions ~ 
	FLAIR_2_10 * FLAIR + FLAIR_2_20 * FLAIR + 
	PD_2_10 * PD + PD_2_20 * PD + 
	T2_2_10 * T2 + T2_2_20 * T2 + 
	T1_2_10 * T1 + T1_2_20 * T1, 
	data = bad_data, 
	family=binomial(),
	control= list(trace = TRUE)
	) 

print(summary(original_model_blur_2))
original_oasis_model_blur2 = keep_mod(original_model_blur_2)
save(original_oasis_model_blur2, 
	file = file.path(hopkins_mri, 
    "original_oasis_model_blur2.rda"),
	compress = "xz")



###################################
# Final Data
###################################
saveRDS(all_data, 
	file = file.path(hopkins_mri, 
		"correct_oasis_data.rds"),
	compress = "xz")

# original data
saveRDS(bad_data, 
	file = file.path(hopkins_mri, 
		"original_oasis_data.rds"),
	compress = "xz")



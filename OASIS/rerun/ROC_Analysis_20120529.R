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
	if (!file.exists('/scratch/temp/esweene')) dir.create('/scratch/temp/esweene')
	name <- tempfile(pattern = "tempfile", tmpdir = '/scratch/temp/esweene', fileext = ".nii")
	file.copy(file,paste(name,'.gz',sep=''),overwrite=TRUE)
	gunzip(paste(name,'.gz',sep=''),name,overwrite=TRUE)
	image <- f.read.nifti.volume(name)
	try(unlink(name))
	return(image)
}


##This function writes and zipts a .nii.gz file FROM A NODE##

writeniigz <- function(image,file){
	if (!file.exists('/scratch/temp/esweene')) dir.create('/scratch/temp/esweene')
	name <- tempfile(pattern = "tempfile", tmpdir = '/scratch/temp/esweene', fileext = ".nii")
	f.write.nifti(image,name,size="float",nii=TRUE)
	gzip(name, paste(name,'.gz',sep=''), remove = TRUE,overwrite=TRUE)
	try.save<-file.copy(paste(name,'.gz',sep=''),file,overwrite=TRUE)
	try(unlink(paste(name,'.gz',sep='')))
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



bigglm.ffdf<-function (formula, data, ..., chunksize = 100000) 
{
    n <- nrow(data)
    cursor <- 0
    datafun <- function(reset = FALSE) {
        if (reset) {
            cursor <<- 0
            return(NULL)
        }
        if (cursor >= n) 
		return(NULL)
        start <- cursor + 1
        cursor <<- cursor + min(chunksize, n - cursor)
        data[start:cursor, ]
    }
    rval <- bigglm(formula = formula, data = datafun, ...)
    rval
}



library(ff)
library(AnalyzeFMRI)
library(ROCR)
library(R.utils)
library(biglm)
library(randomForest)
library(MASS)

##SET SEED##
set.seed(33)

##set working directory##



##set data directorys##
lesion_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/corrected_lesions/')
lesion_files <- dir(lesion_dir,full.names = TRUE)

flair_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Healthy_Brain_Data_process/flair/')
flair_files <- dir(flair_dir,full.names = TRUE)

data_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/Register_to_Navid/NORMALIZED/')
healthy_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Healthy_Brain_Data_process/N3_CoReg_Healthy_with_Spectre/')
blur_10_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/BLUR/divided/10/')
blur_20_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/BLUR/divided/20/')

mask_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/Masks/')



n <- 98

##set training and validation set##
m <- 49
training <- sample(1:n, m, replace = F)
sets <- c(1:n)
valid <-(sets[-training])


h_n <- 33
h_m <- 16
h_training <- sample(1:h_n, h_m, replace = F)
sets <- c(1:h_n)
h_valid <-(sets[-h_training])






load('the_models')






GOLD <- c() 
PREDICTION_FOREST <- c()
PREDICTION_LOGISTIC <- c()
PREDICTION_LDA <- c()

l <- n - m 

for (i in 1:l) {
	
	k <- valid[i]
	
	lesions <- lesion_files[k]
	patient_id <- unlist(strsplit(lesions, "/"))
	h <- length(patient_id)
	patient_id <- unlist(strsplit(patient_id[h], "_")) 
	patient_id <- patient_id[1]
	
#####################################################################################################################
##load in files
#####################################################################################################################
	
	gold_lesions <- readniigz(paste(lesion_dir, patient_id, '_mod.nii.gz', sep =''))
	FLAIR <- readniigz(paste(data_dir, patient_id, '_normalized_flairmni.nii.gz', sep = ''))
	PD <- readniigz(paste(data_dir, patient_id, '_normalized_pdmni.nii.gz', sep = ''))
	T1 <- readniigz(paste(data_dir, patient_id, '_normalized_t2mni.nii.gz', sep = ''))
	T2 <- readniigz(paste(data_dir, patient_id, '_normalized_t1mni.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(mask_dir, patient_id, 'mask.nii.gz', sep = ''))
	
	FLAIR_10 <- readniigz(paste(blur_10_dir, patient_id, 'flair_div_10.nii.gz', sep =''))
	PD_10 <- readniigz(paste(blur_10_dir, patient_id, 'pd_div_10.nii.gz', sep =''))
	T1_10 <- readniigz(paste(blur_10_dir, patient_id, 't1_div_10.nii.gz', sep =''))
	T2_10 <- readniigz(paste(blur_10_dir, patient_id, 't2_div_10.nii.gz', sep =''))
	
	
	FLAIR_20 <- readniigz(paste(blur_20_dir, patient_id, 'flair_div_20.nii.gz', sep =''))
	PD_20 <- readniigz(paste(blur_20_dir, patient_id, 'pd_div_20.nii.gz', sep =''))
	T1_20 <- readniigz(paste(blur_20_dir, patient_id, 't1_div_20.nii.gz', sep =''))
	T2_20 <- readniigz(paste(blur_20_dir, patient_id, 't2_div_20.nii.gz', sep =''))	
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, FLAIR)
	
#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	OASIS <- as.data.frame(cbind(FLAIR, PD, T1, T2, FLAIR_10, PD_10, T1_10, T2_10, FLAIR_20, PD_20, T1_20, T2_20))
	prob_map <- predict(forest_fit , newdata = OASIS, type = c("response"))
	
	prob_map <- array(prob_map,dim=dim(FLAIR))
	prob_map[topvoxels==0]<-0
	
#####################################################################################################################	
##Smooth the image
#####################################################################################################################
	
	sigma.smooth<-diag(3,3)
	k.size<- 5
	prob_map<-GaussSmoothArray(prob_map,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
####################################################################################################################
##write probability map##
#####################################################################################################################

#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	prob_map_1 <- predict(logistic_model, newdata = OASIS, type = c("response"))
	
	prob_map_1 <- array(prob_map_1,dim=dim(FLAIR))
	prob_map_1[topvoxels==0]<-0
	
#####################################################################################################################	
##Smooth the image
#####################################################################################################################
	
	sigma.smooth<-diag(3,3)
	k.size<- 5
	prob_map_1<-GaussSmoothArray(prob_map_1,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
####################################################################################################################
##write probability map##
#####################################################################################################################	
	
	

#####################################################################################################################
##create probability maps
#####################################################################################################################
	prob_map_2 <- predict(lda_fit, newdata = OASIS)
	
	prob_map_2 <- array(as.numeric(prob_map_2$posterior[,2]),dim=dim(FLAIR))
	prob_map_2[topvoxels==0]<-0
	
#####################################################################################################################	
##Smooth the image
#####################################################################################################################
	
	sigma.smooth<-diag(3,3)
	k.size<- 5
	prob_map_2<-GaussSmoothArray(prob_map_2,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
####################################################################################################################
##write probability map##
#####################################################################################################################	
		
	filename <- paste(patient_id, "forest" , sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
	filename <- paste(patient_id, "logisitc" , sep = "_")
	f.write.nifti(prob_map_1, filename ,size="float",nii=TRUE)
	filename <- paste(patient_id, "lda" , sep = "_")
	f.write.nifti(prob_map_2, filename ,size="float",nii=TRUE)
		
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	prob_map_1 <- prob_map_1[dura_mask == 1]
	prob_map_2 <- prob_map_2[dura_mask == 1]	
	
	GOLD <- c(GOLD, gold_lesions) 
	PREDICTION_FOREST <- c(PREDICTION_FOREST, prob_map)
	PREDICTION_LOGISTIC <- c(PREDICTION_LOGISTIC, prob_map_1)
	PREDICTION_LDA <- c(PREDICTION_LDA, prob_map_2)
	
	print(i)
}


save.image('just_MS')

l <- h_n - h_m

for (i in 1:l) {
	
	k <- h_valid[i]
	
	lesions <- flair_files[k]
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
	
		topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, FLAIR)
	
#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	OASIS <- as.data.frame(cbind(FLAIR, PD, T1, T2, FLAIR_10, PD_10, T1_10, T2_10, FLAIR_20, PD_20, T1_20, T2_20))
	prob_map <- predict(forest_fit , newdata = OASIS, type = c("response"))
	
	prob_map <- array(prob_map,dim=dim(FLAIR))
	prob_map[topvoxels==0]<-0
	
#####################################################################################################################	
##Smooth the image
#####################################################################################################################
	
	sigma.smooth<-diag(3,3)
	k.size<- 5
	prob_map<-GaussSmoothArray(prob_map,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
####################################################################################################################
##write probability map##
#####################################################################################################################

#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	prob_map_1 <- predict(logistic_model, newdata = OASIS, type = c("response"))
	
	prob_map_1 <- array(prob_map_1,dim=dim(FLAIR))
	prob_map_1[topvoxels==0]<-0
	
#####################################################################################################################	
##Smooth the image
#####################################################################################################################
	
	sigma.smooth<-diag(3,3)
	k.size<- 5
	prob_map_1<-GaussSmoothArray(prob_map_1,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
####################################################################################################################
##write probability map##
#####################################################################################################################	
	
	

#####################################################################################################################
##create probability maps
#####################################################################################################################
	prob_map_2 <- predict(lda_fit, newdata = OASIS)
	
	prob_map_2 <- array(as.numeric(prob_map_2$posterior[,2]),dim=dim(FLAIR))
	prob_map_2[topvoxels==0]<-0
	
#####################################################################################################################	
##Smooth the image
#####################################################################################################################
	
	sigma.smooth<-diag(3,3)
	k.size<- 5
	prob_map_2<-GaussSmoothArray(prob_map_2,sigma=sigma.smooth,ksize=k.size,mask=dura_mask)
	
####################################################################################################################
##write probability map##
#####################################################################################################################	
		
	filename <- paste(patient_id, "forest_healthy" , sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
	filename <- paste(patient_id, "logisitc_healthy" , sep = "_")
	f.write.nifti(prob_map_1, filename ,size="float",nii=TRUE)
	filename <- paste(patient_id, "lda_healthy" , sep = "_")
	f.write.nifti(prob_map_2, filename ,size="float",nii=TRUE)
		
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	prob_map_1 <- prob_map_1[dura_mask == 1]
	prob_map_2 <- prob_map_2[dura_mask == 1]	
	
	GOLD <- c(GOLD, gold_lesions) 
	PREDICTION_FOREST <- c(PREDICTION_FOREST, prob_map)
	PREDICTION_LOGISTIC <- c(PREDICTION_LOGISTIC, prob_map_1)
	PREDICTION_LDA <- c(PREDICTION_LDA, prob_map_2)
	
	print(i)	
}

save.image('all')
save(GOLD, PREDICTION_LDA, PREDICTION_FOREST, PREDICTION_LOGISTIC, file= 'ROC_stuffs')

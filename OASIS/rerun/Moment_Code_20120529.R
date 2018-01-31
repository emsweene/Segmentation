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

##This function calculates the local moments of an image## 

local_moment <- function(image, window, moment) {
	
	moment_image <- image^moment
	moment <- array(0, dim = dim(moment_image)[1:3])
	indices <- permutations(window, 3, v= (-window/2 + .5):(window/2 - .5), repeats.allowed=TRUE)	
	
	for ( i in 1:dim(indices)[1]){
		shifter <- ashift(moment_image, indices[i,])
		moment <-moment + shifter
	}
 	
	moment <- moment / dim(indices)[1]
	moment <- array(moment, dim = dim(image))
	return(moment)
}


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
library(magic)

##SET SEED##
set.seed(33)

##set working directory##
setwd('/dexter/disk1/smart/msmri/Hopkins_MRI/OASIS_ML/')


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





GOLD_Lesions <- c()
FLAIR <- c()
T2 <- c()
T1 <- c()

FLAIR_10 <- c()
T2_10 <- c()
T1_10 <- c()

FLAIR_20 <- c()
T2_20 <- c()
T1_20 <- c()

FLAIR_FIRST_3 <- c()
FLAIR_SECOND_3 <- c()
FLAIR_THIRD_3 <- c()

T2_FIRST_3 <- c()
T2_SECOND_3 <- c()
T2_THIRD_3 <- c()

T1_FIRST_3 <- c()
T1_SECOND_3 <- c()
T1_THIRD_3 <- c()


#####################################################################################################################
#####################################################################################################################
# Train the model
#####################################################################################################################
#####################################################################################################################


#####################################################################################################################
# Reads in the data for the MS subjects
#####################################################################################################################


for (i in 1:m) {
	
	k <- training[i]
	
	lesions <- lesion_files[k]
	patient_id <- unlist(strsplit(lesions, "/"))
	h <- length(patient_id)
	patient_id <- unlist(strsplit(patient_id[h], "_")) 
	patient_id <- patient_id[1]
	
	#####################################################################################################################
	##load in files
	#####################################################################################################################
	

	gold_lesions <- readniigz(paste(lesion_dir, patient_id, '_mod.nii.gz', sep =''))
	flair <- readniigz(paste(data_dir, patient_id, '_normalized_flairmni.nii.gz', sep = ''))
	t1 <- readniigz(paste(data_dir, patient_id, '_normalized_t2mni.nii.gz', sep = ''))
	t2 <- readniigz(paste(data_dir, patient_id, '_normalized_t1mni.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(mask_dir, patient_id, 'mask.nii.gz', sep = ''))
	
	flair_10 <- readniigz(paste(blur_10_dir, patient_id, 'flair_div_10.nii.gz', sep =''))
	t1_10 <- readniigz(paste(blur_10_dir, patient_id, 't1_div_10.nii.gz', sep =''))
	t2_10 <- readniigz(paste(blur_10_dir, patient_id, 't2_div_10.nii.gz', sep =''))
	
	flair_20 <- readniigz(paste(blur_20_dir, patient_id, 'flair_div_20.nii.gz', sep =''))
	t1_20 <- readniigz(paste(blur_20_dir, patient_id, 't1_div_20.nii.gz', sep =''))
	t2_20 <- readniigz(paste(blur_20_dir, patient_id, 't2_div_20.nii.gz', sep =''))	
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	

	flair_10 <- flair_10[topvoxels == 1]
	t2_10 <- t2_10[topvoxels == 1]
	t1_10 <- t1_10[topvoxels == 1]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(t2_10)
	rm(t1_10)
	
	flair_20 <- flair_20[topvoxels == 1]
	t2_20 <- t2_20[topvoxels == 1]
	t1_20 <- t1_20[topvoxels == 1]	

	FLAIR_20 <- c(FLAIR_20, flair_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(t2_20)
	rm(t1_20)

	flair_first_3 <- local_moment(flair, 3, 1)
	flair_first_3 <- flair_first_3[topvoxels == 1]
	flair_second_3 <- local_moment(flair, 3, 2)
	flair_second_3 <- flair_second_3[topvoxels == 1]
	flair_third_3 <- local_moment(flair, 3, 3)
	flair_third_3 <- flair_third_3[topvoxels == 1]
	
	t1_first_3 <- local_moment(t1, 3, 1)
	t1_first_3 <- t1_first_3[topvoxels == 1]
	t1_second_3 <- local_moment(t1, 3, 2)
	t1_second_3 <- t1_second_3[topvoxels == 1]
	t1_third_3 <- local_moment(t1, 3, 3)
	t1_third_3 <- t1_third_3[topvoxels == 1]
	
	t2_first_3 <- local_moment(t2, 3, 1)
	t2_first_3 <- t2_first_3[topvoxels == 1]
	t2_second_3 <- local_moment(t2, 3, 2)
	t2_second_3 <- t2_second_3[topvoxels == 1] 
	t2_third_3 <- local_moment(t2, 3, 3)
	t2_third_3 <- t2_third_3[topvoxels == 1]
	
	FLAIR_FIRST_3 <- c(FLAIR_FIRST_3, flair_first_3)
	FLAIR_SECOND_3 <- c(FLAIR_SECOND_3, flair_second_3)
	FLAIR_THIRD_3 <- c(FLAIR_THIRD_3, flair_third_3)
	
	T2_FIRST_3 <- c(T2_FIRST_3, t2_first_3)
	T2_SECOND_3 <- c(T2_SECOND_3, t2_second_3)
	T2_THIRD_3 <- c(T2_THIRD_3, t2_third_3)
	
	T1_FIRST_3 <- c(T1_FIRST_3, t1_first_3)
	T1_SECOND_3 <- c(T1_SECOND_3, t1_second_3)
	T1_THIRD_3 <- c(T1_THIRD_3, t1_third_3)
	
	gold_lesions <- gold_lesions[topvoxels == 1]
	flair <- flair[topvoxels == 1]
	t2 <- t2[topvoxels == 1]
	t1 <- t1[topvoxels == 1]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	
	print(i)
	print(patient_id)
	}


for (i in 1:h_m) {
	
	k <- h_training[i]
	
	lesions <- flair_files[k]
	patient_id <- unlist(strsplit(lesions, "/"))
	h <- length(patient_id)
	patient_id <- unlist(strsplit(patient_id[h], "f")) 
	patient_id <- patient_id[1]
	
#####################################################################################################################
##load in files
#####################################################################################################################

	flair <- readniigz(paste(healthy_dir, patient_id, 'FLAIRnorm.nii.gz', sep = ''))
	t1 <- readniigz(paste(healthy_dir, patient_id, 'T1norm.nii.gz', sep = ''))
	t2 <- readniigz(paste(healthy_dir, patient_id, 'T2norm.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(healthy_dir, patient_id, 'csfmask.nii.gz', sep = ''))
	
	flair_10 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur1_10_div.nii.gz', sep =''))
	t1_10 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur1_10_div.nii.gz', sep =''))
	t2_10 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur1_10_div.nii.gz', sep =''))
	
	
	flair_20 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur1_20_div.nii.gz', sep =''))
	t1_20 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur1_20_div.nii.gz', sep =''))
	t2_20 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur1_20_div.nii.gz', sep =''))	
	
	
	gold_lesions <- array(0,dim=dim(flair))
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	
	flair_first_3 <- local_moment(flair, 3, 1)
	flair_first_3 <- flair_first_3[topvoxels == 1]
	flair_second_3 <- local_moment(flair, 3, 2)
	flair_second_3 <- flair_second_3[topvoxels == 1]
	flair_third_3 <- local_moment(flair, 3, 3)
	flair_third_3 <- flair_third_3[topvoxels == 1]
	
	t1_first_3 <- local_moment(t1, 3, 1)
	t1_first_3 <- t1_first_3[topvoxels == 1]
	t1_second_3 <- local_moment(t1, 3, 2)
	t1_second_3 <- t1_second_3[topvoxels == 1]
	t1_third_3 <- local_moment(t1, 3, 3)
	t1_third_3 <- t1_third_3[topvoxels == 1]
	
	t2_first_3 <- local_moment(t2, 3, 1)
	t2_first_3 <- t2_first_3[topvoxels == 1]
	t2_second_3 <- local_moment(t2, 3, 2)
	t2_second_3 <- t2_second_3[topvoxels == 1] 
	t2_third_3 <- local_moment(t2, 3, 3)
	t2_third_3 <- t2_third_3[topvoxels == 1]
	
	FLAIR_FIRST_3 <- c(FLAIR_FIRST_3, flair_first_3)
	FLAIR_SECOND_3 <- c(FLAIR_SECOND_3, flair_second_3)
	FLAIR_THIRD_3 <- c(FLAIR_THIRD_3, flair_third_3)
	
	T2_FIRST_3 <- c(T2_FIRST_3, t2_first_3)
	T2_SECOND_3 <- c(T2_SECOND_3, t2_second_3)
	T2_THIRD_3 <- c(T2_THIRD_3, t2_third_3)
	
	T1_FIRST_3 <- c(T1_FIRST_3, t1_first_3)
	T1_SECOND_3 <- c(T1_SECOND_3, t1_second_3)
	T1_THIRD_3 <- c(T1_THIRD_3, t1_third_3)	
	
	flair_10 <- flair_10[topvoxels == 1]
	t2_10 <- t2_10[topvoxels == 1]
	t1_10 <- t1_10[topvoxels == 1]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(t2_10)
	rm(t1_10)
	
	flair_20 <- flair_20[topvoxels == 1]
	t2_20 <- t2_20[topvoxels == 1]
	t1_20 <- t1_20[topvoxels == 1]	
	
	FLAIR_20 <- c(FLAIR_20, flair_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(t2_20)
	rm(t1_20)

	gold_lesions <- gold_lesions[topvoxels == 1]
	flair <- flair[topvoxels == 1]
	t2 <- t2[topvoxels == 1]
	t1 <- t1[topvoxels == 1]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)

print(i)
	
}

save.image('train_vector')


# Logistic Regression Model
logistic_model<-glm(GOLD_Lesions ~  FLAIR_10 *FLAIR +  FLAIR_20*FLAIR + FLAIR_FIRST_3*FLAIR + FLAIR_SECOND_3*FLAIR + FLAIR_THIRD_3*FLAIR
+ T2_10 *T2  + T2_20 *T2 + T2_FIRST_3*T2 + T2_SECOND_3*T2 + T2_THIRD_3*T2
+ T1_10 *T1  + T1_20 *T1 + T1_FIRST_3*T1 + T1_SECOND_3*T1 + T1_THIRD_3*T1, family=binomial()) 

save(logistic_model, file='logistic_model_with_moments.RData')



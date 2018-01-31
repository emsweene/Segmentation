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


library(AnalyzeFMRI)
library(ROCR)
library(R.utils)

##SET SEED##
##set.seed(55)


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
	

	flair_10 <- flair_10[topvoxels == 1]
	pd_10 <- pd_10[topvoxels == 1]
	t2_10 <- t2_10[topvoxels == 1]
	t1_10 <- t1_10[topvoxels == 1]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	PD_10 <- c(PD_10, pd_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(pd_10)
	rm(t2_10)
	rm(t1_10)
	
	
	flair_20 <- flair_20[topvoxels == 1]
	pd_20 <- pd_20[topvoxels == 1]
	t2_20 <- t2_20[topvoxels == 1]
	t1_20 <- t1_20[topvoxels == 1]	

	FLAIR_20 <- c(FLAIR_20, flair_20)
	PD_20 <- c(PD_20, pd_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(pd_20)
	rm(t2_20)
	rm(t1_20)

	gold_lesions <- gold_lesions[topvoxels == 1]
	flair <- flair[topvoxels == 1]
	pd <- pd[topvoxels == 1]
	t2 <- t2[topvoxels == 1]
	t1 <- t1[topvoxels == 1]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	PD <- c(PD, pd)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	
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
	
	
	flair_10 <- flair_10[topvoxels == 1]
	pd_10 <- pd_10[topvoxels == 1]
	t2_10 <- t2_10[topvoxels == 1]
	t1_10 <- t1_10[topvoxels == 1]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	PD_10 <- c(PD_10, pd_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(pd_10)
	rm(t2_10)
	rm(t1_10)
	
	
	flair_20 <- flair_20[topvoxels == 1]
	pd_20 <- pd_20[topvoxels == 1]
	t2_20 <- t2_20[topvoxels == 1]
	t1_20 <- t1_20[topvoxels == 1]	
	
	FLAIR_20 <- c(FLAIR_20, flair_20)
	PD_20 <- c(PD_20, pd_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(pd_20)
	rm(t2_20)
	rm(t1_20)
	
	gold_lesions <- gold_lesions[topvoxels == 1]
	flair <- flair[topvoxels == 1]
	pd <- pd[topvoxels == 1]
	t2 <- t2[topvoxels == 1]
	t1 <- t1[topvoxels == 1]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	PD <- c(PD, pd)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	
}

model<-glm(GOLD_Lesions ~  FLAIR_10 *FLAIR +  FLAIR_20*FLAIR + PD_10 *PD  + PD_20 *PD 
+ T2_10 *T2  + T2_20 *T2 + T1_10 *T1  + T1_20 *T1 , family=binomial()) 


b0 <- coef(model)[1]
b1 <-  coef(model)[2]
b2 <-  coef(model)[3]
b3 <- coef(model)[4]
b4 <-  coef(model)[5]	
b5 <-  coef(model)[6]
b6 <-  coef(model)[7]
b7 <- coef(model)[8]
b8 <-  coef(model)[9]	
b9 <- coef(model)[10]
b10 <-  coef(model)[11]
b11 <-  coef(model)[12]
b12 <- coef(model)[13]
b13 <-  coef(model)[14]	
b14 <-  coef(model)[15]
b15 <-  coef(model)[16]
b16 <- coef(model)[17]
b17 <-  coef(model)[18]	
b18 <- coef(model)[19]
b19 <-  coef(model)[20]	
b20 <-  coef(model)[21]

print(summary(model))

save(b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, 
file='coefficient_1_output.RData')

GOLD <- c()
PREDICTION <- c()


for (i in 1:98) {
	
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
	
#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	model_prediction <- b0 + b1 *flair_10 + b2 * flair +  b3 * flair_20 + 
	b4 * pd_10 + b5 * pd + b6 * pd_20 + 
	b7 * t2_10 + b8 * t2 + b9 * t2_20 + 
	b10 * t1_10 + b11 * t1 + b12 * t1_20 + 
	b13 * flair_10 * flair + b14 * flair * flair_20
	b15 * pd_10 * pd + b16 * pd * pd_20
	b17 * t2_10 * t2 + b18 * t2 * t2_20
	b19 * t1_10 * t1 + b20 * t1 * t1_20
	
#####################################################################################################################
##transform predictions into probailities##
#####################################################################################################################
	
	prob_map <- (exp(model_prediction)/(1 + exp(model_prediction)))
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
	filename <- paste(patient_id, "prob_map_1",  sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
	
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	
	GOLD <- c(GOLD, gold_lesions) 
	PREDICTION <- c(PREDICTION, prob_map)
	
}




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
	
	
#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	model_prediction <- b0 + b1 *flair_10 + b2 * flair +  b3 * flair_20 + 
	b4 * pd_10 + b5 * pd + b6 * pd_20 + 
	b7 * t2_10 + b8 * t2 + b9 * t2_20 + 
	b10 * t1_10 + b11 * t1 + b12 * t1_20 + 
	b13 * flair_10 * flair + b14 * flair * flair_20
	b15 * pd_10 * pd + b16 * pd * pd_20
	b17 * t2_10 * t2 + b18 * t2 * t2_20
	b19 * t1_10 * t1 + b20 * t1 * t1_20
	
	
#####################################################################################################################
##transform predictions into probailities##
#####################################################################################################################
	
	prob_map <- (exp(model_prediction)/(1 + exp(model_prediction)))
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
	filename <- paste(patient_id, "prob_map_1",  sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
	
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	
	GOLD <- c(GOLD, gold_lesions) 
	PREDICTION <- c(PREDICTION, prob_map)
}



voxel_level_prediction <-prediction(PREDICTION , GOLD)
voxel_level_performance<-performance(voxel_level_prediction, "tpr",  "fpr")
auc_voxel<-as.numeric(performance(voxel_level_prediction, measure = "auc")@y.values)


y <- as.numeric(voxel_level_performance@y.values[[1]])
x <-  as.numeric(voxel_level_performance@x.values[[1]])
alpha <- as.numeric(voxel_level_performance@alpha.values[[1]])

Blur_x_values <- c()
Blur_y_values <- c()
Blur_alpha_values <- c()

fpr_stop <- .01

n <- length(x)


n <- length(x)
r <- n/1000

for (i in 1:r){
	if (x[i*1000] < fpr_stop) 
	{ 
		Blur_x_values[i] <- x[i*1000]
		Blur_y_values[i] <- y[i*1000]
		Blur_alpha_values[i] <- alpha[i*1000]
	}
	else{
		stop
	}
}


movie_alpha <- min(Blur_alpha_values)
sensitivity <- max(Blur_y_values)
specificity <- max(Blur_x_values)

print(movie_alpha)
print(sensitivity)
print(specificity)

save(movie_alpha, sensitivity, specificity, file= 'first_cut_output.RData')

# all this shit is on doubly - blurred data

lesion_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/corrected_lesions/')
lesion_files <- dir(lesion_dir,full.names = TRUE)

flair_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Healthy_Brain_Data_process/flair/')
flair_files <- dir(flair_dir,full.names = TRUE)

data_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/Register_to_Navid/NORMALIZED/')
healthy_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Healthy_Brain_Data_process/N3_CoReg_Healthy_with_Spectre/')

blur_10_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_final/BLUR/divided/10/')
blur_20_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_final/BLUR/divided/20/')

mask_dir <- c('/dexter/disk1/smart/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/Masks/')

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
# Reads in the data, creates and applies voxel selection masks to images, puts intensities into ff object
#####################################################################################################################


for (i in 1:98) {
	

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
	

	flair_10 <- flair_10[topvoxels == 1]
	pd_10 <- pd_10[topvoxels == 1]
	t2_10 <- t2_10[topvoxels == 1]
	t1_10 <- t1_10[topvoxels == 1]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	PD_10 <- c(PD_10, pd_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(pd_10)
	rm(t2_10)
	rm(t1_10)
	
	
	flair_20 <- flair_20[topvoxels == 1]
	pd_20 <- pd_20[topvoxels == 1]
	t2_20 <- t2_20[topvoxels == 1]
	t1_20 <- t1_20[topvoxels == 1]	

	FLAIR_20 <- c(FLAIR_20, flair_20)
	PD_20 <- c(PD_20, pd_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(pd_20)
	rm(t2_20)
	rm(t1_20)

	gold_lesions <- gold_lesions[topvoxels == 1]
	flair <- flair[topvoxels == 1]
	pd <- pd[topvoxels == 1]
	t2 <- t2[topvoxels == 1]
	t1 <- t1[topvoxels == 1]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	PD <- c(PD, pd)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	
	}



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
	
	flair_10 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur2_10_div.nii.gz', sep =''))
	pd_10 <- readniigz(paste(healthy_dir, patient_id, 'PD_blur2_10_div.nii.gz', sep =''))
	t1_10 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur2_10_div.nii.gz', sep =''))
	t2_10 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur2_10_div.nii.gz', sep =''))
	
	
	flair_20 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur2_20_div.nii.gz', sep =''))
	pd_20 <- readniigz(paste(healthy_dir, patient_id, 'PD_blur2_20_div.nii.gz', sep =''))
	t1_20 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur2_20_div.nii.gz', sep =''))
	t2_20 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur2_20_div.nii.gz', sep =''))	
	
	
	gold_lesions <- array(0,dim=dim(flair))
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	
	
	flair_10 <- flair_10[topvoxels == 1]
	pd_10 <- pd_10[topvoxels == 1]
	t2_10 <- t2_10[topvoxels == 1]
	t1_10 <- t1_10[topvoxels == 1]
	
	FLAIR_10 <- c(FLAIR_10, flair_10)
	PD_10 <- c(PD_10, pd_10)
	T2_10 <- c(T2_10, t2_10)
	T1_10 <- c(T1_10, t1_10)
	
	rm(flair_10)
	rm(pd_10)
	rm(t2_10)
	rm(t1_10)
	
	
	flair_20 <- flair_20[topvoxels == 1]
	pd_20 <- pd_20[topvoxels == 1]
	t2_20 <- t2_20[topvoxels == 1]
	t1_20 <- t1_20[topvoxels == 1]	
	
	FLAIR_20 <- c(FLAIR_20, flair_20)
	PD_20 <- c(PD_20, pd_20)
	T2_20 <- c(T2_20, t2_20)
	T1_20 <- c(T1_20, t1_20)
	
	rm(flair_20)
	rm(pd_20)
	rm(t2_20)
	rm(t1_20)
	
	gold_lesions <- gold_lesions[topvoxels == 1]
	flair <- flair[topvoxels == 1]
	pd <- pd[topvoxels == 1]
	t2 <- t2[topvoxels == 1]
	t1 <- t1[topvoxels == 1]
	
	GOLD_Lesions <- c(GOLD_Lesions,gold_lesions)
	FLAIR <- c(FLAIR, flair)
	PD <- c(PD, pd)
	T2 <- c(T2, t2)
	T1 <- c(T1, t1)
	
}

model<-glm(GOLD_Lesions ~  FLAIR_10 *FLAIR +  FLAIR_20*FLAIR + PD_10 *PD  + PD_20 *PD 
+ T2_10 *T2  + T2_20 *T2 + T1_10 *T1  + T1_20 *T1 , family=binomial()) 

print(summary(model))


	

#####################################################################################################################
##enter the coefficients from the regression##
#####################################################################################################################


b0_2 <- coef(model)[1]
b1_2 <-  coef(model)[2]
b2_2 <-  coef(model)[3]
b3_2 <- coef(model)[4]
b4_2 <-  coef(model)[5]	
b5_2 <-  coef(model)[6]
b6_2 <-  coef(model)[7]
b7_2 <- coef(model)[8]
b8_2 <-  coef(model)[9]	
b9_2 <- coef(model)[10]
b10_2 <-  coef(model)[11]
b11_2 <-  coef(model)[12]
b12_2 <- coef(model)[13]
b13_2 <-  coef(model)[14]	
b14_2 <-  coef(model)[15]
b15_2 <-  coef(model)[16]
b16_2 <- coef(model)[17]
b17_2 <-  coef(model)[18]	
b18_2 <- coef(model)[19]
b19_2 <-  coef(model)[20]	
b20_2 <-  coef(model)[21]


save(b0_2, b1_2, b2_2, b3_2, b4_2, b5_2, b6_2, b7_2, b8_2, b9_2, b10_2, b11_2, b12_2, b13_2, b14_2, b15_2, 
b16_2, b17_2, b18_2, b19_2, b20_2, file='coefficient_2_output.RData')

GOLD <- c()
PREDICTION <- c()


for (i in 1:98) {
	
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
	
#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	model_prediction <- b0_2 + b1_2 *flair_10 + b2_2 * flair +  b3_2 * flair_20 + 
	b4_2 * pd_10 + b5_2 * pd + b6_2 * pd_20 + 
	b7_2 * t2_10 + b8_2 * t2 + b9_2 * t2_20 + 
	b10_2 * t1_10 + b11_2 * t1 + b12_2 * t1_20 + 
	b13_2 * flair_10 * flair + b14_2 * flair * flair_20
	b15_2 * pd_10 * pd + b16_2 * pd * pd_20
	b17_2 * t2_10 * t2 + b18_2 * t2 * t2_20
	b19_2 * t1_10 * t1 + b20_2 * t1 * t1_20
	
#####################################################################################################################
##transform predictions into probailities##
#####################################################################################################################
	
	prob_map <- (exp(model_prediction)/(1 + exp(model_prediction)))
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
	filename <- paste(patient_id, "prob_map_2",  sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
	
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	
	GOLD <- c(GOLD, gold_lesions) 
	PREDICTION <- c(PREDICTION, prob_map)
	
}




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
	
	flair_10 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur2_10_div.nii.gz', sep =''))
	pd_10 <- readniigz(paste(healthy_dir, patient_id, 'PD_blur2_10_div.nii.gz', sep =''))
	t1_10 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur2_10_div.nii.gz', sep =''))
	t2_10 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur2_10_div.nii.gz', sep =''))
	
	
	flair_20 <- readniigz(paste(healthy_dir, patient_id, 'FLAIR_blur2_20_div.nii.gz', sep =''))
	pd_20 <- readniigz(paste(healthy_dir, patient_id, 'PD_blur2_20_div.nii.gz', sep =''))
	t1_20 <- readniigz(paste(healthy_dir, patient_id, 'T1_blur2_20_div.nii.gz', sep =''))
	t2_20 <- readniigz(paste(healthy_dir, patient_id, 'T2_blur2_20_div.nii.gz', sep =''))	
	
	
	gold_lesions <- array(0,dim=dim(flair))
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
	
	
#####################################################################################################################
##create probability maps
#####################################################################################################################
	
	model_prediction <- b0_2 + b1_2 *flair_10 + b2_2 * flair +  b3_2 * flair_20 + 
	b4_2 * pd_10 + b5_2 * pd + b6_2 * pd_20 + 
	b7_2 * t2_10 + b8_2 * t2 + b9_2 * t2_20 + 
	b10_2 * t1_10 + b11_2 * t1 + b12_2 * t1_20 + 
	b13_2 * flair_10 * flair + b14_2 * flair * flair_20
	b15_2 * pd_10 * pd + b16_2 * pd * pd_20
	b17_2 * t2_10 * t2 + b18_2 * t2 * t2_20
	b19_2 * t1_10 * t1 + b20_2 * t1 * t1_20
	
	
#####################################################################################################################
##transform predictions into probailities##
#####################################################################################################################
	
	prob_map <- (exp(model_prediction)/(1 + exp(model_prediction)))
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
	filename <- paste(patient_id, "prob_map_2",  sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
	
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	
	GOLD <- c(GOLD, gold_lesions) 
	PREDICTION <- c(PREDICTION, prob_map)

}


voxel_level_prediction <-prediction(PREDICTION , GOLD)
voxel_level_performance<-performance(voxel_level_prediction, "tpr",  "fpr")
auc_voxel_2<-as.numeric(performance(voxel_level_prediction, measure = "auc")@y.values)


y <- as.numeric(voxel_level_performance@y.values[[1]])
x <-  as.numeric(voxel_level_performance@x.values[[1]])
alpha <- as.numeric(voxel_level_performance@alpha.values[[1]])

Blur_x_values <- c()
Blur_y_values <- c()
Blur_alpha_values <- c()

fpr_stop <- .005

n <- length(x)


n <- length(x)
r <- n/1000

for (i in 1:r){
	if (x[i*1000] < fpr_stop) 
	{ 
		Blur_x_values[i] <- x[i*1000]
		Blur_y_values[i] <- y[i*1000]
		Blur_alpha_values[i] <- alpha[i*1000]
	}
	else{
		stop
	}
}


movie_alpha_005 <- min(Blur_alpha_values)
sensitivity_005 <- max(Blur_y_values)
specificity_005 <- max(Blur_x_values)

print(movie_alpha_005)
print(sensitivity_005)
print(specificity_005)


Blur_x_values <- c()
Blur_y_values <- c()
Blur_alpha_values <- c()

fpr_stop <- .0025                

n <- length(x)

n <- length(x)
r <- n/1000

for (i in 1:r){
	if (x[i*1000] < fpr_stop) 
	{ 
		Blur_x_values[i] <- x[i*1000]
		Blur_y_values[i] <- y[i*1000]
		Blur_alpha_values[i] <- alpha[i*1000]
	}
	else{
		stop
	}
}


movie_alpha_0025 <- min(Blur_alpha_values)
sensitivity_0025 <- max(Blur_y_values)
specificity_0025 <- max(Blur_x_values)

print(movie_alpha_0025)
print(sensitivity_0025)
print(specificity_0025)





Blur_x_values <- c()
Blur_y_values <- c()
Blur_alpha_values <- c()

fpr_stop <- .0075                

n <- length(x)

n <- length(x)
r <- n/1000

for (i in 1:r){
	if (x[i*1000] < fpr_stop) 
	{ 
		Blur_x_values[i] <- x[i*1000]
		Blur_y_values[i] <- y[i*1000]
		Blur_alpha_values[i] <- alpha[i*1000]
	}
	else{
		stop
	}
}


movie_alpha_0075 <- min(Blur_alpha_values)
sensitivity_0075 <- max(Blur_y_values)
specificity_0075 <- max(Blur_x_values)

print(movie_alpha_0075)
print(sensitivity_0075)
print(specificity_0075)




Blur_x_values <- c()
Blur_y_values <- c()
Blur_alpha_values <- c()

fpr_stop <- .01                

n <- length(x)

n <- length(x)
r <- n/1000

for (i in 1:r){
	if (x[i*1000] < fpr_stop) 
	{ 
		Blur_x_values[i] <- x[i*1000]
		Blur_y_values[i] <- y[i*1000]
		Blur_alpha_values[i] <- alpha[i*1000]
	}
	else{
		stop
	}
}


movie_alpha_01 <- min(Blur_alpha_values)
sensitivity_01 <- max(Blur_y_values)
specificity_01 <- max(Blur_x_values)

print(movie_alpha_01)
print(sensitivity_01)
print(specificity_01)



save(movie_alpha_01, sensitivity_01, specificity_01, movie_alpha_005, sensitivity_005 ,specificity_005, movie_alpha_0025, sensitivity_0025,
specificity_0025, movie_alpha_0075, sensitivity_0075, specificity_0075, file='cutoff_2_output.RData')
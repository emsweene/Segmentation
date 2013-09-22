##########################################################################################################################
##########################################################################################################################
# This file creates probability maps of brain lesions using OASIS. The file performs logsitc 
# regression after voxel selection masks  have been made and applied to the brain.  Manual  
# segmentations are used as the Gold Standard. It estimates and saves ROC curves.
#
# Elizabeth Sweeney
# November 1, 2011
##########################################################################################################################
##########################################################################################################################
##CLUSTER commands
##qrsh -l  mem_free=19G,h_vmem=20G
##########################################################################################################################


#####################################################################################################################
#####################################################################################################################
## Fit Model a Second Time
#####################################################################################################################
#####################################################################################################################


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


##load libraries##
library(AnalyzeFMRI)
library(ROCR)
library(R.utils)


##SET SEED##
set.seed(55)

##set working directory##
setwd('/nexsan/bst5/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/OASIS_maps_run1/')


##set data directorys##
lesion_dir <- c('/nexsan/bst5/msmri/Hopkins_MRI/Navid_Data/corrected_lesions/')
lesion_files <- dir(lesion_dir,full.names = TRUE)
data_dir <- c('/nexsan/bst5/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/Register_to_Navid/NORMALIZED/')

blur_10_dir <- c('/nexsan/bst5/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/BLUR_2/divided/10/')
blur_20_dir <- c('/nexsan/bst5/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/BLUR_2/divided/20/')

mask_dir <- c('/nexsan/bst5/msmri/Hopkins_MRI/Navid_Data/Elizabeth_Processed_MS_Data/OASIS_RePro/Masks/')


n <- 98

##set training and validation set##
m <- 5
training <- sample(1:n, m, replace = F)
sets <- c(1:n)
valid <-(sets[-training])



GOLD_Lesions <- c()
FLAIR <- c()
PD <- c()
T2 <- c()
T1 <- c()

FLAIR_10 <- c()
PD_10 <- c()
T2_10 <- c()
T1_10 <- c()

FLAIR_15 <- c()
PD_15 <- c()
T2_15 <- c()
T1_15 <- c()

FLAIR_20 <- c()
PD_20 <- c()
T2_20 <- c()
T1_20 <- c()

#####################################################################################################################
#####################################################################################################################
# Train the model
#####################################################################################################################
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
	pd <- readniigz(paste(data_dir, patient_id, '_normalized_pdmni.nii.gz', sep = ''))
	t1 <- readniigz(paste(data_dir, patient_id, '_normalized_t2mni.nii.gz', sep = ''))
	t2 <- readniigz(paste(data_dir, patient_id, '_normalized_t1mni.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(mask_dir, patient_id, 'mask.nii.gz', sep = ''))
	
	flair_10 <- readniigz(paste(blur_10_dir, patient_id, 'flair_div_10.nii.gz', sep =''))
	pd_10 <- readniigz(paste(blur_10_dir, patient_id, 'pd_div_10.nii.gz', sep =''))
	t1_10 <- readniigz(paste(blur_10_dir, patient_id, 't1_div_10.nii.gz', sep =''))
	t2_10 <- readniigz(paste(blur_10_dir, patient_id, 't2_div_10.nii.gz', sep =''))
	
	flair_15 <- readniigz(paste(blur_15_dir, patient_id, 'flair_div_15.nii.gz', sep =''))
	pd_15 <- readniigz(paste(blur_15_dir, patient_id, 'pd_div_15.nii.gz', sep =''))
	t1_15 <- readniigz(paste(blur_15_dir, patient_id, 't1_div_15.nii.gz', sep =''))
	t2_15 <- readniigz(paste(blur_15_dir, patient_id, 't2_div_15.nii.gz', sep =''))
	
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
	
	flair_15 <- flair_15[topvoxels == 1]
	pd_15 <- pd_15[topvoxels == 1]
	t2_15 <- t2_15[topvoxels == 1]
	t1_15 <- t1_15[topvoxels == 1]
	
	FLAIR_15 <- c(FLAIR_15, flair_15)
	PD_15 <- c(PD_15, pd_15)
	T2_15 <- c(T2_15, t2_15)
	T1_15 <- c(T1_15, t1_15)
	
	rm(flair_15)
	rm(pd_15)
	rm(t2_15)
	rm(t1_15)
	
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
	
	print(i)
	
}

save.image('train_set')


model<-glm(GOLD_Lesions ~  FLAIR_10 *FLAIR  + FLAIR_15 *FLAIR +  FLAIR_20*FLAIR + PD_10 *PD  + PD_15 *PD + PD_20 *PD
+ T2_10 *T2 + T2_15 *T2 + T2_20 *T2 + T1_10 *T1  + T1_15 *T1  + T1_20 *T1 , family=binomial())

summary(model)

save(model,file="summary_2")	


model.table <- xtable(model)
print(model.table)


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
b21 <-  coef(model)[22]
b22 <- coef(model)[23]
b23 <-  coef(model)[24]
b24 <- coef(model)[25]
b25 <-  coef(model)[26]
b26 <-  coef(model)[27]
b27 <-  coef(model)[28]
b28 <- coef(model)[29]



GOLD <- c()
PREDICTION <- c()

for (i in 1:(n-m)) {
	
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
	flair <- readniigz(paste(data_dir, patient_id, '_normalized_flairmni.nii.gz', sep = ''))
	pd <- readniigz(paste(data_dir, patient_id, '_normalized_pdmni.nii.gz', sep = ''))
	t1 <- readniigz(paste(data_dir, patient_id, '_normalized_t2mni.nii.gz', sep = ''))
	t2 <- readniigz(paste(data_dir, patient_id, '_normalized_t1mni.nii.gz', sep = ''))
	dura_mask <- readniigz(paste(mask_dir, patient_id, 'mask.nii.gz', sep = ''))
	
	flair_10 <- readniigz(paste(blur_10_dir, patient_id, 'flair_div_10.nii.gz', sep =''))
	pd_10 <- readniigz(paste(blur_10_dir, patient_id, 'pd_div_10.nii.gz', sep =''))
	t1_10 <- readniigz(paste(blur_10_dir, patient_id, 't1_div_10.nii.gz', sep =''))
	t2_10 <- readniigz(paste(blur_10_dir, patient_id, 't2_div_10.nii.gz', sep =''))
	
	flair_15 <- readniigz(paste(blur_15_dir, patient_id, 'flair_div_15.nii.gz', sep =''))
	pd_15 <- readniigz(paste(blur_15_dir, patient_id, 'pd_div_15.nii.gz', sep =''))
	t1_15 <- readniigz(paste(blur_15_dir, patient_id, 't1_div_15.nii.gz', sep =''))
	t2_15 <- readniigz(paste(blur_15_dir, patient_id, 't2_div_15.nii.gz', sep =''))
	
	flair_20 <- readniigz(paste(blur_20_dir, patient_id, 'flair_div_20.nii.gz', sep =''))
	pd_20 <- readniigz(paste(blur_20_dir, patient_id, 'pd_div_20.nii.gz', sep =''))
	t1_20 <- readniigz(paste(blur_20_dir, patient_id, 't1_div_20.nii.gz', sep =''))
	t2_20 <- readniigz(paste(blur_20_dir, patient_id, 't2_div_20.nii.gz', sep =''))
	
	topvoxels <- candidate_voxels_one_modal(1, .85 , dura_mask, flair)
    
	#####################################################################################################################
	##create probability maps
	#####################################################################################################################
	
	model_prediction <- b0 + b1 *flair_10 + b2 * flair + b3 * flair_15 + b4 * flair_20 +
    b5 * pd_10 + b6 * pd + b7 * pd_15 + b8 * pd_20 +
    b9 * t2_10 + b10 * t2 + b11 * t2_15 + b12 * t2_20 +
    b13 * t1_10 + b14 * t1 + b15 * t1_15 + b16 * t1_20 +
    b17 * flair_10 * flair + b18 * flair * flair_15 + b19 * flair * flair_20
    b20 * pd_10 * pd + b21 * pd * pd_15 + b22 * pd * pd_20
    b23 * t2_10 * t2 + b24 * t2 * t2_15 + b25 * t2 * t2_20
    b26 * t1_10 * t1 + b27 * t1 * t1_15 + b28 * t1 * t1_20
	
	
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
	
	filename <- paste(patient_id, "prob_map" , sep = "_")
	f.write.nifti(prob_map, filename ,size="float",nii=TRUE)
    
	gold_lesions <- gold_lesions[dura_mask == 1]
	prob_map <- prob_map[dura_mask == 1]
	
	GOLD <- c(GOLD, gold_lesions)
	PREDICTION <- c(PREDICTION, prob_map)
	
	print(i)
	print(patient_id)
}

save.image('valid_set_final')




voxel_level_prediction <-prediction(PREDICTION , GOLD)
voxel_level_performance<-performance(voxel_level_prediction, "tpr",  "fpr")
auc_voxel<-as.numeric(performance(voxel_level_prediction, measure = "auc")@y.values)


pdf('pROC_OASIS.pdf')

logi_x <- as.numeric(voxel_level_performance@x.values[[1]])
sample_it <- seq(1, length(logi_x), by = length(logi_x)/1000)
logi_x <- logi_x[sample_it]
logi_y <- as.numeric(voxel_level_performance@y.values[[1]])[sort(sample_it)]

plot(logi_y ~ logi_x, col = 'red' , main = "Partial ROC Curve", xlab = "False Positive Rate", ylab = "True Positive Rate",
xlim = c(0,.01), ylim = c(0,1), lwd = 3, lty = 1, bg = 338)

dev.off()




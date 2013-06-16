#!/usr/bin/env Rscript 

temp <- commandArgs(TRUE)

flair.filename<- temp[1] 
pd.filename <- temp[2]
t2.filename <- temp[3]
t1.filename <- temp[4]
csf.mask <- temp[5]
flair.blur_10 <- temp[6]
flair.blur_20 <- temp[7]
pd.blur_10 <- temp[8]
pd.blur_20 <- temp[9]
t2.blur_10 <- temp[10]
t2.blur_20 <- temp[11]
t1.blur_10 <- temp[12]
t1.blur_20 <- temp[13]



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



library(ff)
library(AnalyzeFMRI)
library(ROCR)
library(R.utils)


   
b0 <-  -6.03230   
b1 <-  -0.64586    
b2 <-  2.57678 
b3 <-  8.92298   
b4 <-  -0.50567 
b5 <- -0.36789   
b6 <- -1.85971   
b7 <- -7.43754    
b8 <- 1.09178   
b9 <-  13.39033    
b10 <-  6.63801    
b11 <- 1.78357    
b12 <- -13.95689   
b13 <-  2.43494   
b14 <-  -13.82208  
b15 <- -1.75647  
b16 <- 1.40212  
b17 <- 2.56332   
b18 <- -6.42123   
b19 <- 0.01067    
b20 <-  -1.27630 


## fpr rate of not sure ##
movie_alpha_1 <-  0.16
movie_alpha_2 <-  0.23
movie_alpha_3 <-  0.30
movie_alpha_4 <-  0.35
movie_alpha_5 <-  0.40
	

#####################################################################################################################
##load in files
#####################################################################################################################
flair <- readniigz(flair.filename)
pd <-  readniigz(pd.filename)
t2 <- readniigz(t2.filename)
t1 <- readniigz(t1.filename)
dura_mask <- readniigz(csf.mask)

flair_10 <- readniigz(flair.blur_10)
pd_10 <- readniigz(pd.blur_10)
t2_10 <- readniigz(t2.blur_10)
t1_10 <- readniigz(t1.blur_10)

flair_20 <- readniigz(flair.blur_20)
pd_20 <- readniigz(pd.blur_20)
t2_20 <- readniigz(t2.blur_20)
t1_20 <- readniigz(t1.blur_20)
	
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
	
	prob_map_1 <- prob_map
	prob_map_2 <- prob_map
	prob_map_3 <- prob_map
	prob_map_4 <- prob_map
	prob_map_5 <- prob_map
	
####################################################################################################################
##write probability map##
#####################################################################################################################
	
##binary segmentation of the image
	prob_map_1[prob_map_1  > movie_alpha_1] <- 1
	prob_map_1[prob_map_1  <= movie_alpha_1] <- 0	
	
    filename <- c('_segment_run2_1')
	f.write.nifti(prob_map_1, filename ,size="float",nii=TRUE)

	prob_map_2[prob_map_2  > movie_alpha_2] <- 1
	prob_map_2[prob_map_2  <= movie_alpha_2] <- 0	

	filename <- c('_segment_run2_2')
	f.write.nifti(prob_map_2, filename ,size="float",nii=TRUE)

	prob_map_3[prob_map_3  > movie_alpha_3] <- 1
	prob_map_3[prob_map_3  <= movie_alpha_3] <- 0	

	filename <- c('_segment_run2_3')
	f.write.nifti(prob_map_3, filename ,size="float",nii=TRUE)

	prob_map_4[prob_map_4  > movie_alpha_4] <- 1
	prob_map_4[prob_map_4  <= movie_alpha_4] <- 0	

	filename <- c('_segment_run2_4')
	f.write.nifti(prob_map_4, filename ,size="float",nii=TRUE)	

	prob_map_5[prob_map_5  > movie_alpha_5] <- 1
	prob_map_5[prob_map_5  <= movie_alpha_5] <- 0	

	filename <- c('_segment_run2_5')
	f.write.nifti(prob_map_5, filename ,size="float",nii=TRUE)


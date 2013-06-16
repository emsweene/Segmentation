##########################################################################################################################
##########################################################################################################################
#
# Elizabeth Sweeney
# May 7, 2012
##########################################################################################################################
##########################################################################################################################


##This function unzips and reads a .nii.gz file.##
unzipit <- function(file){
	name <- tempfile(pattern = "file", tmpdir = '/dexter/disk1/smart/msmri/SuBLIME_Data_2/TEMP/', fileext = ".nii")
	gunzip(file, name, remove = FALSE)
	image <- f.read.nifti.volume(name)
	return(image)
}

##This function select slices from images (slices 50, 70, 90, 110 and 120) and returns the intesnity from the selected## 
##slices as a column vector##
sliceit <- function(file) {
	for(k in 0:4){
		x <- cbind(file[1,,20*k+50,])
		for(i in 2:182){
			y <- cbind(file[i,,20*k+50,])
			y <- rbind(x,y)
			x <- y
		}
		assign(paste("column",k , sep = ""), x)		
	}
	z <- rbind(column0, column1, column2, column3, column4)
	return(z)
}


#####################################################################################################################
#####################################################################################################################
##Preperation to Read in the Data
#####################################################################################################################


library(ff)
library(biglm)
library(AnalyzeFMRI)
library(R.utils)
library(ROCR)

##SET SEED##
set.seed(33)

##set working directory##
setwd('/dexter/disk1/smart/msmri/SuBLIME_Data_2/Incidence_Masks_20120507/')

##set data directory##
dir <- c('/dexter/disk1/smart/msmri/SuBLIME_Data_2/')

##Subject Numbers##
subject <-  c("01","02","09", "10")
n <- length(subject)

##Create file directory for each patient##
fileDir <- rep(0, n)
for (i in 1:n) {
	directory <- paste(dir, subject[i], sep = "")
	fileDir[i] <- directory
}


b0 <-  -7.7420 
b1 <-  0.7412 
b2 <-  0.4099  
b3 <-  -0.3226 
b4 <-  0.7807 
b5 <-  0.1841  
b6 <- 0.5383  
b7 <- 0.8546 
b8 <-  -0.9016 


###fpr rate = .0005
## sensitivity 
## 0.7824324
## movie_alpha <-  0.04941387


###fpr rate = .00025
## sensitivity 
## 0.6885135
## move_alpha <- 0.07949942


####fpr rate = .0001
##sensitivity
## 0.4898649
##movie_alpha <-  0.1503217


 movie_alpha <- 0.07949942

#################################################################################################################################
## Make SuBLIME segmentations
#################################################################################################################################

for(i in 1:4) {
	
	data.dir <- data.dir <- fileDir[i]
	files <- dir(data.dir,full.names = TRUE)
	k <- length(files) - 1
	
	for(j in 1:k) {
		data.path1 <- files[j]
		
#####################################################################################################################
## calculates the time between scans
######################################################################################################################
		
		data.path2 <- files[j+1]
		
##get date from scan j+1##
		get.date <- unlist(strsplit(data.path2, "-"))
		get.date <- unlist(strsplit(get.date[1], "/")) 
		n <- length(get.date)
		date2 <- as.numeric(get.date[n])
				

		
#####################################################################################################################
## reads FLAIR, PD, T2, and T1 images from time point j and j+1
#####################################################################################################################
		
		FLAIR_norm_1<-unzipit(paste(data.path1, 'FLAIRnorm.nii.gz', sep = '/'))
		PD_norm_1 <- unzipit(paste(data.path1, 'PDnorm.nii.gz', sep = '/'))
		T2_norm_1 <- unzipit(paste(data.path1, 'T2norm.nii.gz', sep = '/'))
		T1_norm_1 <- unzipit(paste(data.path1, 'VolumetricT1norm.nii.gz', sep = '/'))
		
		FLAIR_norm_2 <-unzipit(paste(data.path2, 'FLAIRnorm.nii.gz', sep = '/'))
		PD_norm_2 <- unzipit(paste(data.path2, 'PDnorm.nii.gz', sep = '/'))
		T2_norm_2 <- unzipit(paste(data.path2, 'T2norm.nii.gz', sep = '/'))
		T1_norm_2 <- unzipit(paste(data.path2, 'VolumetricT1norm.nii.gz', sep = '/'))
		
#####################################################################################################################
## reads brain mask from time point j and j+1
######################################################################################################################
		
		fov <- PD_norm_1[1,1,1,]
		PD_1 <- PD_norm_1
		print(fov)
		PD_1[PD_1 > fov] <- 1
		PD_1[PD_1 <= fov] <- 0
		
		fov <- PD_norm_2[1,1,1,]
		PD_2 <- PD_norm_2
		print(fov)
		PD_2[PD_2 > fov] <- 1
		PD_2[PD_2 <= fov] <- 0
		
		fov <- FLAIR_norm_1[1,1,1,]
		FLAIR_1 <- FLAIR_norm_1
		print(fov)
		FLAIR_1[FLAIR_1 > fov] <- 1
		FLAIR_1[FLAIR_1 <= fov] <- 0
		
		fov <- FLAIR_norm_2[1,1,1,]
		FLAIR_2 <- FLAIR_norm_2
		print(fov)
		FLAIR_2[FLAIR_2 > fov] <- 1
		FLAIR_2[FLAIR_2 <= fov] <- 0
		
		fov <- T1_norm_1[1,1,1,]
		T1_1 <- T1_norm_1
		print(fov)
		T1_1[T1_1 > fov] <- 1
		T1_1[T1_1 <= fov] <- 0
		
		fov <- T1_norm_2[1,1,1,]
		T1_2 <- T1_norm_2
		print(fov)
		T1_2[T1_2 > fov] <- 1
		T1_2[T1_2 <= fov] <- 0
		
		fov <- T2_norm_1[1,1,1,]
		T2_1 <- T2_norm_1
		print(fov)
		T2_1[T2_1 > fov] <- 1
		T2_1[T2_1 <= fov] <- 0
		
		fov <- T2_norm_2[1,1,1,]
		T2_2 <- T2_norm_2
		print(fov)
		T2_2[T2_2 > fov] <- 1
		T2_2[T2_2 <= fov] <- 0
		
		FOV_mask <- PD_1 * PD_2 * FLAIR_1 * FLAIR_2 * T1_1 * T1_2 * T2_1 * T2_2
		Class2 <-unzipit(paste(data.path2, 'Class2.nii.gz', sep = '/'))
		Class2[Class2 == 5] <- 0
		Class2[Class2 == 14] <- 0
		Class2[Class2 == 1] <- 0
		Class2[Class2 > 0] <- 1
		
		mask_brain <- FOV_mask * Class2
##		filename <- paste('fov', subject[i], date2  , sep = "_")
##		f.write.nifti(mask_brain, filename ,size="float",nii=TRUE)
		
#####################################################################################################################
## voxel selection method
######################################################################################################################
		
		
		
##difference of the images##
		diff1 <- T2_norm_2 - T2_norm_1
		diff2 <- FLAIR_norm_2 - FLAIR_norm_1
		diff3 <- PD_norm_2 - PD_norm_1
		
##Smooth the image##
		brain.mask.01<-1*(mask_brain>0)[,,,1]
		sigma.smooth<-diag(3,3)
		k.size<-5
		diff1 <-GaussSmoothArray(diff1 ,sigma=sigma.smooth,ksize=k.size,mask=brain.mask.01)
		diff2 <-GaussSmoothArray(diff2 ,sigma=sigma.smooth,ksize=k.size,mask=brain.mask.01)
		diff3 <-GaussSmoothArray(diff3 ,sigma=sigma.smooth,ksize=k.size,mask=brain.mask.01)
		
##convert to binary mask, cutoff point c##
		c1 <- c(sd(diff1))
		diff1[diff1 > c1] <- 1
		diff1[diff1 <= c1] <- 0
		
		c2 <- c(sd(diff2))
		diff2[diff2 > c2] <- 1
		diff2[diff2 <= c2] <- 0
		
		c3 <- c(sd(diff3))
		diff3[diff3 > c3] <- 1
		diff3[diff3 <= c3] <- 0
		
		quant.diff <- diff1 * diff2 * diff3
		voxel_selection <- mask_brain * quant.diff
		
		
#####################################################################################################################
## Predictions from the SuBLIME model
######################################################################################################################
		
		FLAIR_norm_diff <- FLAIR_norm_2 - FLAIR_norm_1
		PD_norm_diff <- PD_norm_2 - PD_norm_1
		T2_norm_diff <- T2_norm_2 - T2_norm_1
		T1_norm_diff <- T1_norm_2 - T1_norm_1
		
		model_prediction_1 <- b0 + b1*FLAIR_norm_2 + b2*PD_norm_2 + b3*T2_norm_2 + b4*T1_norm_2 + b5*FLAIR_norm_diff + 
		b6*PD_norm_diff + b7*T2_norm_diff + b8*T1_norm_diff
		
		prob_map_1 <- (exp(model_prediction_1)/(1 + exp(model_prediction_1)))
		prob_map_1[voxel_selection ==0]<-0
		
		
##Smooth the image
		
		brain.mask.01<-1*(mask_brain>0)[,,,1]
		sigma.smooth<-diag(3,3)
		k.size<- 3
		
		prob_map_1<-GaussSmoothArray(prob_map_1,sigma=sigma.smooth,ksize=k.size,mask=brain.mask.01)

		
##write probability map##
		print(j)
		
		filename <- paste(subject[i], date2  , 'prob_map' , sep = "_")
		f.write.nifti(prob_map_1, filename ,size="float",nii=TRUE)
	}

	
		
}

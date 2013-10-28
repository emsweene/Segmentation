library(rgl)
library(car)

setwd('/Users/ane/Dropbox/Elizabeth_Sweeney_Documents/Current_Projects/ML/Visualization_20130814/')

load('OASIS_vectors_5_subj_all_visualization_seed_33.RData')

##########################################################################################################################################
##Create masks of voxels neighboring the lesion in fsl (dilate) and use these
##to plot voxels in the neighborhood of a lesion
##########################################################################################################################################


NEIGHBOR_5 <- NEIGHBOR_5 - NEIGHBOR_3
NEIGHBOR_3 <- NEIGHBOR_3 - GOLD_Lesions  



##########################################################################################################################################
##Sample *insert your favorite value below* voxels for the plots
##########################################################################################################################################

set.seed(33)
n <- 10000
MS_subjects <- as.numeric(SUBJECT[duplicated(SUBJECT) == FALSE])
x <- c()
numbers <- seq(1:length(SUBJECT))

##########################################################################################################################################
##Sample *insert your favorite value below* voxels from each subject 
##########################################################################################################################################

for(i in 1:5){
  y <- SUBJECT == MS_subjects[i]
  x1 <- numbers[y]
  x <- c(x, sample(x1, n, replace = FALSE))
}

##########################################################################################################################################
##Downsample original vectors
##########################################################################################################################################

S_GOLD_Lesions <- GOLD_Lesions[x]
S_TOPVOXELS <- TOPVOXELS[x]
S_NO_VS_NO_NORM_FLAIR <- NO_VS_NO_NORM_FLAIR[x]
S_NO_VS_NO_NORM_T2 <- NO_VS_NO_NORM_T2[x]
S_NO_VS_NO_NORM_T1 <- NO_VS_NO_NORM_T1[x]
S_NO_VS_NORM_FLAIR <- NO_VS_NORM_FLAIR[x]
S_NO_VS_NORM_T2 <- NO_VS_NORM_T2[x]
S_NO_VS_NORM_T1 <- NO_VS_NORM_T1[x]
S_NO_VS_FLAIR_10 <- NO_VS_FLAIR_10[x]
S_NO_VS_T2_10 <- NO_VS_T2_10[x]
S_NO_VS_T1_10 <- NO_VS_T1_10[x]
S_NO_VS_FLAIR_20 <- NO_VS_FLAIR_20[x]
S_NO_VS_T2_20 <- NO_VS_T2_20[x]
S_NO_VS_T1_20 <- NO_VS_T1_20[x]
S_FLAIR_FIRST_3 <- FLAIR_FIRST_3[x]
S_FLAIR_SECOND_3 <- FLAIR_SECOND_3[x]
S_FLAIR_THIRD_3 <- FLAIR_THIRD_3[x]
S_T2_FIRST_3 <- T2_FIRST_3[x]
S_T2_SECOND_3 <- T2_SECOND_3[x]
S_T2_THIRD_3 <- T2_THIRD_3[x]
S_T1_FIRST_3 <- T1_FIRST_3[x]
S_T1_SECOND_3 <- T1_SECOND_3[x]
S_T1_THIRD_3 <- T1_THIRD_3[x]
S_FLAIR_FIRST_5 <- FLAIR_FIRST_5[x]
S_FLAIR_SECOND_5 <- FLAIR_SECOND_5[x]
S_FLAIR_THIRD_5 <- FLAIR_THIRD_5[x]
S_T2_FIRST_5 <- T2_FIRST_5[x]
S_T2_SECOND_5 <- T2_SECOND_5[x]
S_T2_THIRD_5 <- T2_THIRD_5[x]
S_T1_FIRST_5 <- T1_FIRST_5[x]
S_T1_SECOND_5 <- T1_SECOND_5[x]
S_T1_THIRD_5 <- T1_THIRD_5[x]
S_NEIGHBOR_3 <- NEIGHBOR_3[x]
S_NEIGHBOR_5 <- NEIGHBOR_5[x]
S_SUBJECT <- SUBJECT[x]

label <- c(S_GOLD_Lesions + 2 * S_NEIGHBOR_3 + 3 * S_NEIGHBOR_5)
subject_it <- as.numeric(levels(as.factor(SUBJECT)))


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##
## Plots for the Un-Normalized images, FLAIR, T2, and T1
##
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##########################################################################################################################################
##Histgram of Lesions and Non-Lesions for each image individually, seperated by subject
##########################################################################################################################################


##FLAIR Image

par(mfrow=c(2,1))
plot(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 0 & SUBJECT == as.numeric(subject_it[1])]),col = "gray70", main = "FLAIR, No Lesions", xlim = c(min(NO_VS_NO_NORM_FLAIR), max(NO_VS_NO_NORM_FLAIR)), ylim = c(0, .00025))
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 0 & SUBJECT == subject_it[2]]),col = "gray55", lty = 2) 
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 0 & SUBJECT == subject_it[3]]),col = "gray40", lty = 3) 
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 0 & SUBJECT == subject_it[4]]),col = "gray25", lty = 4) 
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 0 &  SUBJECT == subject_it[5]]),col = "gray10", lty = 5) 

plot(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 1 & SUBJECT == subject_it[1]]),  col = "blue", main = "FLAIR, Lesions", xlim = c(min(NO_VS_NO_NORM_FLAIR), max(NO_VS_NO_NORM_FLAIR)), ylim = c(0, .00025))
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 1 & SUBJECT == subject_it[2]]),  col = "blue1", lty = 2)
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 1 & SUBJECT == subject_it[3]]),  col = "blue2", lty = 3)
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 1 & SUBJECT == subject_it[4]]),  col = "blue3", lty = 4)
lines(density(NO_VS_NO_NORM_FLAIR[GOLD_Lesions== 1 & SUBJECT == subject_it[5]]),  col = "blue4", lty = 5)


##T2 Image

par(mfrow=c(2,1))
plot(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 0 & SUBJECT == as.numeric(subject_it[1])]),col = "gray70", main = "T2, No Lesions", xlim = c(min(NO_VS_NO_NORM_T2), max(NO_VS_NO_NORM_T2)), ylim = c(0, .005))
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 0 & SUBJECT == subject_it[2]]),col = "gray55", lty = 2) 
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 0 & SUBJECT == subject_it[3]]),col = "gray40", lty = 3) 
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 0 & SUBJECT == subject_it[4]]),col = "gray25", lty = 4) 
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 0 &  SUBJECT == subject_it[5]]),col = "gray10", lty = 5) 

plot(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 1 & SUBJECT == subject_it[1]]),  col = "blue", main = "T2, Lesions", xlim = c(min(NO_VS_NO_NORM_T2), max(NO_VS_NO_NORM_T2)), ylim = c(0, .005))
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 1 & SUBJECT == subject_it[2]]),  col = "blue1", lty = 2)
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 1 & SUBJECT == subject_it[3]]),  col = "blue2", lty = 3)
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 1 & SUBJECT == subject_it[4]]),  col = "blue3", lty = 4)
lines(density(NO_VS_NO_NORM_T2[GOLD_Lesions== 1 & SUBJECT == subject_it[5]]),  col = "blue4", lty = 5)


##T1 Image

par(mfrow=c(2,1))
plot(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 0 & SUBJECT == as.numeric(subject_it[1])]),col = "gray70", main = "T1, No Lesions", xlim = c(min(NO_VS_NO_NORM_T1), max(NO_VS_NO_NORM_T1)), ylim = c(0, .00001))
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 0 & SUBJECT == subject_it[2]]),col = "gray55", lty = 2) 
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 0 & SUBJECT == subject_it[3]]),col = "gray40", lty = 3) 
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 0 & SUBJECT == subject_it[4]]),col = "gray25", lty = 4) 
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 0 &  SUBJECT == subject_it[5]]),col = "gray10", lty = 5) 

plot(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 1 & SUBJECT == subject_it[1]]),  col = "blue", main = "T1, Lesions", xlim = c(min(NO_VS_NO_NORM_T1), max(NO_VS_NO_NORM_T1)), ylim = c(0, .00001))
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 1 & SUBJECT == subject_it[2]]),  col = "blue1", lty = 2)
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 1 & SUBJECT == subject_it[3]]),  col = "blue2", lty = 3)
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 1 & SUBJECT == subject_it[4]]),  col = "blue3", lty = 4)
lines(density(NO_VS_NO_NORM_T1[GOLD_Lesions== 1 & SUBJECT == subject_it[5]]),  col = "blue4", lty = 5)


##########################################################################################################################################
##2D Scatterplots, FLAIR, T2 and T1
##########################################################################################################################################

tiff(file = '2_way_not_normal_all.tif', antialias = "none")
scatterplotMatrix(~S_NO_VS_NO_NORM_FLAIR +S_NO_VS_NO_NORM_T1 +S_NO_VS_NO_NORM_T2|label, var.labels = c('FLAIR','T1', 'T2'), smooth = FALSE,     
                  col = c('gray', 'blue', 'green', 'purple'), diagonal = 'histogram', id.col = c('non lesion', 'lesion', 'neighbor 1', 'neighbor 2'), 
                  pch = c(20,20,20,20), cex = .25, legend.plot = FALSE, reg.line = FALSE, main = c("FLAIR, T1, and T2 Image Intesities"), )
dev.off()

##########################################################################################################################################
##3D Scatterplots, FLAIR, T2 and T1
##########################################################################################################################################

plot3d(S_NO_VS_NO_NORM_FLAIR[S_GOLD_Lesions == 1], S_NO_VS_NO_NORM_T2[S_GOLD_Lesions == 1], S_NO_VS_NO_NORM_T1[S_GOLD_Lesions == 1], col = 'blue', xlab = 'FLAIR', ylab = 'T2', zlab = 'T1', box = FALSE)
plot3d(S_NO_VS_NO_NORM_FLAIR[S_NEIGHBOR_3 == 0 & S_NEIGHBOR_5 == 0 & S_GOLD_Lesions == 0], S_NO_VS_NO_NORM_T2[S_NEIGHBOR_3 == 0 & 
                                                                                                                S_NEIGHBOR_5 == 0 & S_GOLD_Lesions == 0], S_NO_VS_NO_NORM_T1[S_NEIGHBOR_3 == 0 & S_NEIGHBOR_5 == 0 & S_GOLD_Lesions == 0], col = 'gray', add = T, alpha = .4)
plot3d(S_NO_VS_NO_NORM_FLAIR[S_NEIGHBOR_5 == 1], S_NO_VS_NO_NORM_T2[S_NEIGHBOR_5 == 1], S_NO_VS_NO_NORM_T1[S_NEIGHBOR_5 == 1], col = 'purple', add = TRUE)
plot3d(S_NO_VS_NO_NORM_FLAIR[S_NEIGHBOR_3 == 1], S_NO_VS_NO_NORM_T2[S_NEIGHBOR_3 == 1], S_NO_VS_NO_NORM_T1[S_NEIGHBOR_3 == 1], col = 'green', add = TRUE)


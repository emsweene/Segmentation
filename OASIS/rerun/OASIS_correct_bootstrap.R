#!/usr/bin/env Rscript 
########################
########################
# Predict from the updated models
########################
########################
##CLUSTER commands
##qrsh -l  mem_free=39G,h_vmem=40G
########################
readniigz <- function(file){
  name = gunzip(file, 
                temporary = TRUE, 
                overwrite = TRUE,
                remove = FALSE)
  image <- f.read.nifti.volume(name)
  file.remove(name)
  return(image)
}



candidate_voxels_one_modal <- function(
  kernel, cutoff, 
  brain_mask, modal) {
  
  dura_mask <-1*(brain_mask>0)[,,,1]	
  sigma.smooth<-diag(3,3)
  
  k.size<-kernel
  blur_modal<-GaussSmoothArray(modal ,
                               sigma=sigma.smooth,
                               ksize=k.size,
                               mask=dura_mask)
  
  brain_modal <- blur_modal[brain_mask == 1]
  
  cutpoint<-quantile(brain_modal, cutoff)
  modal_voxels<-1*(blur_modal>cutpoint)
  
  return(modal_voxels)
}

##load libraries## 
library(AnalyzeFMRI)
library(R.utils)
library(pbapply)
library(dplyr)
library(readr)

##SET SEED##
##set.seed(55)


##set data directorys##
n_boot = 1000
n_ms = 98
n_healthy = 33
n_boot_ms = 15
n_boot_healthy = 5

boots = lapply(seq(n_boot), function(seed) {
  set.seed(seed)

  subject_boot <- sample(n_ms, n_ms, 
    replace = TRUE) 

  ##set training and validation set##
  sets <- c(1:n_ms)
  sets_t <- sample(n_ms, n_boot_ms, replace = FALSE)
  training <- subject_boot[sets_t]
  valid <- (sets[-sets_t])
  valid <- (subject_boot[valid])

  h_subject_boot <- sample(n_healthy, n_healthy, 
    replace = TRUE) 

  ##set training and validation set##
  h_sets <- c(1:n_healthy)
  h_sets_t <- sample(1:n_healthy, n_boot_healthy, 
    replace = FALSE)
  h_training <- h_subject_boot[h_sets_t]
  h_valid <- (h_sets[-h_sets_t])
  h_valid <- (h_subject_boot[h_valid])
  L = list(
    training = training,
    valid = valid,
    h_training = h_training,
    h_valid = h_valid)
})

hopkins_mri = file.path(
  "/dcl01/smart/data/structural",
  "msmri/Hopkins_MRI/")
navid_data = paste0(
  hopkins_mri,
  "Navid_Data")

e_data = file.path(navid_data, 
                   "Elizabeth_Processed_MS_Data")
re_data = file.path(e_data, "OASIS_RePro")
##set working directory##

##set data directorys##
lesion_dir <- file.path(navid_data, 
                        'corrected_lesions/')
lesion_files <- dir(lesion_dir,full.names = TRUE)
data_dir <- file.path(e_data, "Register_to_Navid",
                      "NORMALIZED/")

blur_10_dir <- file.path(re_data, 'BLUR/divided/10/')
blur_20_dir <- file.path(re_data, 'BLUR/divided/20/')

mask_dir <- file.path(re_data, 'Masks/')


h_dir = file.path(hopkins_mri, 
                  "Healthy_Brain_Data_process")
flair_dir <- file.path(h_dir, "flair/")
flair_files <- dir(flair_dir,full.names = TRUE)

healthy_dir <- file.path(
  h_dir, 
  "N3_CoReg_Healthy_with_Spectre/")

dir.exists(lesion_dir)
dir.exists(data_dir)
dir.exists(blur_10_dir)
dir.exists(blur_20_dir)
dir.exists(mask_dir)


outfile = file.path(hopkins_mri, 
                    "full_prediction_data.rds")
  
  data_file = file.path(hopkins_mri, 
                        "correct_oasis_data.rds")
  data = read_rds(data_file)
  
  
  # data$oasis_p = predict(oasis_model,
  # 	newdata = data, type = "response")
  
  bad_data = data
  bad_data = data %>% 
    mutate(
      T1_orig = T1,
      T2_orig = T2,
      ms = grepl("MS", id),
      T1 = ifelse(ms, T2_orig, T1_orig),
      T2 = ifelse(ms, T1_orig, T2_orig)
    ) %>% 
    select(-T1_orig, -T2_orig, -ms)
  
  # bad_data_file = file.path(hopkins_mri, 
  # 		"original_oasis_data.rds")
  # bad_data = readRDS(bad_data_file)
  
  ###################
  ###################
  # Train the model
  ###################
  ###################
  
  
  ###################
  # Reads in the data for the MS subjects
  ###################
  
  all_df = vector(mode = "list", length = 98)
  i = 1
  
  for (i in 1:98) {
    
    print(i)
    
    lesions <- lesion_files[i]
    patient_id <- basename(lesions)
    patient_id <- sub("(.*)_mod.*", "\\1", 
                      patient_id)
    
    run_id = paste0("MS_", i)
    
    ind = which(data$id == run_id)
    df = data[ind, ]
    ### predict with correct model
    df$oasis_p = predict(
      oasis_model,
      newdata = df, type = "response")
    
    bad_df = bad_data[ind, ]
    ### predict with incorrect model
    df$bad_oasis_p = predict(
      original_oasis_model,
      newdata = bad_df, type = "response")
    
    
    ### predict with refined model
    df$oasis_p2 = predict(
      oasis_model_blur2,
      newdata = df, type = "response")  
    ### predict with incorrect refined model
    df$bad_oasis_p2 = predict(
      original_oasis_model_blur2,
      newdata = bad_df, type = "response")    
    
    rm(list = "bad_df")
    
    ################################
    ##load in files
    ################################
    
    gold_lesions <- readniigz(
      paste0(lesion_dir, patient_id, 
             '_mod.nii.gz'))
    dura_mask <- readniigz(
      paste0(mask_dir, patient_id, 
             'mask.nii.gz'))
    
    flair <- readniigz(
      paste0(data_dir, patient_id, 
             '_normalized_flairmni.nii.gz'))
    topvoxels <- candidate_voxels_one_modal(
      1, .85 , 
      dura_mask, flair)
    
    flair = drop(flair)
    dura_mask = drop(dura_mask)
    topvoxels = drop(topvoxels)
    gold_lesions = drop(gold_lesions)
    
    mask_index = dura_mask == 1
    mask_index = which(mask_index)  
    gold_lesions <- gold_lesions[mask_index]
    
    mat = data_frame(y = gold_lesions)
    
    index = which(topvoxels == 1)
    
    
    sigma.smooth <- diag(3,3)
    k.size <- 5
    
    cn = c("oasis_p", 
           "bad_oasis_p",
           "oasis_p2", 
           "bad_oasis_p2")
    icn = cn[1]
    
    xarr = array(0, dim = dim(flair)[1:3])
    
    for (icn in cn) {
      arr = xarr
      arr[index] = df[, icn, drop = TRUE]
      
      arr <- GaussSmoothArray(arr,
                              sigma = sigma.smooth,
                              ksize = k.size,
                              mask = dura_mask)  
      arr = arr[mask_index]
      stopifnot(!any(is.na(arr)))
      
      arr = data_frame(arr)
      colnames(arr) = icn
      mat = bind_cols(mat, arr)
    }
    
    all_df[[i]] = mat
    
  }
  
  
  
  healthy_df = vector(mode = "list", length = 33)
  
  for (i in 1:33) {
    
    lesions <- flair_files[i]
    patient_id <- basename(lesions)
    patient_id <- sub("(.*)flair.*", "\\1", 
                      patient_id)
    
    ###################
    ##load in files
    ###################
    run_id = paste0("Control_", i)
    
    ind = which(data$id == run_id)
    df = data[ind, ]
    ### predict with correct model
    df$oasis_p = predict(
      oasis_model,
      newdata = df, type = "response")
    
    bad_df = bad_data[ind, ]
    ### predict with incorrect model
    df$bad_oasis_p = predict(
      original_oasis_model,
      newdata = bad_df, type = "response")
    
    
    ### predict with refined model
    df$oasis_p2 = predict(
      oasis_model_blur2,
      newdata = df, type = "response")  
    ### predict with incorrect refined model
    df$bad_oasis_p2 = predict(
      original_oasis_model_blur2,
      newdata = bad_df, type = "response")    
    
    rm(list = "bad_df")
    
    ###################
    ##load in files
    ###################
    
    
    csf_fname = paste0(healthy_dir, patient_id, 
                       'csfmask.nii.gz')
    dura_mask <- readniigz(csf_fname)
    flair <- readniigz(paste0(healthy_dir, 
                              patient_id, 
                              'FLAIRnorm.nii.gz'))
    
    topvoxels <- candidate_voxels_one_modal(
      1, .85 , 
      dura_mask, flair)
    
    flair = drop(flair)
    topvoxels = drop(topvoxels)
    dura_mask = drop(dura_mask)
    
    index = which(topvoxels==1)
    
    xarr = array(0, dim = dim(flair)[1:3])
    
    mask_index = dura_mask == 1
    mask_index = which(mask_index)
    gold_lesions <- rep(0, 
      length = length(mask_index))
    
    mat = data_frame(y = gold_lesions)
    
    for (icn in cn) {
      arr = xarr
      arr[index] = df[, icn, drop = TRUE]
      
      arr <- GaussSmoothArray(arr,
                              sigma = sigma.smooth,
                              ksize = k.size,
                              mask = dura_mask)  
      arr = arr[mask_index]
      stopifnot(!any(is.na(arr)))
      
      arr = data_frame(arr)
      colnames(arr) = icn
      mat = bind_cols(mat, arr)
    }
    
    healthy_df[[i]] = mat
    print(i)
    rm(list = c("flair", "mat", "gold_lesions"))
    
  }
  
  rm(list = c("bad_data", "data")); gc()
  gc()
  names(all_df) = paste0("MS_", 
                         seq_along(all_df))
  names(healthy_df) = paste0("Control_", 
                             seq_along(healthy_df))
  
  all_data = bind_rows(all_df, 
                       healthy_df, .id = "id")
  print(dim(all_data))
  
  rm(list = c("all_df", "healthy_df")); gc(); gc()
  
  write_rds(all_data, 
            path = outfile)
} else {
  all_data = read_rds(outfile)
}
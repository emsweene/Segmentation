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
library(ROCR)
library(R.utils)
library(pbapply)
library(dplyr)
library(readr)
library(data.table)
library(tidyr)
##SET SEED##
##set.seed(55)


##set data directorys##


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

if (!file.exists(outfile)) {
  load(file.path(hopkins_mri, 
                 "original_oasis_model.rda"))
  
  load(file.path(hopkins_mri, 
                 "oasis_model.rda"))
  
  load(file.path(hopkins_mri, 
                 "oasis_model_blur2.rda"))
  
  load(file.path(hopkins_mri, 
                 "original_oasis_model_blur2.rda"))
  
  
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

xdata = all_data %>% 
  select(-id) %>% 
  gather(key = pred, value = x, -y)
rm(all_data); gc()

setDT(xdata)
setorder(xdata, pred, -x)
setDF(xdata)

xdata = as_data_frame(xdata)
xdata = split(xdata, f = xdata$pred)
xnames = sapply(xdata, function(x) {{
  unique(x$pred)
}})
names(xdata) = xnames

df = xdata[[1]]

preds = pblapply(xdata, function(df){ 


  df$x = round(as.numeric(df$x), 5)
  df$y = as.integer(df$y)
  predictions = df$x
  labels = df$y

  n.pos = sum(df$y)
  N = nrow(df)
  n.neg = N - n.pos
  df$xrow = seq.int(N)
  df$tp = cumsum(df$y)
  df$fp = df$xrow - df$tp
  df$xrow = NULL
  # notdups <- !rev(duplicated(rev(df$x)))
  notdups = c(abs(diff(df$x)) < 1e-5, TRUE)

  df = df[notdups, ]
  df = rename(df, cutoffs = x)
  df = df[, c("tp", "fp", "cutoffs")]
  df$fn = n.pos - df$tp
  df$tn = n.neg - df$fp

  df$n.pos.pred = df$tp + df$fp
  df$n.neg.pred = df$tn + df$fn

  df = as.list(df)
  df$n.pos = n.pos
  df$n.neg = n.neg
  df$predictions = predictions
  df$labels = labels

  make_pred = function(...) {
    new("prediction", ...)
  }
  df = lapply(df, list)
  pred = do.call(make_pred, df)
  pred
  # df
})

names(preds) = xnames

print_auc = function(pred, fpr.stop = 1) {
  
  xx = performance(pred, 
    measure = "auc", 
    fpr.stop = fpr.stop[1])

  xx = as.numeric(xx@y.values)/fpr.stop[1]
  print(xx)
  return(xx)
}

quick_auc = function(pred, 
  fpr.stop = c(0.005, 0.01, 1)
  ) {

  xvals = pred@fp[[1]] / pred@n.neg[[1]]
  yvals = pred@tp[[1]] / pred@n.pos[[1]]
  finite.bool <- is.finite(xvals) & is.finite(yvals)
  xvals = xvals[finite.bool]
  yvals = yvals[finite.bool]
  cum_add = function(y) {
    n = length(y)
    y[2:n] + y[1:(n-1)]
  }  
  aucs = lapply(fpr.stop, function(r) {
    ind = max(which(xvals <= r))
    if (ind == length(xvals)) {
      tpr.stop = yvals[ind]
    } else {
    tpr.stop <- approxfun(
      xvals[ind:(ind + 1)], 
      yvals[ind:(ind + 1)])(r)
    }
    xx <- c(xvals[1:ind], r)
    yy <- c(yvals[1:ind], tpr.stop)
    auc = sum(diff(xx) * cum_add(yy) * 0.5)
    auc = auc / r
    return(auc)
  })
  names(aucs) = fpr.stop
  print(aucs)
  return(aucs)
}

auc = pblapply(preds, quick_auc)

# names(preds)
# auc = sapply(preds, print_auc,
#   fpr.stop = 1)
# pauc_0.01 = sapply(preds, print_auc, 
#   fpr.stop = 0.01)
# pauc_0.005 = sapply(preds, print_auc, 
#   fpr.stop = 0.005)



# xdata = xdata %>% mutate_all(round, 5)

# ydata = sapply(seq(ncol(xdata)), function(z) {
#   data_frame()
# })


# outfile = file.path(hopkins_mri, 
# 	"vox_pred.rds")

vox_pred <-prediction(xdata, ydata)

# saveRDS(vox_pred, 
# 	file = outfile)

vox_perf<-performance(
  vox_pred, "tpr",  "fpr")





# y <- as.numeric(vox_perf@y.values[[1]])
# x <-  as.numeric(vox_perf@x.values[[1]])
# alpha <- as.numeric(
# 	vox_perf@alpha.values[[1]])

# Blur_x_values <- c()
# Blur_y_values <- c()
# Blur_alpha_values <- c()

# fpr_stop <- .01

# n <- length(x)


# n <- length(x)
# r <- n/1000

# for (i in 1:r){
# 	if (x[i*1000] < fpr_stop) 
# 	{ 
# 		Blur_x_values[i] <- x[i*1000]
# 		Blur_y_values[i] <- y[i*1000]
# 		Blur_alpha_values[i] <- alpha[i*1000]
# 	}
# 	else{
# 		stop
# 	}
# }


# movie_alpha <- min(Blur_alpha_values)
# sensitivity <- max(Blur_y_values)
# specificity <- max(Blur_x_values)

# print(movie_alpha)
# print(sensitivity)
# print(specificity)





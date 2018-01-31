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
library(dplyr)
library(readr)
library(ROCR)
library(pbapply)

##SET SEED##
##set.seed(55)

smoother = function(arr, dura_mask) {
  sigma.smooth <- diag(3,3)
  k.size <- 5
  arr <- GaussSmoothArray(
    arr,
    sigma = sigma.smooth,
    ksize = k.size,
    mask = dura_mask)  
}

make_mask = function(df, col = "dura_mask") {
  dimg = attributes(df)$dimg
  df = as.data.frame(df)
  arr = array(df[, col], dim = dimg)
  return(arr)
}

smooth_data = function(df, p, p2) {
  df$p[ df$topvoxels ] = p
  df$p2[ df$topvoxels ] = p2
  dura_mask = make_mask(df, "dura_mask")
  p = make_mask(df, "p")
  p2 = make_mask(df, "p2")
  sm_p = smoother(p, dura_mask)
  sm_p2 = smoother(p2, dura_mask)

  mat = data.frame(p = c(p),
    p2 = c(p2),
    sm_p = c(sm_p),
    sm_p2 = c(sm_p2)
    )
  mat$gold_lesions = df$gold_lesions
  mat = mat[ df$dura_mask, ]
  rownames(mat) = NULL
  return(mat)
}


##set data directorys##
set.seed(20180131)
n_boot = 1000
n_ms = 98
n_healthy = 33
n_boot_ms = 15
n_boot_healthy = 5

boots = lapply(seq(n_boot), function(x) {
  ids = paste0("MS_", seq(n_ms))
  ms = sample(ids)
  ids = paste0("Control_", seq(n_healthy))
  healthy = sample(ids)  
  df = data_frame(
    id = c(ms, healthy),
    train = c(
      rep(TRUE, n_boot_ms),
      rep(FALSE, n_ms - n_boot_ms),
      rep(TRUE, n_boot_healthy),
      rep(FALSE, n_healthy - n_boot_healthy)
      ),
    patient_group = sub("_.*", "", id),
    boot = x     
  )
  df = df %>% 
    arrange(desc(patient_group), 
      desc(train), id)
  return(df) 
})


hopkins_mri = file.path(
  "/dcl01/smart/data/structural",
  "msmri/Hopkins_MRI/")
outfile = file.path(hopkins_mri, 
                    "full_prediction_data.rds")
  
voxsel_file = file.path(
  hopkins_mri, 
  "full_voxsel_data.rds")
voxsel = read_rds(voxsel_file)

voxsel = lapply(voxsel, function(x) {
  x$p = x$p2 = 0
  x
})

data_file = file.path(hopkins_mri, 
                      "correct_oasis_data.rds")
data = read_rds(data_file)
xcn = colnames(data)

iboot = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iboot)) {
  iboot = 1
}

bad = FALSE
if (bad) {
# mess up the data like in the original paper
  data = data %>% 
  mutate(
    T1_orig = T1,
    T2_orig = T2,
    ms = grepl("MS", id),
    T1 = ifelse(ms, T2_orig, T1_orig),
    T2 = ifelse(ms, T1_orig, T2_orig)
    ) %>% 
  select(-T1_orig, -T2_orig, -ms)
}


boot_df = boots[[iboot]]
data = data[, xcn]

# getting the training data
data = left_join(data, boot_df, by = "id")

train_df = data %>% 
  filter(train)
valid_df = data %>% 
  filter(!train)    

rm(data)

#########################
# Should we rerun the bad data?
#########################  


###################################
# Correct Full Model
###################################
model<- glm(GOLD_Lesions ~ 
  FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  PD_10 * PD + PD_20 * PD + 
  T2_10 * T2 + T2_20 * T2 + 
  T1_10 * T1 + T1_20 * T1 , 
  data = train_df,
  family = binomial(),
  control= list(trace = TRUE)) 

# double smooth model
model2 <- glm(GOLD_Lesions ~ 
  FLAIR_2_10 * FLAIR + FLAIR_2_20 * FLAIR + 
  PD_2_10 * PD + PD_2_20 * PD + 
  T2_2_10 * T2 + T2_2_20 * T2 + 
  T1_2_10 * T1 + T1_2_20 * T1, 
  data = train_df, 
  family=binomial(),
  control= list(trace = TRUE)
  )


# train_df$p = predict(model, type = "response")
valid_df$p = predict(model, 
  newdata = valid_df,
  type = "response")

# train_df$p2 = predict(model2, type = "response")
valid_df$p2 = predict(model2, 
  newdata = valid_df,
  type = "response")

rm(train_df)
model = keep_mod(model)
model2 = keep_mod(model2)

valid_df = split(valid_df, valid_df$id)
ids = names(valid_df)


iid = 2
id = ids[iid]
res = pblapply(ids, 
  function(id) {
  smooth_data(
    df = voxsel[[id]],
    p = valid_df[[id]]$p,
    p2 = valid_df[[id]]$p2
    )
  })
rm(valid_df)

res = bind_rows(res)


p = prediction(res$p, res$gold_lesions)
p2 = prediction(res$p2, res$gold_lesions)
sm_p = prediction(res$sm_p, res$gold_lesions)
sm_p2 = prediction(res$sm_p2, res$gold_lesions)



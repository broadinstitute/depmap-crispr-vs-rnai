
library(purrr)
library(data.table)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
sparkle_params <- args[1]
feat_suffix <- args[2]
pred_suffix <- args[3]
tmp_file_dir <- args[4]
data_name <- args[5]

param_table <- read.csv(sparkle_params,stringsAsFactors = F)
param_table$featfile <- paste0(param_table$model,"_",param_table$start,"_",param_table$end,"_",feat_suffix)
param_table$predfile <- paste0(param_table$model,"_",param_table$start,"_",param_table$end,"_",pred_suffix)

#### HACK for compiling with missing task files
# existing_files <- list.files(tmp_file_dir)
# if(!all(param_table$featfile %in% existing_files)){
#   example <- param_table$featfile[!(param_table$featfile %in% existing_files)]
#   stop(paste0("Missing feature.csv files -- ",example[1]))
# }
# param_table %<>% subset(.,featfile %in% existing_files)
# param_table %<>% subset(.,predfile %in% existing_files)
# 
# write_csv(param_table,sparkle_params)


for (m in unique(param_table$model)){
  model_files <- subset(param_table,model == m)
  pred_files <- model_files$predfile
  feat_files <- model_files$featfile
  
  #If rownames are consistant across all files, use cbind instead of full join
  pred_list <- list()
  for (i in 1:length(pred_files)){
    tmp_df <- fread(file.path(tmp_file_dir,pred_files[i])) %>% as.data.frame(.)
    if (ncol(tmp_df) > 1){
      pred_list[[as.character(i)]] <- tmp_df
    }
  }
  
  t <- table(unlist(lapply(pred_list,nrow)))
  common_length <- as.numeric(names(t)[which.max(t)])
  bind_pred_list <- Filter(function(x){nrow(x) == common_length},pred_list)
  join_pred_list <- Filter(function(x){nrow(x) != common_length},pred_list)
  rm(pred_list)
  
  bind_pred_list <- lapply(bind_pred_list,function(x){x %<>% column_to_rownames(.,var="Row.name")})
  cl_order <- lapply(bind_pred_list,rownames)
  cl_order <- Reduce(intersect,cl_order)
  if (length(cl_order) == common_length){
    bind_pred_list <- lapply(bind_pred_list,function(x){x[cl_order,]})
    bind_pred_list <- bind_cols(bind_pred_list)
    rownames(bind_pred_list) <- cl_order
    bind_pred_list %<>% rownames_to_column(.,var="Row.name")
  } else {
    stop("Cannot bind columns of predictions across common length")
  }
  
  if (length(join_pred_list) > 1){
    join_pred_list <- join_pred_list %>% reduce(full_join, by="Row.name")
    model_preds <- full_join(bind_pred_list,join_pred_list,by="Row.name")
    write_csv(model_preds,paste0(data_name,"-",m,"-",pred_suffix))
  } else {
    write_csv(bind_pred_list,paste0(data_name,"-",m,"-",pred_suffix))
  }
  
  feat_list <- list()
  for (i in 1:length(feat_files)){
    tmp_df <- fread(file.path(tmp_file_dir,feat_files[i]))
    if (nrow(tmp_df) > 1){
      feat_list[[as.character(i)]] <- tmp_df
    }
  }
  model_feats <- bind_rows(feat_list)
  write_csv(model_feats,paste0(data_name,"-",m,"-",feat_suffix))
  
}

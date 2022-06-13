source("src/packages_paths.R")

#Accuracy cutoff based on CRISPR and RNAi null distributions 
accuracy_cutoff <- .5

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  select(.,gene,model,pearson,feature0)

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  select(.,gene,model,pearson,feature0)

related <- fread(file.path(data_raw,"gene-set-related-features.csv"))
paralogs <- subset(related,source == "paralog")

check_self <- function(target_CDS, feature_ensemble){
  
  feat_type <- strsplit(feature_ensemble,"_")[[1]]
  feat_type <- feat_type[length(feat_type)]
  
  result = paste0("SL_",feat_type)
  
  if (feat_type %in% c("RNAseq","CN","MutDam","MutMis","MutHot","RRBS")){
    
    feat_name = gsub(paste0("_",feat_type),"",feature_ensemble)
    if (feat_name == target_CDS){
      result = paste0("Self_",feat_type)
    } 
    
  } else if (feat_type %in% c("RPPA","proteomics")){
    
    feat_symbol = gsub(" .*","",feature_ensemble)
    target_symbol = gsub(" .*","",target_CDS)
    
    if (feat_symbol == target_symbol){
      result = paste0("Self_",feat_type)
    }
    
  } else if (feat_type %in% c("Fusion")) {
    
    target_symbol = gsub(" .*","",target_CDS)
    feat_name = gsub(paste0("_",feat_type),"",feature_ensemble)
    feat_symbols = strsplit(feat_name,"_")[[1]]
    
    if (target_symbol %in% feat_symbols){
      result = paste0("Self_",feat_type)
    } 
  }
  
  return(result)
}

check_paralog <- function(target_CDS, feature_ensemble, related_table){
  
  related_vec <- subset(related_table,target == target_CDS)$partner
  
  feat_type <- strsplit(feature_ensemble,"_")[[1]]
  feat_type <- feat_type[length(feat_type)]
  
  result = FALSE
  
  if (feat_type %in% c("RNAseq","CN","MutDam","MutMis","MutHot","RRBS")){
    
    feat_name = gsub(paste0("_",feat_type),"",feature_ensemble)
    if (feat_name %in% related_vec){
      result = T
    } 
  } else if (feat_type %in% c("RPPA","proteomics")){
    
    feat_symbol = gsub(" .*","",feature_ensemble)
    related_symbols = gsub(" .*","",related_vec)
    
    if (feat_symbol %in% related_symbols){
      result = T
    }
    
  }
  
  return(result)
  
}

crispr_res %<>% rename(.,CDS_ID=gene,ensemble_pearson=pearson) %>% 
  mutate(.,tech="CRISPR")
rnai_res %<>% rename(.,CDS_ID=gene,ensemble_pearson=pearson) %>% 
  mutate(.,tech="RNAi")

mark_paralogs <- function(results_summary){
  results_summary$self_feat <- "SL"
  results_summary$paralog_feat <- FALSE
  
  for (i in 1:nrow(results_summary)){
    self_res = check_self(target_CDS=results_summary$CDS_ID[i], feature_ensemble=results_summary$feature0[i])
    if (grepl("^SL_",self_res)){
      results_summary$paralog_feat[i] <- check_paralog(target_CDS=results_summary$CDS_ID[i], feature_ensemble=results_summary$feature0[i], paralogs)
    }
    results_summary$self_feat[i] <- self_res
  }
  
  return(results_summary)
  
}

crispr_res <- mark_paralogs(crispr_res)
crispr_res %<>% subset(.,paralog_feat)
crispr_res <- crispr_res[order(crispr_res$ensemble_pearson,decreasing=T),]
crispr_res %<>% subset(.,!duplicated(crispr_res$CDS_ID))
crispr_paralog <- dplyr::select(crispr_res,gene=CDS_ID,model) %>%
  mutate(.,class="paralog",type="CRISPR")

rnai_res <- mark_paralogs(rnai_res)
rnai_res %<>% subset(.,paralog_feat)
rnai_res <- rnai_res[order(rnai_res$ensemble_pearson,decreasing=T),]
rnai_res %<>% subset(.,!duplicated(rnai_res$CDS_ID))
rnai_paralog <- dplyr::select(rnai_res,gene=CDS_ID,model) %>%
  mutate(.,class="paralog",type="RNAi")

paralog_df <- rbind(crispr_paralog,rnai_paralog)

write_csv(paralog_df,file.path(data_processed,"predictive_marker_class_paralog.csv"))


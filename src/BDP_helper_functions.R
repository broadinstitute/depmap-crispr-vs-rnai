split_features <- function(result_summary){
  keep_cols <- c("gene",paste0("feature",0:9))
  result_summary <- as.data.frame(result_summary)
  result_summary <- result_summary[,colnames(result_summary) %in% keep_cols]
  return(result_summary)
}

split_importance <- function(result_summary){
  keep_cols <- c("gene",paste0("feature",0:9,"_importance"))
  result_summary <- as.data.frame(result_summary)
  result_summary <- result_summary[,colnames(result_summary) %in% keep_cols]
  return(result_summary)
}

split_imp_features <- function(result_summary,imp_cutoff){
  feats <- split_features(result_summary)
  imps <- split_importance(result_summary)
  
  stopifnot(all(feats$gene == imps$gene))
  
  genes <- feats$gene
  feats %<>% dplyr::select(.,-gene)
  imps %<>% dplyr::select(.,-gene) %>% as.matrix(.)
  feats[imps < imp_cutoff] <- NA
  
  feats %<>% add_column(gene=genes,.before=1)
  return(feats)
}

filter_feature_type <- function(feature_df,feat_type){
  genes <- feature_df$gene
  feature_df %<>% dplyr::select(.,-gene)
  for (j in 1:ncol(feature_df)){
    feature_df[!grepl(paste0("_",feat_type),feature_df[,j]),j] <- NA
  }
  
  feature_df %<>% add_column(.,gene=genes,.before=1)
  return(feature_df)
}


check_self_predictor <- function(feature_df,top_n=10,reduce=T){
  feature_df$self <- F
  predictors <- dplyr::select(feature_df,-gene,-self)
  feature_df %<>% dplyr::select(.,gene,self)
  
  target_lab <- gsub(" .*","",feature_df$gene)
  for (i in 1:nrow(feature_df)){
    feature_df$self[i] <- any(grepl(target_lab[i],predictors[i,1:top_n]))
  }
  
  if(reduce){
    hits <- unique(subset(feature_df,self)$gene)
    return(hits)
  }
  
  return(feature_df)
}

check_list_predictor <- function(feature_df,cds_id){
  feature_df$in_group <- F
  
  predictors <- dplyr::select(feature_df,-gene,-in_group)
  feature_df %<>% dplyr::select(.,gene,in_group)
  
  feat_symbols <- gsub(" .*","",cds_id)
  target_lab <- gsub(" .*","",feature_df$gene)
  for (i in 1:nrow(feature_df)){
    tmp_feat_symbols <- feat_symbols[!(feat_symbols %in% target_lab[i])]
    tmp_predictors <- gsub(" .*","",predictors[i,]) 
    feature_df$in_group[i] <- any(tmp_feat_symbols %in% tmp_predictors)
  }
  return(feature_df)
}

get_cor_dir <- function(dep_mat,feat_mat,target_ID,feat_ID,feat_type){
  
  cls = intersect(rownames(dep_mat),rownames(feat_mat))
  
  feat_name = gsub(paste0("_",feat_type),"",feat_ID)
  
  return(cor(dep_mat[cls,target_ID],feat_mat[cls,feat_name],use="pairwise.complete"))
  
}

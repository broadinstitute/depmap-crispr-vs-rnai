source("src/packages_paths.R")

#Accuracy cutoff based on CRISPR and RNAi null distributions 
accuracy_cutoff <- .5

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff)

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff)

self_expr_filter <- function(result_summary){
  result_summary %<>% subset(.,grepl("_RNAseq",result_summary$feature0))
  result_summary %<>% subset(.,gsub(" .*","",result_summary$feature0) == gsub(" .*","",result_summary$gene))
  return(result_summary)
}


crispr_res <- self_expr_filter(crispr_res)
rnai_res <- self_expr_filter(rnai_res)

feat_list <- list("RNAseq"="depmap-omics-expression-rnaseq-tpm-18Q4.csv")
feat_list <- lapply(feat_list,function(x){fread(file.path(data_raw,x)) %>% 
    column_to_rownames(.,var="Row.name") %>% 
    as.matrix(.)})

dep_list <- list("CRISPR"="gene-effect-scaled-crispr-matched.csv",
                 "RNAi"="gene-effect-scaled-rnai-matched.csv")
dep_list <- lapply(dep_list,function(x){fread(file.path(data_raw,x)) %>% 
    column_to_rownames(.,var="Row.name") %>% 
    as.matrix(.)})


stopifnot(all(rownames(dep_list[["CRISPR"]]) == rownames(dep_list[["RNAi"]])))
feat_list <- lapply(feat_list,function(x){x[rownames(x) %in% rownames(dep_list[["CRISPR"]]),]})

expr_feat_cor <- function(result_summary,pert_type){
  
  result_summary$topfeat_cor <- NA
  for (i in 1:nrow(result_summary)){
    result_summary$topfeat_cor[i] <- get_cor_dir(dep_mat=dep_list[[pert_type]],
                                                 feat_mat=feat_list[["RNAseq"]],
                                                 target_ID=result_summary$gene[i],
                                                 feat_ID=result_summary$feature0[i],
                                                 feat_type="RNAseq")
    
    
  }
  
  result_summary %<>% subset(.,topfeat_cor < 0)
  result_summary <- result_summary[order(result_summary$pearson,decreasing=T),]
  result_summary %<>% subset(.,!duplicated(result_summary$gene))
  return(result_summary)
}

crispr_expAddict <- expr_feat_cor(crispr_res,"CRISPR")
rnai_expAddict <- expr_feat_cor(rnai_res,"RNAi")

crispr_expAddict %<>% dplyr::select(.,gene,model) %>%
  mutate(.,class="expression addiction",type="CRISPR")

rnai_expAddict %<>% dplyr::select(.,gene,model) %>%
  mutate(.,class="expression addiction",type="RNAi")

expAddict_df <- rbind(crispr_expAddict,rnai_expAddict)

write_csv(expAddict_df,file.path(data_processed,"predictive_marker_class_expression_addiction.csv"))


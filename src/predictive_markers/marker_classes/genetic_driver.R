source("src/packages_paths.R")

#Accuracy cutoff based on CRISPR and RNAi null distributions 
accuracy_cutoff <- .5

#Filter for genes with accuracy > .5 and skewness < 0
t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2 %<>% add_column(.,CDS_ID=entrez_to_cds(t2$entrez_id,hgnc))

rnai_left_skew <- subset(t2,RNAi_Skewness < 0)
crispr_left_skew <- subset(t2,CRISPR_Skewness < 0)

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,gene %in% crispr_left_skew$CDS_ID)

crispr_good <- subset(crispr_res,pearson >= accuracy_cutoff)
crispr_res %<>% subset(.,gene %in% crispr_good$gene)

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,gene %in% rnai_left_skew$CDS_ID)

rnai_good <- subset(rnai_res,pearson >= accuracy_cutoff)
rnai_res %<>% subset(.,gene %in% rnai_good$gene)

#### Filter features for self hotspot, missense, fusion, or CN with importance > .05

get_top_driver_feature <- function(result_summary){
  feature_df <- split_imp_features(result_summary,imp_cutoff=.05)
  feature_df %<>% dplyr::select(.,-gene)
  
  target_lab <- gsub(" .*","",result_summary$gene)
  for (i in 1:nrow(feature_df)){
    feature_df[i,!grepl(paste0("^",target_lab[i]," "),feature_df[i,])] <- NA
  }
  
  valid_feats <- c("MutHot","MutMis","Fusion","CN")
  genes <- result_summary$gene
  for (j in 1:ncol(feature_df)){
    feat_vec <- gsub(".+_","",feature_df[,j])
    feature_df[!(feat_vec %in% valid_feats),j] <- NA
  }
  
  feature_df$top_feat <- NA
  for (i in 1:nrow(feature_df)){
    top_feat <- feature_df[i,]
    top_feat <- top_feat[!is.na(top_feat)]
    if (length(top_feat) > 0){
      feature_df$top_feat[i] <- top_feat[1]
    }
  }
  
  feature_df <- add_column(feature_df,gene=genes,.before=1)
  feature_df %<>% dplyr::select(.,gene,top_feat)
  
  result_summary <- cbind(feature_df,result_summary)
  result_summary %<>% subset(.,!is.na(top_feat))
  result_summary %<>% dplyr::select(.,-gene.1)
  return(result_summary)
}

crispr_res <- get_top_driver_feature(crispr_res)
rnai_res <- get_top_driver_feature(rnai_res)

#### Check that top features that use CN are negative correlations

feat_list <- list("RNAseq"="depmap-omics-expression-rnaseq-tpm-18Q4.csv",
                  "CN"="depmap-omics-cn-gene-internal-18q4.csv")
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

cn_cor <- function(result_summary,pert_type){
  
  result_summary$cn_cor <- NA
  for (i in 1:nrow(result_summary)){
    result_summary$cn_cor[i] <- get_cor_dir(dep_mat=dep_list[[pert_type]],
                                            feat_mat=feat_list[["CN"]],
                                            target_ID=result_summary$gene[i],
                                            feat_ID=result_summary$gene[i],
                                            feat_type="CN")
    
    
  }
  
  return(result_summary)
}

crispr_res <- cn_cor(crispr_res,"CRISPR")
crispr_res$mutation <- !grepl("_CN",crispr_res$top_feat)
crispr_res %<>% subset(.,mutation | (cn_cor < 0))
crispr_res <- crispr_res[order(crispr_res$pearson,decreasing=T),]
crispr_res %<>% subset(.,!duplicated(crispr_res$gene))
crispr_driver <- dplyr::select(crispr_res,gene,model) %>%
  mutate(.,class="genetic driver",type="CRISPR")

rnai_res <- cn_cor(rnai_res,"RNAi")
rnai_res$mutation <- !grepl("_CN",rnai_res$top_feat)
rnai_res %<>% subset(.,mutation | (cn_cor < 0))
rnai_res <- rnai_res[order(rnai_res$pearson,decreasing=T),]
rnai_res %<>% subset(.,!duplicated(rnai_res$gene))
rnai_driver <- dplyr::select(rnai_res,gene,model) %>%
  mutate(.,class="genetic driver",type="RNAi")

driver_df <- rbind(crispr_driver,rnai_driver)

write_csv(driver_df,file.path(data_processed,"predictive_marker_class_genetic_driver.csv"))



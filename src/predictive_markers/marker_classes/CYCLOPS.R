source("src/packages_paths.R")

#Accuracy cutoff based on CRISPR and RNAi null distributions 
accuracy_cutoff <- .5

#### Filter models for left skewed genes and predictive accuracy
t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2 %<>% add_column(.,CDS_ID=entrez_to_cds(t2$entrez_id,hgnc))

rnai_left_skew <- subset(t2,RNAi_Skewness < 0)
crispr_left_skew <- subset(t2,CRISPR_Skewness < 0)

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  subset(.,gene %in% crispr_left_skew$CDS_ID)

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  subset(.,gene %in% rnai_left_skew$CDS_ID)

#### Filter for models with CN or expression as top feature
top_cn_filter <- function(result_summary){
  result_summary %<>% subset(.,grepl("_CN",result_summary$feature0) | grepl("_RNAseq",result_summary$feature0))
  return(result_summary)
}

crispr_res <- top_cn_filter(crispr_res)
rnai_res <- top_cn_filter(rnai_res)

#### Filter for models where top feature is on the same chromosome arm as target gene
gene_loc <- fread(file.path("data/processed","hgnc_gene_arm_location.csv"))

arm_filter <- function(result_summary,loc_map){
  
  result_summary$feature0_symbol <- gsub("_CN","",result_summary$feature0)
  result_summary$feature0_symbol <- gsub("_RNAseq","",result_summary$feature0_symbol)
  
  colnames(loc_map) <- c("gene","gene_arm")
  result_summary %<>% left_join(.,loc_map,by="gene")
  colnames(loc_map) <- c("feature0_symbol","feature_arm")
  result_summary %<>% left_join(.,loc_map,by="feature0_symbol")
  result_summary$same_arm <- result_summary$gene_arm == result_summary$feature_arm
  
  result_summary %<>% subset(.,same_arm)
  return(result_summary)
  
}

crispr_res <- arm_filter(crispr_res,gene_loc)
rnai_res <- arm_filter(rnai_res,gene_loc)

#### Uses self CN or RNAseq feature with importance > .05

filter_cyclops_features <- function(result_summary){
  feature_df <- split_imp_features(result_summary,imp_cutoff=.05)
  
  genes <- feature_df$gene
  feature_df <- dplyr::select(feature_df,-gene)
  for (j in 1:ncol(feature_df)){
    feature_df[!(grepl("_CN",feature_df[,j]) | grepl("_RNAseq",feature_df[,j])),j] <- NA
  }
  feature_df <- add_column(feature_df,gene=genes,.before=1)
  
  self_res <- check_self_predictor(feature_df,top_n=10,reduce=F)
  result_summary %<>% subset(.,self_res$self)
  return(result_summary)
}

crispr_res <- filter_cyclops_features(crispr_res)
rnai_res <- filter_cyclops_features(rnai_res)

#### Correlation of top feature is positive

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

top_feat_cor <- function(result_summary,pert_type){
  
  result_summary$topfeat_type <- gsub(".+_","",result_summary$feature0)
  result_summary$topfeat_id <- gsub("_.*","",result_summary$feature0)
  result_summary$topfeat_cor <- NA
  for (i in 1:nrow(result_summary)){
    result_summary$topfeat_cor[i] <- get_cor_dir(dep_mat=dep_list[[pert_type]],
                                                 feat_mat=feat_list[[result_summary$topfeat_type[i]]],
                                                 target_ID=result_summary$gene[i],
                                                 feat_ID=result_summary$topfeat_id[i],
                                                 feat_type=result_summary$topfeat_type[i])
    
    
  }
  
  result_summary %<>% subset(.,topfeat_cor > 0)
  return(result_summary)
}

crispr_res <- top_feat_cor(crispr_res,"CRISPR")
rnai_res <- top_feat_cor(rnai_res,"RNAi")

#### Correlation to self CN is positive

cn_cor <- function(result_summary,pert_type){
  
  result_summary$cn_cor <- NA
  for (i in 1:nrow(result_summary)){
    result_summary$cn_cor[i] <- get_cor_dir(dep_mat=dep_list[[pert_type]],
                                            feat_mat=feat_list[["CN"]],
                                            target_ID=result_summary$gene[i],
                                            feat_ID=result_summary$gene[i],
                                            feat_type="CN")
    
    
  }
  
  result_summary %<>% subset(.,cn_cor > 0)
  return(result_summary)
}

crispr_res <- cn_cor(crispr_res,"CRISPR")
crispr_res <- crispr_res[order(crispr_res$pearson,decreasing=T),]
crispr_res %<>% subset(.,!duplicated(crispr_res$gene))
crispr_cyclops <- dplyr::select(crispr_res,gene,model) %>%
  mutate(.,class="CYCLOPS",type="CRISPR")

rnai_res <- cn_cor(rnai_res,"RNAi")
rnai_res <- rnai_res[order(rnai_res$pearson,decreasing=T),]
rnai_res %<>% subset(.,!duplicated(rnai_res$gene))
rnai_cyclops <- dplyr::select(rnai_res,gene,model) %>%
  mutate(.,class="CYCLOPS",type="RNAi")

cyclops_df <- rbind(crispr_cyclops,rnai_cyclops)

write_csv(cyclops_df,file.path(data_processed,"predictive_marker_class_CYCLOPS.csv"))



source("src/packages_paths.R")

accuracy_cutoff <- 0

#### Filter models for left skewed genes and predictive accuracy

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% add_column(.,CDS_ID=paste0(t2$symbol," (",t2$entrez_id,")"),.before=1)

rnai_left_skew <- subset(t2,RNAi_Skewness < 0)

rnai_res <- fread(file.path("data/processed","ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,best) %>% 
  subset(.,gene %in% subset(t2,CRISPR_PD)$CDS_ID)

#### Annotate models with CN or expression as top feature
top_cn_filter <- function(result_summary){
  result_summary$CN_or_Exp <- grepl("_CN",result_summary$feature0) | grepl("_RNAseq",result_summary$feature0)
  return(result_summary)
}

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
  result_summary$same_arm[is.na(result_summary$same_arm)] <- F
  
  return(result_summary)
  
}

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
  self_genes <- subset(self_res,self)$gene
  result_summary$self <- result_summary$gene %in% self_genes
  # result_summary %<>% subset(.,self_res$self)
  return(result_summary)
}

rnai_res <- filter_cyclops_features(rnai_res)

#### Correlation of top feature is positive

feat_list <- list("RNAseq"=fread(file.path("data/raw","depmap-omics-expression-rnaseq-tpm-18Q4.csv")),
                  "CN"=fread(file.path("data/raw","depmap-omics-cn-gene-internal-18q4.csv")))
feat_list <- lapply(feat_list,function(x){x %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

dep_list <- list("CRISPR"=fread(file.path("data/raw","gene-effect-scaled-crispr-matched.csv")),
                 "RNAi"=fread(file.path("data/raw","gene-effect-scaled-rnai-matched.csv")))
dep_list <- lapply(dep_list,function(x){x %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

stopifnot(all(rownames(dep_list[["CRISPR"]]) == rownames(dep_list[["RNAi"]])))
feat_list <- lapply(feat_list,function(x){x[rownames(x) %in% rownames(dep_list[["CRISPR"]]),]})

top_feat_cor <- function(result_summary,pert_type){
  
  result_summary$topfeat_type <- gsub(".+_","",result_summary$feature0)
  result_summary$topfeat_id <- gsub("_.*","",result_summary$feature0)
  result_summary$topfeat_cor <- NA
  for (i in 1:nrow(result_summary)){
    if (result_summary$CN_or_Exp[i]) {
      result_summary$topfeat_cor[i] <- get_cor_dir(dep_mat=dep_list[[pert_type]],
                                                   feat_mat=feat_list[[result_summary$topfeat_type[i]]],
                                                   target_ID=result_summary$gene[i],
                                                   feat_ID=result_summary$topfeat_id[i],
                                                   feat_type=result_summary$topfeat_type[i])
    }
    
    
    
  }
  
  return(result_summary)
}

rnai_res <- top_feat_cor(rnai_res,"RNAi")

#### Correlation to self CN is positive

get_cor_dir <- function(dep_mat,feat_mat,target_ID,feat_ID,feat_type){
  
  cls = intersect(rownames(dep_mat),rownames(feat_mat))
  
  feat_name = gsub(paste0("_",feat_type),"",feat_ID)
  
  return(cor(dep_mat[cls,target_ID],feat_mat[cls,feat_name],use="pairwise.complete"))
  
}

cn_cor <- function(result_summary,pert_type){
  
  result_summary$cn_cor <- NA
  for (i in 1:nrow(result_summary)){
    if (result_summary$gene[i] %in% colnames(feat_list[["CN"]])){
      result_summary$cn_cor[i] <- get_cor_dir(dep_mat=dep_list[[pert_type]],
                                              feat_mat=feat_list[["CN"]],
                                              target_ID=result_summary$gene[i],
                                              feat_ID=result_summary$gene[i],
                                              feat_type="CN")
    }
    
  }
  
  return(result_summary)
}

rnai_res <- cn_cor(rnai_res,"RNAi")
rnai_res$left_skew <- rnai_res$gene %in% rnai_left_skew$CDS_ID

self_exp_pos <- subset(rnai_res,self & (topfeat_type == "RNAseq") & (topfeat_cor > 0) & left_skew)$gene
rnai_res$self_exp_pos <- rnai_res$gene %in% self_exp_pos

local_cn_pos <- subset(rnai_res, same_arm & (topfeat_type == "CN") & (topfeat_cor > 0))$gene
rnai_res$local_cn_pos <- rnai_res$gene %in% local_cn_pos

write_csv(rnai_res,file.path("data/processed","CYCLOPS_annotations.csv"))


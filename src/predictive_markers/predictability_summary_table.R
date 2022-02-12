

source("src/packages_paths.R")

crispr_pred <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,!is.na(gene))
crispr_pred$ID <- paste0(crispr_pred$gene,"::",crispr_pred$model)

crispr_pred %<>% dplyr::select(.,ID,CRISPR_accuracy=pearson,CRISPR_top_feature=feature0,CRISPR_top_feature_importance=feature0_importance)

rnai_pred <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,!is.na(gene))
rnai_pred$ID <- paste0(rnai_pred$gene,"::",rnai_pred$model)
  
rnai_pred %<>% dplyr::select(.,ID,RNAi_accuracy=pearson,RNAi_top_feature=feature0,RNAi_top_feature_importance=feature0_importance)

summary_pred <- inner_join(crispr_pred,rnai_pred,by="ID")
summary_pred %<>% add_column(., Gene=gsub("::.*","",summary_pred$ID),.before=1)
summary_pred %<>% add_column(., Model=gsub(".+::","",summary_pred$ID),.after=1)
summary_pred %<>% dplyr::select(.,-ID)

summary_pred %<>% add_column(.,entrez_id=extract_entrez(summary_pred$Gene),.before=1)
hgnc <- fread(file.path("data/raw","hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
summary_pred %<>% add_column(.,symbol=entrez_to_symbol(summary_pred$entrez_id,hgnc),.before=1)

summary_pred %<>% rename(.,target_id=Gene)
write_csv(summary_pred,file.path("tables","Supplemental-Table-3.csv"))

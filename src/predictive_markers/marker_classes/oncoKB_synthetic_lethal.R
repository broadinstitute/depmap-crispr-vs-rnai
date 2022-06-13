source("src/packages_paths.R")

#Accuracy cutoff based on CRISPR and RNAi null distributions 
accuracy_cutoff <- .5

onco <- fread(file.path(data_raw,"control-oncoKB.csv")) %>% 
  subset(.,source=="oncogene")
tsg <- fread(file.path(data_raw,"control-oncoKB.csv")) %>% 
  subset(.,source=="TSG")

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2 %<>% add_column(.,CDS_ID=entrez_to_cds(t2$entrez_id,hgnc))

rnai_left_skew <- subset(t2,RNAi_Skewness < 0)
crispr_left_skew <- subset(t2,CRISPR_Skewness < 0)


##### Oncogene synthetic lethal

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  subset(.,gene %in% crispr_left_skew$CDS_ID)

crispr_oncoSL <- check_list_predictor(split_imp_features(crispr_res,.05),onco$CDS_ID)
crispr_res <- crispr_res[crispr_oncoSL$in_group]
crispr_res <- crispr_res[order(crispr_res$pearson,decreasing=T),]
crispr_res %<>% subset(.,!duplicated(crispr_res$gene))

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  subset(.,gene %in% rnai_left_skew$CDS_ID)

rnai_oncoSL <- check_list_predictor(split_imp_features(rnai_res,.05),onco$CDS_ID)
rnai_res <- rnai_res[rnai_oncoSL$in_group] 
rnai_res <- rnai_res[order(rnai_res$pearson,decreasing=T),]
rnai_res %<>% subset(.,!duplicated(rnai_res$gene))

crispr_res %<>% dplyr::select(.,gene,model) %>%
  mutate(.,class="synthetic lethal (oncogene)",type="CRISPR")

rnai_res %<>% dplyr::select(.,gene,model) %>%
  mutate(.,class="synthetic lethal (oncogene)",type="RNAi")

oncoSL_df <- rbind(crispr_res,rnai_res)


##### TSG synthetic lethal

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  subset(.,gene %in% crispr_left_skew$CDS_ID)

crispr_tsgSL <- check_list_predictor(split_imp_features(crispr_res,.05),tsg$CDS_ID)
crispr_res <- crispr_res[crispr_tsgSL$in_group]
crispr_res <- crispr_res[order(crispr_res$pearson,decreasing=T),]
crispr_res %<>% subset(.,!duplicated(crispr_res$gene))

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,pearson >= accuracy_cutoff) %>%
  subset(.,gene %in% rnai_left_skew$CDS_ID)

rnai_tsgSL <- check_list_predictor(split_imp_features(rnai_res,.05),tsg$CDS_ID)
rnai_res <- rnai_res[rnai_tsgSL$in_group] 
rnai_res <- rnai_res[order(rnai_res$pearson,decreasing=T),]
rnai_res %<>% subset(.,!duplicated(rnai_res$gene))

crispr_res %<>% dplyr::select(.,gene,model) %>%
  mutate(.,class="synthetic lethal (TSG)",type="CRISPR")

rnai_res %<>% dplyr::select(.,gene,model) %>%
  mutate(.,class="synthetic lethal (TSG)",type="RNAi")

tsgSL_df <- rbind(crispr_res,rnai_res)

oncoTSG_SL <- rbind(oncoSL_df,tsgSL_df)

write_csv(oncoTSG_SL,file.path(data_processed,"predictive_marker_class_oncoKB_synthetic_lethal.csv"))


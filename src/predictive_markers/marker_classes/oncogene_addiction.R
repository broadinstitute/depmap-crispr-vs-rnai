source("src/packages_paths.R")

#Accuracy cutoff based on CRISPR and RNAi null distributions 
accuracy_cutoff <- .5

onco <- fread(file.path(data_raw,"control-oncoKB.csv")) %>% 
  subset(.,source=="oncogene")

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2 %<>% add_column(.,CDS_ID=entrez_to_cds(t2$entrez_id,hgnc))

rnai_left_skew <- subset(t2,RNAi_Skewness < 0)
crispr_left_skew <- subset(t2,CRISPR_Skewness < 0)

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,gene %in% onco$CDS_ID) %>%
  subset(.,gene %in% crispr_left_skew$CDS_ID)

crispr_good <- subset(crispr_res,pearson >= accuracy_cutoff)
crispr_res %<>% subset(.,gene %in% crispr_good$gene)

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,gene %in% onco$CDS_ID) %>%
  subset(.,gene %in% rnai_left_skew$CDS_ID)

rnai_good <- subset(rnai_res,pearson >= accuracy_cutoff)
rnai_res %<>% subset(.,gene %in% rnai_good$gene)

crispr_onco <- check_self_predictor(split_features(crispr_res),top_n=10)
rnai_onco <- check_self_predictor(split_features(rnai_res),top_n=10)

oncoAddict_df <- data.frame(gene=crispr_onco,class="oncogene addiction",type="CRISPR",stringsAsFactors = F)
oncoAddict_df <- rbind(oncoAddict_df,
                       data.frame(gene=rnai_onco,class="oncogene addiction",type="RNAi",stringsAsFactors = F))
write_csv(oncoAddict_df,file.path(data_processed,"predictive_marker_class_oncogene_addiction.csv"))


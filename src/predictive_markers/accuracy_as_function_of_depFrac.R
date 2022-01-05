
source("src/packages_paths.R")

crispr_res <- fread(file.path(data_processed,"gene-effect-ensemble-regression-crispr-matched.csv")) %>%
  subset(.,best) %>% 
  dplyr::select(.,gene,CRISPR_ensemble_pearson=pearson)

rnai_res <- fread(file.path(data_processed,"gene-effect-ensemble-regression-rnai-matched.csv")) %>%
  subset(.,best) %>% 
  dplyr::select(.,gene,RNAi_ensemble_pearson=pearson)

regression_res <- full_join(crispr_res,rnai_res,by="gene")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

master <- add_column(t2,CDS_ID=entrez_to_cds(t2$entrez_id,hgnc),.before=1)

master %<>% dplyr::select(.,CDS_ID,`CRISPR_numDeps(50-100)`,`CRISPR_numDeps(0-50)`,`RNAi_numDeps(50-100)`,`RNAi_numDeps(0-50)`)
master %<>% left_join(.,regression_res,by=c("CDS_ID"="gene"))

crispr_df <- dplyr::select(master,CDS_ID,CRISPR_ensemble_pearson,`CRISPR_numDeps(50-100)`,`CRISPR_numDeps(0-50)`) %>%
  rename(.,ensemble_pearson=CRISPR_ensemble_pearson,numDeps=`CRISPR_numDeps(50-100)`,nonDeps=`CRISPR_numDeps(0-50)`) %>%
  mutate(.,tech="CRISPR")
crispr_df$total_deps <- crispr_df$numDeps + crispr_df$nonDeps
crispr_df$depFrac <- crispr_df$numDeps / crispr_df$total_deps

rnai_df <- dplyr::select(master,CDS_ID,RNAi_ensemble_pearson,`CRISPR_numDeps(50-100)`,`CRISPR_numDeps(0-50)`) %>%
  rename(.,ensemble_pearson=RNAi_ensemble_pearson,numDeps=`CRISPR_numDeps(50-100)`,nonDeps=`CRISPR_numDeps(0-50)`) %>%
  mutate(.,tech="RNAi")
rnai_df$total_deps <- rnai_df$numDeps + rnai_df$nonDeps
rnai_df$depFrac <- rnai_df$numDeps / rnai_df$total_deps

plot_df <- bind_rows(crispr_df,rnai_df)

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

ggplot(plot_df,aes(x=depFrac,y=ensemble_pearson,color=tech)) +
  geom_smooth() +
  theme_bw(base_size=11) +
  ylab("Predictive Accuracy (r)") +
  xlab("CRISPR Dep. Cell Line Frac.") +
  scale_color_manual(values=tech_pal,name="") +
  theme(legend.position=c(.8,.2))
ggsave(file.path("figures","predictive_markers_smoothed_accuracy_funcOf_depFrac.pdf"),width=2.25,height=2)



dep_classes <- fread(file.path("tables","Supplemental-Table-2.csv")) %>% 
  dplyr::select(.,entrez_id,symbol,CRISPR_class,RNAi_class)
dep_classes$entrez_id %<>% as.character(.)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
dep_classes %<>% add_column(.,CDS_ID=entrez_to_cds(dep_classes$entrez_id,hgnc))

bdp_classes <- list(cyclops="CYCLOPS.csv",
                    expAddict="expression_addiction.csv",
                    driver="genetic_driver.csv",
                    onco="oncogene_addiction.csv",
                    oncoSL="oncoKB_synthetic_lethal.csv",
                    paralog="paralog.csv")
bdp_classes <- lapply(bdp_classes,function(x){fread(file.path(data_processed,paste0("predictive_marker_class_",x)))})


bdp_classes <- bind_rows(bdp_classes)

##### Check fraction of accurate models that are BDP
crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,best) %>% 
  mutate(.,perturbation="CRISPR")

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,best) %>% 
  mutate(.,perturbation="RNAi")

pred_res <- bind_rows(crispr_res,rnai_res)
pred_res$confounder <- grepl("_Confounders",pred_res$feature0)

pred_res %<>% subset(., (pearson > .5) & !confounder)

crispr_bdps <- unique(subset(bdp_classes,type == "CRISPR")$gene)
length(crispr_bdps)
nrow(subset(pred_res,(perturbation == "CRISPR") & (gene %in% crispr_bdps) ))  / nrow(subset(pred_res,(perturbation == "CRISPR")))

rnai_bdps <- unique(subset(bdp_classes,type == "RNAi")$gene)
length(rnai_bdps)
nrow(subset(pred_res,(perturbation == "RNAi") & (gene %in% rnai_bdps) ))  / nrow(subset(pred_res,(perturbation == "RNAi")))
##### end

#### Genome-wide plot

plot_supp <- bdp_classes

plot_supp$class <- gsub(" \\(oncogene\\)","",plot_supp$class)
plot_supp$class <- gsub(" \\(TSG\\)","",plot_supp$class)

plot_supp$class <- factor(plot_supp$class,levels=c("genetic driver","expression addiction","paralog","CYCLOPS","oncogene addiction","synthetic lethal"))

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

ggplot(plot_supp,aes(x=class,fill=type)) +
  geom_bar(stat="count",position="dodge",color="black",width=.7) +
  theme_bw(base_size=11) +
  # theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = c("none")) +
  # facet_grid(facet ~ ., scales = "free") +
  scale_fill_manual(values=tech_pal) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45,vjust=.5))
ggsave(file.path("figures","predictive_markers_class_barplot.pdf"),width=2.8,height=3)





source("src/packages_paths.R")

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,best) %>% 
  mutate(.,perturbation="CRISPR")

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,best) %>% 
  mutate(.,perturbation="RNAi")

pred_res <- bind_rows(crispr_res,rnai_res)
pred_res$confounder <- grepl("_Confounders",pred_res$feature0)

pred_res %<>% subset(., (pearson > .5) & !confounder)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
pandeps <- subset(t2,CRISPR_PD)$entrez_id

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
pandeps <- entrez_to_cds(pandeps,hgnc)
pred_res$pandep <- pred_res$gene %in% pandeps

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")
ggplot(subset(pred_res,pearson > .5),aes(x=pandep,fill=perturbation)) +
  geom_bar(stat="count",position="dodge",width=.8,color="black") +
  theme_classic(base_size=11) +
  scale_fill_manual(values=tech_pal) +
  theme(legend.position = "none")
ggsave(file.path("figures","predictive_markers_gw_accurate_model_count_barplot.pdf"),height=2,width=2.25)


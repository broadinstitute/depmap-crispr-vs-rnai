
# CYCLOPS histogram

source("src/packages_paths.R")

rnai_res <- fread("data/processed/ensemble-prediction-summary-rnai-matched.csv") %>%
  subset(.,best) %>% 
  dplyr::select(.,gene,model,pearson,feature0)

rnai_res$status <- "Other"
rnai_res$status[grepl("_Confounder",rnai_res$feature0)] <- "Confounder"

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc)

rnai_res$CRISPR_pandep <- rnai_res$gene %in% subset(t2,CRISPR_PD)$CDS_ID

main_data <- subset(rnai_res,CRISPR_pandep)
main_data <- subset(rnai_res,CRISPR_pandep & (status != "Confounder"))

cyclops_res <- fread(file.path("data/processed","CYCLOPS_annotations.csv"))
cyclops_res$CYCLOPS_uni <- cyclops_res$cn_cor >= .5
multi_genes <- subset(cyclops_res,self_exp_pos | local_cn_pos)$gene
cyclops_res$CYCLOPS_multi <- cyclops_res$gene %in% multi_genes

main_data$CYCLOPS <- "Other"
main_data$CYCLOPS[main_data$gene %in% subset(cyclops_res,CYCLOPS_multi)$gene] <- "CYCLOPS related"
main_data$CYCLOPS[main_data$gene %in% subset(cyclops_res,CYCLOPS_uni)$gene] <- "CYCLOPS"

binary_pal <- c("TRUE"="#C71000B2","FALSE"="white")

group_pal<- c("Other"="#046C9A","CYCLOPS"="#D69C4E","CYCLOPS related"="#ABDDDE")

main_data$CYCLOPS <- factor(main_data$CYCLOPS,levels=c("CYCLOPS related","Other","CYCLOPS"))
ggplot(main_data, aes(pearson, fill = CYCLOPS)) +
  geom_histogram(binwidth = .05,boundary=0,color="black") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=group_pal,name="Dep.-- Marker\nRelationship") +
  ylab("CRISPR Pan-dep. Count") +
  geom_vline(xintercept = .5,linetype="dashed") +
  theme(legend.position = "right") + 
  theme(legend.text=element_text(size=7)) +
  theme(legend.title=element_text(size=9))
ggsave(file.path("figures","predictive_markers_CYCLOPS_or_related_hist.pdf"),height=2,width=3.5)

main_data$CYCLOPS <- "Other"
main_data$CYCLOPS[main_data$gene %in% subset(cyclops_res,CYCLOPS_uni)$gene] <- "CYCLOPS"
main_data$CYCLOPS <- factor(main_data$CYCLOPS,levels=c("Other","CYCLOPS"))
group_pal["Other"] <- "white"

ggplot(main_data, aes(pearson, fill = CYCLOPS)) +
  geom_histogram(binwidth = .05,boundary=0,color="black") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=group_pal) +
  geom_vline(xintercept = .5,linetype="dashed") +
  theme(legend.position = "none")

ggsave(file.path("figures","predictive_markers_CYCLOPS_hist.pdf"),height=2,width=2.25)

nrow(subset(main_data,(pearson > .5) & (CYCLOPS == "CYCLOPS"))) / nrow(subset(main_data,(pearson > .5)))

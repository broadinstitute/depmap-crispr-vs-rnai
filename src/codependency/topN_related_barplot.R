source("src/packages_paths.R")

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

ce_genes <- subset(t2,(CRISPR_class == "Pan-dependency"))$CDS_ID

#Annotate related codep results
results<- fread(file.path(data_processed,"codependency_network_enrichment.csv"))
results$pandep <- results$target %in% ce_genes

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF","SNF"="#669966")

results$pandep_label <- "Other"
results$pandep_label[results$pandep] <- "CRISPR Pan-dep."

top_N <- list("1"=subset(results,related_min_rank_abs <= 1),
              "5"=subset(results,related_min_rank_abs <= 5),
              "10"=subset(results,related_min_rank_abs <= 10),
              "100"=subset(results,related_min_rank_abs <= 100))
for (n in names(top_N)){
  tmp_top = top_N[[n]]
  plot_summary <- expand.grid(tech=unique(results$tech),pandep_label=unique(results$pandep_label))
  plot_summary$count <- 0
  for (i in 1:nrow(plot_summary)){
    plot_summary$count[i] <- nrow(subset(tmp_top,(tech == plot_summary$tech[i]) & (pandep_label == plot_summary$pandep_label[i])))
  }
  plot_summary$top_N <- n
  top_N[[n]] <- plot_summary
}
plot_summary <- bind_rows(top_N)
plot_summary$top_N <- factor(plot_summary$top_N,levels=names(top_N))

ggplot(data=plot_summary,aes(x=top_N,y=count,fill=tech)) +
  geom_bar(position="dodge",stat="identity",color="black") +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  scale_fill_manual(values=tech_pal) +
  scale_color_manual(values=tech_pal) +
  ylab("Number of Genes") +
  xlab("Highest ranking related co-dependency") +
  facet_grid(pandep_label ~ .)
ggsave(file.path("figures","codependency_CRISPR-RNAi-SNF_topN_barplot.pdf"),width=3.5,height=2.8)


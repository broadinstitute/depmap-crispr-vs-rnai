
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
results <- fread(file.path(data_processed,"codependency_network_enrichment.csv"))
results$pandep <- results$target %in% ce_genes

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF","SNF"="#669966")

exponents <- c(1,2,4,8,16)
b <- 10^-(exponents)
l <- paste0("10^{-",exponents,"}")
b <- log10(-log2(b)+1)

results$ks_score <- log10(-log2(results$ks_ABS_pvalue)+1)

results$pandep_label <- "Other"
results$pandep_label[results$pandep] <- "CRISPR Pan-dep."

ggplot(data=results,aes(x=pandep_label,y=ks_score,fill=tech)) +
  geom_boxplot(position="dodge",width=0.8) +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  scale_fill_manual(values=tech_pal) +
  scale_color_manual(values=tech_pal) +
  scale_y_continuous(breaks=b,labels=parse(text=l)) +
  ylab("Co-dependency Enrichment (P-Value)") +
  xlab("")
ggsave(file.path("figures","codependency_enrichment_boxplot_pandep_facet_summary.pdf"),width=2.5,height=2.6)
# 
# wilcox.test(subset(results,pandep & (tech == "CRISPR" ))$ks_score,
#             subset(results,pandep & (tech == "RNAi" ))$ks_score)
# 
# wilcox.test(subset(results,!pandep & (tech == "CRISPR" ))$ks_score,
#             subset(results,!pandep & (tech == "RNAi" ))$ks_score)



source("src/packages_paths.R")

univariate_res <- fread(file.path(data_processed,"drug-screen-genetic-targets-correlations.csv"))
summary <- subset(univariate_res,drug_dataset == "ctd2")
summary %<>% subset(.,best_drug_target_cor)

#Plot the actual pearson correlation for all targets-drug pairs in CRISPR vs RNAi
crispr <- subset(summary,genetic_perturbation == "CRISPR") %>%
  dplyr::select(.,target_drug_pair,pearson_r,rank,dose)
colnames(crispr) <- paste0("CRISPR_",colnames(crispr))
crispr %<>% rename(.,target_drug_pair=CRISPR_target_drug_pair)

rnai <- subset(summary,genetic_perturbation == "RNAi") %>%
  dplyr::select(.,target_drug_pair,pearson_r,rank,dose)
colnames(rnai) <- paste0("RNAi_",colnames(rnai))
rnai %<>% rename(.,target_drug_pair=RNAi_target_drug_pair)

plot_df <- left_join(rnai,crispr,by="target_drug_pair")
plot_df$diff <- plot_df$RNAi_pearson_r - plot_df$CRISPR_pearson_r

plot_df$CDS_ID <- gsub("::.*","",plot_df$target_drug_pair)
plot_df$gene_target <- gsub(" .*","",plot_df$CDS_ID)
crispr_top <- subset(plot_df,(CRISPR_rank <= 1) & (CRISPR_pearson_r > .4))
rnai_top <- subset(plot_df,(RNAi_rank <= 1) &  (RNAi_pearson_r > .4))
top_cors <- intersect(crispr_top$gene_target,rnai_top$gene_target)
top_cors

rnai_highlights <- c("PSMA3","BCL2L1","WEE1","XPO1","CHEK1","AURKA","FGFR2","CDK2")
crispr_highlights <- c("MDM2","NAE1","MET","RARA","JAK2","MAPK14","ABL1","PI4KB")
both_highlights <- c("BRAF")

highlights <- c("BRAF","EGFR","IGF1R","FLT3","ERBB2")
plot_df$label <- plot_df$gene_target
plot_df$label[!(plot_df$label %in% highlights)] <- "other"

dot_pal <- c("other"="grey","BRAF"="#35274A","EGFR"="#0B775E","IGF1R"="#5BBCD6","FLT3"="#F2AD00","ERBB2"="#FF0000")

upper_limit <- ceiling(max(plot_df$RNAi_pearson_r,plot_df$CRISPR_pearson_r) * 10) / 10
lower_limit <- floor(min(plot_df$RNAi_pearson_r,plot_df$CRISPR_pearson_r) * 10) / 10

pandep_pal <- c("TRUE"="red","FALSE"="grey50")

p <- ggplot(plot_df, aes(x=CRISPR_pearson_r,y=RNAi_pearson_r,color=label)) +
  geom_point(data=subset(plot_df,label == "other"),size=1,alpha=.7) +
  geom_point(data=subset(plot_df,label != "other"),size=1,alpha=1) +
  geom_abline(slope=1,linetype="dashed",color="grey50") +
  geom_text_repel(data=subset(plot_df,(diff < -.08) & (gene_target %in% crispr_highlights)),
                  aes(x=CRISPR_pearson_r,y=RNAi_pearson_r,label=gene_target),
                  min.segment.length = 0,
                  color="black",
                  size=2,
                  max.overlaps=20) +
  geom_text_repel(data=subset(plot_df,(diff > .05) & (gene_target %in% rnai_highlights)),
                  aes(x=CRISPR_pearson_r,y=RNAi_pearson_r,label=gene_target),
                  min.segment.length = 0,
                  color="black",
                  size=2,
                  max.overlaps=20) +
  geom_text_repel(data=subset(plot_df,(CRISPR_pearson_r > .7) & (gene_target %in% both_highlights)),
                  aes(x=CRISPR_pearson_r,y=RNAi_pearson_r,label=gene_target),
                  min.segment.length = 0,
                  color="black",
                  size=2,
                  max.overlaps=20) +
  scale_color_manual(values=dot_pal) +
  theme_classic(base_size=11) +
  xlab("CRISPR\nGene-drug Correlation") +
  ylab("RNAi\nGene-drug Correlation") +
  xlim(c(lower_limit,upper_limit))  +
  ylim(c(lower_limit,upper_limit))
ggsave(plot=p+theme(legend.position = "none"),file.path("figures","drug_response_CTD2_correlation_scatter.pdf"),width=3.4,height=3.4)

legend <- cowplot::get_legend(p)
grid::grid.newpage()

pdf(file.path("figures","drug_response_CTD2_correlation_scatter_legend.pdf"))
grid::grid.draw(legend) 
dev.off()



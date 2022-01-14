
source("src/packages_paths.R")

univariate_res <- fread(file.path(data_processed,"drug-screen-genetic-targets-correlations.csv"))

for (dset in c("gdsc","prism")){
  
  summary <- subset(univariate_res,drug_dataset == dset)
  summary %<>% subset(.,best_drug_target_cor)
  
  print(paste0(dset,": ",length(unique(summary$target))," gene targets, ",length(unique(summary$broad_id))," drugs"))
  
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
  
  highlights <- c("MDM2","PIK3CA","BRAF","EGFR","CHEK1","UBA3","BIRC2") #PRISM + GDSC
  plot_df$gene_target <- gsub(" .*","",plot_df$target_drug_pair)
  plot_df$label <- plot_df$gene_target
  plot_df$label[!(plot_df$label %in% highlights)] <- "other"
  
  dot_pal <- c("other"="grey","MDM2"="#35274A","PIK3CA"="#0B775E",
               "BRAF"="#5BBCD6","EGFR"="#F2AD00","CHEK1"="#FF0000",
               "UBA3"="blue","BIRC2"="purple")
  
  plot_df$text_label <- plot_df$gene_target
  plot_df$text_label[abs(plot_df$diff) < .2] <- NA
  plot_df$text_label[plot_df$label != "other"] <- NA
  
  upper_limit <- ceiling(max(plot_df$RNAi_pearson_r,plot_df$CRISPR_pearson_r) * 10) / 10
  lower_limit <- floor(min(plot_df$RNAi_pearson_r,plot_df$CRISPR_pearson_r) * 10) / 10
  
  p <- ggplot(plot_df, aes(x=CRISPR_pearson_r,y=RNAi_pearson_r,color=label)) +
    geom_point(data=subset(plot_df,label == "other"),size=.7,alpha=.7) +
    geom_abline(slope=1,linetype="dashed",color="grey50") +
    geom_point(data=subset(plot_df,label != "other"),size=1) +
    geom_text_repel(data=subset(plot_df,!is.na(plot_df$text_label)),
                    aes(x=CRISPR_pearson_r,y=RNAi_pearson_r,label=text_label),
                    min.segment.length = 0,
                    color="black",
                    size=2,
                    max.overlaps=20) +
    scale_color_manual(values=dot_pal,name="Gene Target") +
    theme_classic(base_size=11) +
    xlab("CRISPR\nGene-drug Correlation") +
    ylab("RNAi\nGene-drug Correlation") +
    xlim(c(lower_limit,upper_limit)) +
    ylim(c(lower_limit,upper_limit))
  ggsave(plot=p + theme(legend.position = "none"),file.path("figures",paste0("drug_response_",dset,"_correlation_scatter.pdf")),width=2.5,height=2.5)
  
  legend <- cowplot::get_legend(p)
  grid::grid.newpage()
  
  pdf(file.path("figures",paste0("drug_response_",dset,"_correlation_scatter_legend.pdf")))
  grid::grid.draw(legend) 
  dev.off()
  
}

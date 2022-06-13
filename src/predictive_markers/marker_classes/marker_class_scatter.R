

bdp_classes <- list(cyclops="CYCLOPS.csv",
                    expAddict="expression_addiction.csv",
                    driver="genetic_driver.csv",
                    onco="oncogene_addiction.csv",
                    oncoSL="oncoKB_synthetic_lethal.csv",
                    paralog="paralog.csv")
bdp_classes <- lapply(bdp_classes,function(x){fread(file.path(data_processed,paste0("predictive_marker_class_",x)))})

bdp_classes <- bind_rows(bdp_classes)
bdp_classes$class <- gsub(" \\(oncogene\\)","",bdp_classes$class)
bdp_classes$class <- gsub(" \\(TSG\\)","",bdp_classes$class)

hit_list <- list("genetic driver"=subset(bdp_classes,class == "genetic driver")$gene,
                 "expression addiction"=subset(bdp_classes,class == "expression addiction")$gene,
                 "paralog"=subset(bdp_classes,class == "paralog")$gene,
                 "CYCLOPS"=subset(bdp_classes,class == "CYCLOPS")$gene,
                 "oncogene addiction"=subset(bdp_classes,class == "oncogene addiction")$gene,
                 "synthetic lethal"=subset(bdp_classes,class == "synthetic lethal")$gene)

crispr_res <- fread(file.path(data_processed,"ensemble-prediction-summary-crispr-matched.csv")) %>%
  subset(.,best) %>% 
  dplyr::select(.,gene,pearson) %>%
  rename(CRISPR_pearson = pearson)

rnai_res <- fread(file.path(data_processed,"ensemble-prediction-summary-rnai-matched.csv")) %>%
  subset(.,best) %>% 
  dplyr::select(.,gene,pearson) %>%
  rename(RNAi_pearson = pearson)

ensemble <- left_join(crispr_res,rnai_res,by="gene")

for (g in names(hit_list)){
  plot_df <- subset(ensemble,gene %in% hit_list[[g]])
  plot_df$label <- gsub(" .*","",plot_df$gene)
  plot_df$CRISPR_top <- rank(-plot_df$CRISPR_pearson)
  plot_df$RNAi_top <- rank(-plot_df$RNAi_pearson)
  
  plot_df$CRISPR_bottom <- rank(plot_df$CRISPR_pearson)
  plot_df$RNAi_bottom <- rank(plot_df$RNAi_pearson)
  
  plot_df$CRISPR_label <- (plot_df$CRISPR_top <= 5) | (plot_df$CRISPR_bottom <= 5)
  plot_df$RNAi_label <- (plot_df$RNAi_top <= 5) | (plot_df$RNAi_bottom <= 5)
  
  plot_df$label[!(plot_df$CRISPR_label | plot_df$RNAi_label)] <- NA
  
  min_limit <- min(c(plot_df$CRISPR_pearson,plot_df$RNAi_pearson))
  max_limit <- max(c(plot_df$CRISPR_pearson,plot_df$RNAi_pearson))
  
  p <- ggplot(plot_df, aes(x=CRISPR_pearson,y=RNAi_pearson)) +
    geom_abline(intercept=0,slope=1,linetype="dashed",color="grey50") +
    geom_point(color="darkblue") +
    xlim(c(min_limit,max_limit)) +
    ylim(c(min_limit,max_limit)) +
    geom_text_repel(data=subset(plot_df,!is.na(label)),aes(x=CRISPR_pearson,y=RNAi_pearson,label=label),
                    size=3,
                    min.segment.length=0) +
    theme_bw(base_size=11) +
    xlab("CRISPR Predictive Accuracy (r)") +
    ylab("RNAi Predictive Accuracy (r)") +
    ggtitle(g)
  ggsave(plot=p,file.path("figures",paste0("predictive_markers_class","_scatter_",gsub(" ","_",g),".pdf")),width=3.3,height=2.7)
}


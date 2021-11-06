library(ggrepel)
#requires matrixStats

data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")

results_list <- list()

recall <- fread(file.path(data_processed,"CRISPR-Avana_vs_CRISPR-KY_cor_recall_vals.csv"))
recall$entrez_id %<>% as.character(.)
comp_results <- recall
comp_results$rank <- matrixStats::rowMins(as.matrix(select(comp_results,`rank_in_CRISPR-Avana`,`rank_in_CRISPR-KY`)),na.rm=T)
comp_results %<>% select(.,entrez_id,rank)
results_list[["CRISPR"]] <- comp_results

recall <- fread(file.path(data_processed,"RNAi-Achilles_vs_RNAi-DRIVE_cor_recall_vals.csv"))
recall$entrez_id %<>% as.character(.)
comp_results <- recall
comp_results$rank <- matrixStats::rowMins(as.matrix(select(comp_results,`rank_in_RNAi-Achilles`,`rank_in_RNAi-DRIVE`)),na.rm=T)
comp_results %<>% select(.,entrez_id,rank)
results_list[["RNAi"]] <- comp_results


#### Get genes in union of each dependency class

plot_df <- list()
for (comp in names(results_list)){
  data <- results_list[[comp]]
  
  all_plot_dfs <- list()
  tmp_data <- data
  h <- table(tmp_data$rank)
  tmp_df <- data.frame(breaks=as.numeric(names(h)),counts=as.numeric(h),tech=comp,stringsAsFactors = F)
  tmp_df$Fraction <- tmp_df$counts / sum(tmp_df$counts)
  tmp_df$cumfrac <- cumsum(tmp_df$Fraction)
  last <- tmp_df[nrow(tmp_df),]
  last[1,"breaks"] <- max(data$rank)
  last[1,"cumfrac"] <- 1
  tmp_df <- rbind(tmp_df,last)

  plot_df[[comp]] <- tmp_df
}
  plot_df <- bind_rows(plot_df)

  #plot_df$ranks <- log10(log2(plot_df$breaks)+1)
  plot_df$ranks <- log10((plot_df$breaks))
  
  b <- c(1,10,100,1000,10000)
  l <- c("1","10","10^{2}","10^{3}","10^{4}")  
  b <- log10((b))
  
  #plot_df$Group <- factor(plot_df$Group,levels=c("SSD","Selective","Pan-dependency","No dependency"))
  
  mypal = c("No dependency"="#00204DFF","SSD"="#145A32","Pan-dependency"="#BB3754FF","Selective"="Gray")
  
  crispr_df <- subset(plot_df,tech == "CRISPR")
  crispr_total <- subset(crispr_df,breaks == 1)$counts
  
  p <- ggplot(crispr_df, aes(x=ranks,y=cumfrac)) +
    geom_line(size=.85,color="gray") +
    geom_point(data=subset(crispr_df,ranks == 0),size=1,color="#145A32") +
    geom_label_repel(data=subset(crispr_df,ranks == 0),label=paste0(crispr_total," genes"),xlim=c(.5,NA),ylim=c(.75,NA)) +
    scale_x_continuous(breaks=b,labels=parse(text=l)) +
    theme_classic(base_size=9) +
    xlab("Rank") +
    ylab("Gene Frac.") +
    theme(legend.position = "none") 
  ggsave(plot=p,file.path("figures","high_confidence_crispr_total_recall_curve.pdf"),width=1.25,height=1.25)
  
  rnai_df <- subset(plot_df,tech == "RNAi")
  rnai_total <- subset(rnai_df,breaks == 1)$counts
  
  p <- ggplot(rnai_df, aes(x=ranks,y=cumfrac)) +
    geom_line(size=.85,color="gray") +
    geom_point(data=subset(rnai_df,ranks == 0),size=1,color="#145A32") +
    geom_label_repel(data=subset(rnai_df,ranks == 0),label=paste0(rnai_total," genes"),xlim=c(.5,NA),ylim=c(.75,NA)) +
    scale_x_continuous(breaks=b,labels=parse(text=l)) +
    theme_classic(base_size=9) +
    xlab("Rank") +
    ylab("Gene Frac.") +
    theme(legend.position = "none") 
  ggsave(plot=p,file.path("figures","high_confidence_rnai_total_recall_curve.pdf"),width=1.25,height=1.25)
  
  

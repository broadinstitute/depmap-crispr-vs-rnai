
library(tidyverse)
library(magrittr)
library(data.table)

library(cowplot)
library(grid)
library(gridExtra)

source("src/id_utility.R")

gene_scores <- list(all=list(crispr_avana = fread("data/raw/gene-effect-scaled-crispr-avana.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                    crispr_ky = fread("data/raw/gene-effect-scaled-crispr-ky.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                    rnai_achilles = fread("data/raw/gene-effect-scaled-rnai-achilles.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                    rnai_drive = fread("data/raw/gene-effect-scaled-rnai-drive.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.)))

ce_percentile <- function(data,quantile_cutoff){

  data <- t(data)
  
  rank_data <- plyr::aaply(data, .margins=2, rank, ties.method="min")
  rank_data %<>% t(.)
  rank_data[is.na(data)] <- NA
  rank_data %<>% t(.)
  rank_data <- rank_data / colSums(!is.na(data))
  percentile <- plyr::aaply(rank_data, .margins=2, quantile, probs=quantile_cutoff, na.rm=T)
  
  dens <- density(percentile)
  df <- data.frame(x=dens$x, y=dens$y)
  df_range <- subset(df,df$x > .1 & df$x < .8)
  threshold <- df_range$x[which.min(df_range$y)]
  
  ce_list <- list("Gene"=names(percentile),"CE_percentile"=percentile)
  ce_list[["Common_Essential"]] <- percentile <= threshold
  ce_list[["threshold"]] <- threshold
  
  return(ce_list)
}



ranks <- list()
for (gp in names(gene_scores)){
  ranks[[gp]] <- list()
  for (dset in names(gene_scores[[gp]])){
    ranks[[gp]][[dset]] <- ce_percentile(data=gene_scores[[gp]][[dset]],quantile_cutoff=.9)
  }
}


make_df <- function(ranks_list,list_name){
  
  df <- data.frame(CDS_ID=ranks_list$Gene,
                   rank=ranks_list$CE_percentile,
                   CE=ranks_list$Common_Essential,
                   rank_threshold=ranks_list$threshold,
                   dataset=list_name[1],
                   stringsAsFactors = F)
  return(df)
}

for (gp in names(ranks)){
  for (dname in names(ranks[[gp]])){
    ranks[[gp]][[dname]] <- make_df(ranks_list=ranks[[gp]][[dname]],list_name=dname)
  }
}

results_all <- bind_rows(ranks[["all"]])
write_csv(results_all,"data/processed/multilib_ce_percentile_results.csv")

##########  Hists for all 4 datasets run independently 

make_ce_hist <- function(ranks_output,tech_name,ylim_upper,bwidth){
  pandep_count <- nrow(subset(ranks_output,CE))
  plot_df <- data.frame(values=ranks_output$rank,tech=tech_name,threshold=ranks_output$rank_threshold,stringsAsFactors = F)
  plot_df %<>% remove_rownames(.)
  ce_bins <- round(ranks_output$rank_threshold[1] / bwidth )
  other_bins <- (1 / bwidth ) - ce_bins + 1
  
  mypal <- c(rep("#bd0026",ce_bins),rep("#fff7ec",other_bins))
  
  p <- ggplot(data=plot_df,aes(values)) +
    geom_histogram(col="black",binwidth = bwidth,fill=mypal) +
    annotate("text",x=.2,y=500,label=pandep_count) +
    scale_x_continuous(breaks=c(0,.2,.4,.6,.8,1)) +
    scale_y_continuous(breaks=c(0,300,600,900,1200,1500,1800),limits=c(0,ylim_upper)) +
    xlab("") +
    ylab("") +
    ggtitle(tech_name) +
    theme_classic(base_size = 10) +
    theme(plot.margin = margin(l=-.5,r=0,unit="cm"))
  
  return(p)
}


plot_pairs <- list("CRISPR"=list("Avana"=ranks$all$crispr_avana,"KY"=ranks$all$crispr_ky),
                   "RNAi"=list("Achilles"=ranks$all$rnai_achilles,"DRIVE"=ranks$all$rnai_drive))

for (p in names(plot_pairs)){
  
  pnames <- names(plot_pairs[[p]])
  
  bwidth <- .025
  ylim <- max(unlist(lapply(c(1,2),function(x){max(hist(plot_pairs[[p]][[pnames[x]]]$rank,breaks=(1/bwidth),plot = F)$counts)}))) + 300
  
  plots <- lapply(c(1,2),function(x){make_ce_hist(plot_pairs[[p]][[pnames[x]]],tech_name = pnames[x],ylim_upper=ylim,bwidth=bwidth)})
  
  p1 <- plot_grid(plots[[1]],plots[[2]],nrow=1)
  p2 <- add_sub(p1, "Ranking within 90th Percentile Least Dependent Cell Line", size = 11, vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5)
  p3 <- grid.arrange(textGrob("Number of Genes",x = unit(.5, "npc"),rot=90,gp=gpar(fontsize=11)), p2, ncol=2,widths=c(1,20))
  ggsave(plot=p3,paste0("figures/pandependency_90th_percentile_ranks_",p,".pdf"),device="pdf",height=2.5,width=4.5)
  
  
}



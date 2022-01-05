set.seed(42)

source("src/packages_paths.R")

library(ggrepel)

# require(fgsea)
# require(cowplot)
# require(gridExtra)

#get LRTs
crispr_lrt <- fread(file.path(data_processed,"crispr-matched-LRT.csv")) %>% rename(CRISPR_LRT=LRT)
rnai_lrt <- fread(file.path(data_processed,"rnai-matched-LRT.csv")) %>% rename(RNAi_LRT=LRT)

master_table <- inner_join(crispr_lrt,rnai_lrt,by="Row.name")
master_table %<>% add_column(entrez_id=extract_entrez(master_table$Row.name),.before=1)
master_table %<>% add_column(symbol=gsub(" .*","",master_table$Row.name),.before=2)

#assign LRT rank
master_table <- master_table[order(master_table$CRISPR_LRT,decreasing=T),]
master_table$CRISPR_LRT_rank <- 1:nrow(master_table) 
master_table <- master_table[order(master_table$RNAi_LRT,decreasing=T),]
master_table$RNAi_LRT_rank <- 1:nrow(master_table) 
master_table$entrez_id %<>% as.character(.)

#get oncogenes and TSGs from oncoKB
oncoKB <- fread(file.path(data_raw,"control-oncoKB.csv"))
oncogene <- subset(oncoKB,source =="oncogene")
oncogene$entrez_id <- extract_entrez(oncogene$CDS_ID)
tsg <- subset(oncoKB,source =="TSG")
tsg$entrez_id <- extract_entrez(tsg$CDS_ID)

#Filter for oncogenes with missense mutation or TSG with damaging mutation
mut_other <- fread(file.path(data_raw,"depmap-omics-mutation-hotspot.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)

#filter for cell lines used in D2-Combined overlap
crispr.gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
mut_other <- mut_other[rownames(mut_other) %in% rownames(crispr.gs),]
colnames(mut_other) <- extract_entrez(colnames(mut_other))
mut_other <- mut_other[,(colnames(mut_other) %in% oncogene$entrez_id)]
mut_other_count <- colSums(mut_other,na.rm=T)
keep_onco <- names(mut_other_count)[mut_other_count > 0]
oncogene %<>% subset(.,entrez_id %in% keep_onco)

mut_dam <- fread(file.path(data_raw,"depmap-omics-mutation-damaging.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
mut_dam <- mut_dam[rownames(mut_dam) %in% rownames(crispr.gs),]
colnames(mut_dam) <- extract_entrez(colnames(mut_dam))
mut_dam <- mut_dam[,(colnames(mut_dam) %in% tsg$entrez_id)]
mut_dam_count <- colSums(mut_dam,na.rm=T)
keep_tsg <- names(mut_dam_count)[mut_dam_count > 0]
tsg %<>% subset(.,entrez_id %in% keep_tsg)

master_table$oncoKB <- "Other"
master_table$oncoKB[master_table$entrez_id %in% oncogene$entrez_id] <- "oncogene"
master_table$oncoKB[master_table$entrez_id %in% tsg$entrez_id] <- "TSG"


onkoKB_sets <- list(oncogenes=oncogene$entrez_id,
                    TSGs=tsg$entrez_id)

customEnrichment <- function(set_ids,stats,labels_df,plot_title,gseaParam=1,ticksSize=1,num_labels=10,label_font=8){
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(set_ids, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 3
  labels_df$y <- -diff
  labels_df <- labels_df[order(labels_df$x),]
  
  #Keep only the specified number of labels
  labels_df <- labels_df[1:num_labels,]
  spacing <- floor(n / num_labels)
  labels_df$label_x <- cumsum(rep(spacing,num_labels)) - round(spacing / 2)
  labels_df$label_y <- -.6
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_hline(yintercept=max(tops), colour="grey", linetype="dashed") +
    geom_segment(mapping=aes(x=min(toPlot$x),xend=max(toPlot$x),y=min(bottoms),yend=min(bottoms)),color="grey") +
    geom_segment(mapping=aes(x=min(toPlot$x),xend=max(toPlot$x),y=-diff,yend=-diff),color="grey") +
    geom_segment(mapping=aes(x=max(toPlot$x),xend=max(toPlot$x),y=-diff,yend=min(bottoms)),color="grey") +
    geom_segment(mapping=aes(x=min(toPlot$x),xend=min(toPlot$x),y=-diff,yend=min(bottoms)),color="grey") +
    
    geom_line(color="pink",size=.5) + 
    theme_classic(base_size=11) +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff,
                             xend=x, yend=0),
                 size=ticksSize) +
    scale_y_continuous(breaks=c(0,.5,1),labels=c("0","0.5","1"),limits=c(-1.5,1)) +
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.title.y = element_blank()
    ) +
    geom_label(data = labels_df,mapping=aes(x=label_x,y=label_y,label=symbol,hjust=1),size=(label_font / ggplot2:::.pt)) +
    geom_segment(data=labels_df,
                 mapping=aes(x=label_x, y=label_y,
                             xend=x, yend=y),
                 size=.25,
                 color="grey50") +
    labs(x="rank", y="enrichment score") +
    coord_flip() +
    scale_x_reverse(position="top") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size=11,hjust = 0.5)) +
    theme(axis.text.x = element_text(size=9),
          axis.title.x = element_text(size=11)) 
  
  
  return(g)
}

####### CRISPR #######
crispr_lrt <- master_table$CRISPR_LRT
names(crispr_lrt) <- master_table$entrez_id
crispr_lrt <- crispr_lrt[!is.na(crispr_lrt)]
crispr_fgsea <- fgsea::fgsea(pathways = onkoKB_sets,
                      stats = crispr_lrt,
                      minSize=15,
                      maxSize=500,
                      nperm=1000000)

toPlot_labels <- subset(master_table,entrez_id %in% onkoKB_sets[["oncogenes"]])
toPlot_labels %<>% dplyr::select(.,x=CRISPR_LRT_rank,symbol)

crispr_onco_plot <- customEnrichment(set_ids=onkoKB_sets[["oncogenes"]],
                                     stats=crispr_lrt,
                                     labels_df=toPlot_labels,
                                     plot_title = "CRISPR",
                                     ticksSize = .5,
                                     num_labels=10)

toPlot_labels <- subset(master_table,entrez_id %in% onkoKB_sets[["TSGs"]])
toPlot_labels %<>% dplyr::select(.,x=CRISPR_LRT_rank,symbol)

crispr_tsg_plot <- customEnrichment(set_ids=onkoKB_sets[["TSGs"]],
                                    stats=crispr_lrt,
                                    labels_df=toPlot_labels,
                                    plot_title = "CRISPR",
                                    ticksSize = .5)
####### RNAi #######
rnai_lrt <- master_table$RNAi_LRT
names(rnai_lrt) <- master_table$entrez_id
rnai_lrt <- rnai_lrt[!is.na(rnai_lrt)]
rnai_fgsea <- fgsea::fgsea(pathways = onkoKB_sets,
                    stats = rnai_lrt,
                    minSize=15,
                    maxSize=500,
                    nperm=1000000)


toPlot_labels <- subset(master_table,entrez_id %in% onkoKB_sets[["oncogenes"]])
toPlot_labels %<>% dplyr::select(.,x=RNAi_LRT_rank,symbol)

rnai_onco_plot <- customEnrichment(set_ids=onkoKB_sets[["oncogenes"]],
                                   stats=rnai_lrt,
                                   labels_df=toPlot_labels,
                                   plot_title = "RNAi",
                                   ticksSize = .5,
                                   num_labels=10)

toPlot_labels <- subset(master_table,entrez_id %in% onkoKB_sets[["TSGs"]])
toPlot_labels %<>% dplyr::select(.,x=RNAi_LRT_rank,symbol)

rnai_tsg_plot <- customEnrichment(set_ids=onkoKB_sets[["TSGs"]],
                                  stats=rnai_lrt,
                                  labels_df=toPlot_labels,
                                  plot_title = "RNAi",
                                  ticksSize = .5)

####### Compile plots ########

onco_row <- cowplot::plot_grid(crispr_onco_plot,rnai_onco_plot,
                      align = 'v',
                      nrow=1,
                      rel_widths=c(1,1),
                      axis="tb",
                      label_size = 11)
onco_row_final <- gridExtra::arrangeGrob(grid::textGrob("LRT Rank",vjust = -1,x = unit(1, "npc"),rot=90,gp=gpar(fontsize=11)), onco_row, ncol=2,widths=c(2,20))
ggsave(plot=onco_row_final,file=file.path("figures","gene_deps_profiles_LRT_onco_enrichment.pdf"),height=3,width=3.5)

tsg_row <- cowplot::plot_grid(crispr_tsg_plot,rnai_tsg_plot,
                     align = 'v',
                     nrow=1,
                     rel_widths=c(1,1),
                     axis="tb",
                     label_size = 11)
tsg_row_final <- gridExtra::arrangeGrob(grid::textGrob("LRT Rank",vjust = -1,x = unit(1, "npc"),rot=90,gp=gpar(fontsize=11)), tsg_row, ncol=2,widths=c(2,20))
ggsave(plot=tsg_row_final,file=file.path("figures","gene_deps_profiles_LRT_tsg_enrichment.pdf"),height=3,width=3.5)



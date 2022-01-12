
# require(scales)
# require(reshape2)

source("src/packages_paths.R")

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc)

gene_sets <- list()

gene_sets[["crispr_pd"]] <- subset(t2,CRISPR_class == "Pan-dependency")$CDS_ID
gene_sets[["rnai_pd"]] <- subset(t2,RNAi_class == "Pan-dependency")$CDS_ID

gene_sets[["crispr_ssd"]] <- subset(t2,CRISPR_class == "Strongly Selective Dependency")$CDS_ID
gene_sets[["rnai_ssd"]] <- subset(t2,RNAi_class == "Strongly Selective Dependency")$CDS_ID

gene_sets[["crispr_hv"]] <- subset(t2,CRISPR_class == "High-variance Dependency")$CDS_ID
gene_sets[["rnai_hv"]] <- subset(t2,RNAi_class == "High-variance Dependency")$CDS_ID

#import probabilities for 19Q1-public and D2-combined overlap
crispr_pr <- fread(file.path(data_processed,"dependency-probability-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name")
rnai_pr <- fread(file.path(data_processed,"dependency-probability-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name")

#Total percent of deps per CL
totals <- rowSums(!is.na(crispr_pr))
crispr_hits <- rowSums(crispr_pr > .5,na.rm=T)
crispr_frac <- crispr_hits / totals
rnai_hits <- rowSums(rnai_pr > .5,na.rm=T)
rnai_frac <- rnai_hits / totals

paste0(round(mean(crispr_frac),digits=3)," +/- ",round(sd(crispr_frac),digits=3))
paste0(round(mean(rnai_frac),digits=3)," +/- ",round(sd(rnai_frac),digits=3))

#Sum the hits per cell line
cls_hits <- list()

for (set_name in names(gene_sets)){
  dset <- strsplit(set_name,"_")[[1]][1]
  gtype <- strsplit(set_name,"_")[[1]][2]
  if (gtype %in% c("pd","ssd","hv")){
    if (dset == "crispr"){
      
      totals <- rowSums(!is.na(crispr_pr))
      dep_totals <- rowSums(crispr_pr > .5,na.rm=T)
      hits <- rowSums(crispr_pr[,gene_sets[[set_name]]] > .5,na.rm=T)
      frac <- hits/totals
      dep_frac <- hits/dep_totals
      
      cls_hits[[set_name]] <- data.frame(CL=rownames(crispr_pr),
                                         hits=frac,
                                         group=set_name,
                                         screened=totals,
                                         dep=dep_totals,
                                         frac_of_deps=dep_frac,
                                         stringsAsFactors = F)
      
    } else {
      
      totals <- rowSums(!is.na(rnai_pr))
      dep_totals <- rowSums(rnai_pr > .5,na.rm=T)
      hits <- rowSums(rnai_pr[,gene_sets[[set_name]]] > .5,na.rm=T)
      frac <- hits/totals
      dep_frac <- hits/dep_totals
      
      cls_hits[[set_name]] <- data.frame(CL=rownames(rnai_pr),
                                         hits=frac,
                                         group=set_name,
                                         screened=totals,
                                         dep=dep_totals,
                                         frac_of_deps=dep_frac,
                                         stringsAsFactors = F)
    }
  } 
  
}



########################### total ########################### 

hits_long <- bind_rows(cls_hits)
hits_wide <- bind_cols(cls_hits)

hits_mat <- hits_wide[,grepl("hits",colnames(hits_wide))] #choice of total frac vs dep frac
colnames(hits_mat) <- names(gene_sets)
rownames(hits_mat) <- rownames(hits_wide)
hits_mat$ce_avg <- rowMeans(as.matrix(hits_mat[,c("crispr_pd","rnai_pd")]))
hits_mat <- hits_mat[order(hits_mat$ce_avg,decreasing=F),]

lineage_info <- fread(file.path(data_raw,"depmap-cell-line-annotations-v131.csv"))
lineage_info %<>% dplyr::select(.,DepMap_ID,`Primary Disease`,Disease)
lineage_info %<>% subset(.,DepMap_ID %in% hits_long$CL)
lineage_info %<>% unique(.)
lineage_info$Groups <- gsub(" Cancer","",lineage_info$`Primary Disease`)
lineage_info$Groups <- gsub("/.*","",lineage_info$Groups)
lineage_info %<>% dplyr::select(.,CL=DepMap_ID,Groups)
lineage_info$Groups %<>% as.character(.)

#lineage CL counts
lin_counts <- table(lineage_info$Groups)
lin_counts_label <- data.frame(Groups=names(lin_counts),count=as.numeric(lin_counts),GroupsN=paste0(names(lin_counts)," (",lin_counts,")"),stringsAsFactors = F)
lin_counts_label <- lin_counts_label[order(lin_counts_label$count,decreasing=T),]
lineage_info %<>% left_join(.,lin_counts_label,by="Groups")

hits_long %<>% left_join(.,lineage_info,by="CL")
hits_long$GroupsN <- factor(hits_long$GroupsN,levels=lin_counts_label$GroupsN)

hits_long$tech <- "RNAi"
hits_long$tech[grep("crispr_",hits_long$group)] <- "CRISPR"

hits_long %<>% subset(.,!grepl("_ne",group))

hits_long$gene_type <- "SSD"
hits_long$gene_type[grep("_pd",hits_long$group)] <- "Pan-dependency"
hits_long$gene_type[grep("_hv",hits_long$group)] <- "High-variance"
hits_long$gene_type <- factor(hits_long$gene_type,levels=c("Pan-dependency","High-variance","SSD"))

hits_long$CL <- factor(hits_long$CL,levels=rownames(hits_mat))

bars_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

total_p <- ggplot(hits_long, aes(CL, hits, color = tech)) +
  geom_line(color='gray',aes(group=CL),alpha=.7,size=.25) +
  geom_point(size=.4) +
  scale_x_discrete(expand=c(0,3),breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  facet_grid(gene_type ~ GroupsN,scales="free",space="free_x",switch="x") +
  theme_bw(base_size=11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x= element_blank()) +
  theme(panel.spacing = unit(.15, "lines")) +
  scale_color_manual(values=bars_pal,name="") + 
  xlab("") +
  ylab("Percent of Total Genes Profiled") +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(angle = 45),
        strip.background=element_blank(),
        plot.margin=unit(c(2,2,-12,2), "mm"))
ggsave(plot=total_p,file.path("figures","gene_deps_profiles_totalfrac_per_CL_byDisease.pdf"),height=4.5,width=7)

#Get the mean and sd of dep frac of total screened 
crispr_ce_stats <- subset(hits_long,(gene_type == "Pan-dependency") & (tech == "CRISPR"))
rnai_ce_stats <- subset(hits_long,(gene_type == "Pan-dependency") & (tech == "RNAi"))
paste0(round(mean(crispr_ce_stats$hits),digits=3)," +/- ",round(sd(crispr_ce_stats$hits),digits=3))
paste0(round(mean(rnai_ce_stats$hits),digits=3)," +/- ",round(sd(rnai_ce_stats$hits),digits=3))

paste0(round(mean(crispr_ce_stats$frac_of_deps),digits=3)," +/- ",round(sd(crispr_ce_stats$frac_of_deps),digits=3))
paste0(round(mean(rnai_ce_stats$frac_of_deps),digits=3)," +/- ",round(sd(rnai_ce_stats$frac_of_deps),digits=3))

#Get the mean and sd of dep frac of total deps
crispr_ssd_stats <- subset(hits_long,(gene_type == "SSD") & (tech == "CRISPR"))
rnai_ssd_stats <- subset(hits_long,(gene_type == "SSD") & (tech == "RNAi"))
paste0(round(mean(crispr_ssd_stats$hits),digits=4)," +/- ",round(sd(crispr_ssd_stats$hits),digits=4))
paste0(round(mean(rnai_ssd_stats$hits),digits=4)," +/- ",round(sd(rnai_ssd_stats$hits),digits=4))

paste0(round(mean(crispr_ssd_stats$frac_of_deps),digits=4)," +/- ",round(sd(crispr_ssd_stats$frac_of_deps),digits=4))
paste0(round(mean(rnai_ssd_stats$frac_of_deps),digits=4)," +/- ",round(sd(rnai_ssd_stats$frac_of_deps),digits=4))

group_frac <- list()
group_frac_of_deps <- list()
for (lin in unique(hits_long$Groups)){
  tmp_crispr <- subset(crispr_ssd_stats,Groups == lin)
  tmp_rnai <- subset(rnai_ssd_stats,Groups == lin)
  
  group_frac[[lin]] <- c(Group=lin,
                         CRISPR_u=mean(tmp_crispr$hits),
                         CRISPR_sd=sd(tmp_crispr$hits),
                         RNAi_u=mean(tmp_rnai$hits),
                         RNAi_sd=sd(tmp_rnai$hits))
  
  group_frac_of_deps[[lin]] <- c(Group=lin,
                                 CRISPR_u=mean(tmp_crispr$frac_of_deps),
                                 CRISPR_sd=sd(tmp_crispr$frac_of_deps),
                                 RNAi_u=mean(tmp_rnai$frac_of_deps),
                                 RNAi_sd=sd(tmp_rnai$frac_of_deps))
}
group_ssd_frac <- bind_rows(group_frac)
colnames(group_ssd_frac) <- c("Group","CRISPR_u","CRISPR_sd","RNAi_u","RNAi_sd")
group_ssd_frac %<>% as.data.frame(.)
group_ssd_frac$CRISPR_u %<>% as.character(.)
group_ssd_frac$RNAi_u %<>% as.character(.)
group_ssd_frac$CRISPR_u %<>% as.numeric(.)
group_ssd_frac$RNAi_u %<>% as.numeric(.)
group_ssd_frac$diff <- group_ssd_frac$CRISPR_u - group_ssd_frac$RNAi_u

########################### dep ########################### 

hits_long <- bind_rows(cls_hits)
hits_wide <- bind_cols(cls_hits)

hits_mat <- hits_wide[,grepl("frac_of_deps",colnames(hits_wide))] #choice of total frac vs dep frac
colnames(hits_mat) <- names(gene_sets)
rownames(hits_mat) <- rownames(hits_wide)
hits_mat$ce_avg <- rowMeans(as.matrix(hits_mat[,c("crispr_pd","rnai_pd")]))
hits_mat <- hits_mat[order(hits_mat$ce_avg,decreasing=F),]

lineage_info <- fread(file.path(data_raw,"depmap-cell-line-annotations.csv"))
lineage_info %<>% dplyr::select(.,DepMap_ID,`Primary Disease`,Disease)
lineage_info %<>% subset(.,DepMap_ID %in% hits_long$CL)
lineage_info %<>% unique(.)
lineage_info$Groups <- gsub(" Cancer","",lineage_info$`Primary Disease`)
lineage_info$Groups <- gsub("/.*","",lineage_info$Groups)
lineage_info %<>% dplyr::select(.,CL=DepMap_ID,Groups)
lineage_info$Groups %<>% as.character(.)

#lineage CL counts
lin_counts <- table(lineage_info$Groups)
lin_counts_label <- data.frame(Groups=names(lin_counts),count=as.numeric(lin_counts),GroupsN=paste0(names(lin_counts)," (",lin_counts,")"),stringsAsFactors = F)
lin_counts_label <- lin_counts_label[order(lin_counts_label$count,decreasing=T),]
lineage_info %<>% left_join(.,lin_counts_label,by="Groups")

hits_long %<>% left_join(.,lineage_info,by="CL")
hits_long$GroupsN <- factor(hits_long$GroupsN,levels=lin_counts_label$GroupsN)

hits_long$tech <- "RNAi"
hits_long$tech[grep("crispr_",hits_long$group)] <- "CRISPR"

hits_long %<>% subset(.,!grepl("_ne",group))

hits_long$gene_type <- "SSD"
hits_long$gene_type[grep("_pd",hits_long$group)] <- "Pan-dependency"
hits_long$gene_type[grep("_hv",hits_long$group)] <- "High-variance"
hits_long$gene_type <- factor(hits_long$gene_type,levels=c("Pan-dependency","High-variance","SSD"))

hits_long$CL <- factor(hits_long$CL,levels=rownames(hits_mat))

bars_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

dep_p <- ggplot(hits_long, aes(x=CL, frac_of_deps, color = tech)) +
  geom_line(color='gray',aes(group=CL),alpha=.7,size=.25) +
  geom_point(size=.4) +
  scale_x_discrete(expand=c(0,3),breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  facet_grid(gene_type ~ GroupsN,scales="free",space="free_x",switch="x") +
  theme_bw(base_size=11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x= element_blank()) +
  theme(panel.spacing = unit(.15, "lines")) +
  scale_color_manual(values=bars_pal,name="") + 
  xlab("") +
  ylab("Percent of Cell Line Dependencies") +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(angle = 45),
        strip.background=element_blank(),
        plot.margin=unit(c(2,2,-12,2), "mm"))
ggsave(plot=dep_p,file.path("figures","gene_deps_profiles_depfrac_per_CL_byDisease.pdf"),height=4.5,width=7)

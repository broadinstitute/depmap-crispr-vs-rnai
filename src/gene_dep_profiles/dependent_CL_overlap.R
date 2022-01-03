
source("src/packages_paths.R")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc=hgnc)

crispr_pr <- fread(file.path(data_processed,"dependency-probability-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
rnai_pr <- fread(file.path(data_processed,"dependency-probability-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)

crispr_hit <- crispr_pr >= .5
rnai_hit <- rnai_pr >= .5

total_rnai_hits <- colSums(rnai_hit,na.rm=T)
total_crispr_hits <- colSums(crispr_hit,na.rm=T)
conf_hit <- crispr_hit & rnai_hit
total_conf_hits <- colSums(conf_hit,na.rm=T)

conf_frac <- total_conf_hits / total_rnai_hits

dep_counts <- data.frame(CDS_ID=names(total_crispr_hits),
                         CRISPR_hit=total_crispr_hits,
                         RNAi_hit=total_rnai_hits,
                         RNAi_conf_hit=total_conf_hits,
                         RNAi_conf_frac=conf_frac,
                         stringsAsFactors = F)

dep_counts %<>% subset(.,RNAi_hit > 0)

master_table <- dplyr::select(t2,CDS_ID,CRISPR_class,RNAi_class)

dep_counts %<>% left_join(.,master_table,by="CDS_ID")

plot_names <- c("Pan-dependency","High-variance Dependency","Strongly Selective Dependency","Weakly Selective Dependency")

plot_list <- list()
for (pname in plot_names){
  sub_df <- subset(dep_counts, (CRISPR_class == pname) | (RNAi_class == pname))
  sub_df %<>% mutate(.,plot = pname)
  
  crispr_spec <- subset(dep_counts,(CRISPR_class == pname) & (RNAi_class != pname))$CDS_ID
  rnai_spec <- subset(dep_counts,(CRISPR_class != pname) & (RNAi_class == pname))$CDS_ID
  shared <- subset(dep_counts,(CRISPR_class == pname) & (RNAi_class == pname))$CDS_ID
  
  sub_df$group <- "none"
  sub_df$group[sub_df$CDS_ID %in% crispr_spec] <- "CRISPR"
  sub_df$group[sub_df$CDS_ID %in% rnai_spec] <- "RNAi"
  sub_df$group[sub_df$CDS_ID %in% shared] <- "shared"
  
  plot_list[[pname]] <- sub_df
}

plot_df <- bind_rows(plot_list)

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF","shared"="grey50")

plot_df$plot <- factor(plot_df$plot,levels=plot_names)

ggplot(subset(plot_df,plot %in% c("High-variance Dependency","Strongly Selective Dependency")),aes(x=RNAi_hit,y=CRISPR_hit,color=group)) +
  geom_point(shape=1,stroke=.85) +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_bw(base_size=11) +
  scale_color_manual(values=tech_pal) +
  facet_grid(plot ~ .) +
  xlab("RNAi Dep. Cell Line Count") +
  ylab("CRISPR Dep. Cell Line Count")
ggsave(file.path("figures","gene_deps_profiles_depCount_scatter_1x2_grid.pdf"),width=4,height=4)

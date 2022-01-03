
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

long_form_gs <- function(gene_list,mat_gs,pert_type){
  tmp_gs <- mat_gs[,colnames(mat_gs) %in% gene_list]
  tmp_gs_long <- melt(as.matrix(tmp_gs))
  tmp_gs_long$ID <- paste0(tmp_gs_long$Var1,"-",tmp_gs_long$Var2)
  tmp_gs_long %<>% rename(.,!!pert_type := value)
  return(tmp_gs_long)
}

tech <- c("CRISPR","RNAi")
se_metric <- c("Strongly Selective Dependency","High-variance Dependency")

all_plots <- list()
for (metric in se_metric){
    all_hits = subset(t2,(RNAi_class == metric) | (CRISPR_class == metric))$CDS_ID
    
    crispr_gs_long <- long_form_gs(all_hits,crispr_pr,"CRISPR")
    rnai_gs_long <- long_form_gs(all_hits,rnai_pr,"RNAi")
    plot_df <- left_join(crispr_gs_long,dplyr::select(rnai_gs_long,RNAi,ID),by="ID")
    
    shared_genes <- subset(t2,(RNAi_class == metric) & (CRISPR_class == metric))$CDS_ID
    plot_df$shared <- plot_df$Var2 %in% shared_genes
    
    distinct_CRISPR <- subset(t2,(RNAi_class != metric) & (CRISPR_class == metric))$CDS_ID
    plot_df$distinct_CRISPR <- plot_df$Var2 %in% distinct_CRISPR
    
    distinct_RNAi <- subset(t2,(RNAi_class == metric) & (CRISPR_class != metric))$CDS_ID
    plot_df$distinct_RNAi <- plot_df$Var2 %in% distinct_RNAi
    
    plot_df$group <- "None"
    plot_df$group[plot_df$shared] <- "shared"
    plot_df$group[plot_df$distinct_CRISPR] <- "CRISPR"
    plot_df$group[plot_df$distinct_RNAi] <- "RNAi"
    
    plot_df$metric <- metric
    all_plots[[metric]] <- plot_df
}

plot_df <- bind_rows(all_plots)

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF","shared"="grey50")

ggplot(plot_df,aes(x=RNAi,y=CRISPR,color=group)) +
  geom_point(alpha=.2,size=.1) +
  geom_hline(yintercept = .5,linetype="dashed") +
  geom_vline(xintercept = .5,linetype="dashed") +
  theme_bw(base_size=11) +
  scale_color_manual(values=tech_pal) +
  xlab("RNAi probability of dependency") +
  ylab("CRISPR probability of dependency") +
  facet_grid(metric ~ group) +
  theme(legend.position = "none")
ggsave(file.path("figures","gene_deps_profiles_probabilities_density2d_2x3_grid.pdf"),width=4,height=2.8)
ggsave(file.path("figures","gene_deps_profiles_probabilities_density2d_2x3_grid.png"),width=4,height=2.8)


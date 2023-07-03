
source("src/packages_paths.R")

# crispr_gs <- load.from.taiga(data.name='dependency-data-9cb5', data.version=2, data.file='crispr_gs')
# rnai_gs <- load.from.taiga(data.name='dependency-data-9cb5', data.version=2, data.file='rnai_gs')

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
rnai_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)

# exp <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=12, data.file='CCLE_depMap_18Q4_TPM_ProteinCoding')
exp <- fread(file.path("data/raw","depmap-omics-expression-rnaseq-tpm-18Q4.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
exp <- exp[rownames(crispr_gs),]

plot_df_crispr <- data.frame(CL=rownames(crispr_gs),
                             RAB6B_Exp=exp[,"RAB6B (51560)"],
                             RAB6A_GS=crispr_gs[,"RAB6A (5870)"],
                             dataset="CRISPR",
                             stringsAsFactors = F)
plot_df_rnai <- data.frame(CL=rownames(rnai_gs),
                           RAB6B_Exp=exp[,"RAB6B (51560)"],
                           RAB6A_GS=rnai_gs[,"RAB6A (5870)"],
                           dataset="RNAi",
                           stringsAsFactors = F)
plot_df <- rbind(plot_df_crispr,plot_df_rnai)
plot_df$RAB6B_Loss <- plot_df$RAB6B_Exp < 1.385

crispr_x <- subset(plot_df,(dataset == "CRISPR") & (RAB6B_Loss))$RAB6A_GS
crispr_y <- subset(plot_df,(dataset == "CRISPR") & (!RAB6B_Loss))$RAB6A_GS

rnai_x <- subset(plot_df,(dataset == "RNAi") & (RAB6B_Loss))$RAB6A_GS
rnai_y <- subset(plot_df,(dataset == "RNAi") & (!RAB6B_Loss))$RAB6A_GS

crispr_pval <- wilcox.test(crispr_x,crispr_y)$p.value
rnai_pval <- wilcox.test(rnai_x,rnai_y)$p.value

wes_pal <- c("TRUE"="black","FALSE"="#C7B19C")

ggplot(plot_df, aes(x=dataset, y=RAB6A_GS, color=RAB6B_Loss)) +
  geom_boxplot(fill=NA,outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = .15),size=.7,alpha=.7) +
  theme_classic(base_size=9) +
  scale_color_manual(values=wes_pal, name="RAB6B exp.") +
  xlab("") +
  ylab("RAB6A Gene Effect") + 
  theme(legend.position = "top") +
  theme(legend.text=element_text(size=7)) +
  theme(legend.title=element_text(size=9)) +
  annotate("text", x = .5, y = .5, hjust = 0,label =paste0("P-value",
                                                          "\nCRISPR: ",signif(crispr_pval,digits=2),
                                                          "\nRNAi: ",signif(rnai_pval,digits=2)),size=2)

ggsave(file.path("figures","predictive_markers_RAB6A_boxplot.pdf"),height=2,width=2.25)


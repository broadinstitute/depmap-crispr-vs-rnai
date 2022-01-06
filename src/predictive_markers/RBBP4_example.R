
source("src/packages_paths.R")

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
rnai_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)

exp <- fread(file.path("data/raw","depmap-omics-expression-rnaseq-tpm-18Q4.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
exp <- exp[rownames(crispr_gs),]

plot_df_crispr <- data.frame(CL=rownames(crispr_gs),
                             RBBP7_Exp=exp[,"RBBP7 (5931)"],
                             RBBP4_GS=crispr_gs[,"RBBP4 (5928)"],
                             dataset="CRISPR",
                             stringsAsFactors = F)
plot_df_rnai <- data.frame(CL=rownames(rnai_gs),
                           RBBP7_Exp=exp[,"RBBP7 (5931)"],
                           RBBP4_GS=rnai_gs[,"RBBP4 (5928)"],
                           dataset="RNAi",
                           stringsAsFactors = F)
plot_df <- rbind(plot_df_crispr,plot_df_rnai)

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

ggplot(plot_df, aes(x=RBBP7_Exp, y=RBBP4_GS, color=dataset)) +
  geom_density_2d() +
  geom_smooth(method="lm") +
  theme_classic(base_size=11) +
  scale_color_manual(values=tech_pal, name="") +
  xlab("RBBP7 Expression\nTPM (log2)") +
  ylab("RBBP4 Gene Effect") + 
  theme(legend.position = c(.2,.8))

ggsave(file.path("figures","predictive_markers_RBBP4_density2d_smooth.pdf"),height=2.2,width=2.25)





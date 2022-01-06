
source("src/packages_paths.R")

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
rnai_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)

cn <- fread(file.path("data/raw","depmap-omics-cn-gene-internal-18q4.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
cn <- cn[rownames(crispr_gs),]


plot_df_crispr <- data.frame(CL=rownames(crispr_gs),
                             MTAP_CN=cn[,"MTAP (4507)"],
                             PRMT5_GS=crispr_gs[,"PRMT5 (10419)"],
                             dataset="CRISPR",
                             stringsAsFactors = F)
plot_df_rnai <- data.frame(CL=rownames(rnai_gs),
                           MTAP_CN=cn[,"MTAP (4507)"],
                           PRMT5_GS=rnai_gs[,"PRMT5 (10419)"],
                           dataset="RNAi",
                           stringsAsFactors = F)
plot_df <- rbind(plot_df_crispr,plot_df_rnai)
plot_df$MTAP_Loss <- plot_df$MTAP_CN < -2

wes_pal <- c("TRUE"="black","FALSE"="#C7B19C")

ggplot(plot_df, aes(x=dataset, y=PRMT5_GS, color=MTAP_Loss)) +
  geom_point(position=position_jitterdodge(jitter.width = .15),size=.7,alpha=.7) +
  geom_boxplot(outlier.shape=NA,fill=NA) +
  theme_classic(base_size=11) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  scale_color_manual(values=wes_pal, name="MTAP CN Loss") +
  xlab("") +
  ylab("PRMT5 Gene Effect") + 
  theme(legend.margin=margin(t = -.5, unit='cm')) +
  theme(legend.position = "bottom")
ggsave(file.path("figures","predictive_markers_PRMT5_boxplot.pdf"),height=2,width=2.25)

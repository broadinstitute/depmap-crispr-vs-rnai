
# require(GSVA)

source("src/packages_paths.R")

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix()

#Exclude blood cell lines
master <- fread(file.path(data_raw,"depmap-cell-line-annotations.csv"))
master %<>% subset(.,`Primary Site` == "haematopoietic_and_lymphoid_tissue")

crispr_gs <- crispr_gs[!(rownames(crispr_gs) %in% master$DepMap_ID),]

exp <- fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-18Q4.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
exp <- exp[rownames(crispr_gs),]

ifn <- fread(file.path(data_raw,"gene-set-library.csv"))
ifn %<>% subset(.,gene_set == "interferon-stimulated gene expression")
ifn <- ifn$CDS_ID
ifn <- ifn[ifn %in% colnames(exp)]

geneSets <- list("ISG"=ifn)

ssgsea_res <- GSVA::gsva(t(exp),geneSets,method="ssgsea")

g1 <- "ADAR (103)"
cor(as.numeric(ssgsea_res),as.numeric(crispr_gs[,g1]),use="pairwise.complete")

plot_df <- data.frame(CL=rownames(crispr_gs),
                      CRISPR_GS=crispr_gs[,g1],
                      TPM=as.numeric(ssgsea_res),
                      stringsAsFactors = F)

binary_pal <- c("TRUE"="#925E9FFF","FALSE"="grey")

p <- ggplot(plot_df,aes(x=CRISPR_GS,y=TPM)) +
  geom_vline(xintercept = 0,linetype = "dashed",color="#625441",size=1) +
  geom_vline(xintercept = -1,linetype = "dashed",color="#C15327",size=1) +
  geom_smooth(method="lm") +
  geom_point(size=.7,alpha=.7,color="grey50") +
  scale_x_continuous(breaks=c(0,-.5,-1,-1.25)) +
  scale_color_manual(values=binary_pal) +
  xlab("ADAR CRISPR Gene Effect") +
  ylab("ISG Signature") +
  theme_classic(base_size=11) +
  theme(legend.position = c(0.8, 0.8))
p <- ggExtra::ggMarginal(p,type = "histogram")
ggsave(plot=p,file.path("figures","gene_deps_profiles_ADAR_IFN.pdf"),width=3,height=2.5)

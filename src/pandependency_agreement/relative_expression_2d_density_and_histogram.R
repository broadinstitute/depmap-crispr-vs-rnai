
# require(matrixStats)

source("src/packages_paths.R")

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc)

t2$pandep_group <- "init"
t2$pandep_group[t2$CRISPR_PD & (!t2$RNAi_PD)] <- "CRISPR-specific"
t2$pandep_group[t2$CRISPR_PD & t2$RNAi_PD] <- "shared"
t2 %<>% dplyr::select(.,CDS_ID,entrez_id,pandep_group,symbol)

rnai_pandep <- subset(t2,pandep_group %in% c("shared","RNAi-specific"))$CDS_ID

exp <- fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-18Q4.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
rnai <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
cls <- intersect(rownames(exp),rownames(rnai))
exp <- exp[cls,]
rnai <- rnai[cls,]

colnames(exp) <- entrez_to_cds(extract_entrez(colnames(exp)),hgnc)
colnames(rnai) <- entrez_to_cds(extract_entrez(colnames(rnai)),hgnc)
genes <- intersect(colnames(exp),colnames(rnai))
rnai_pandep <- intersect(rnai_pandep,genes)

exp <- exp[,rnai_pandep]
rnai <- rnai[,rnai_pandep]

rnai_cors <- list()
for (g in rnai_pandep){
  rnai_cors[[g]] <- cor(rnai[,g],exp[,g],use="pairwise.complete")
}

basic_stats <- data.frame(Gene=colnames(exp),exp_mean=colMeans(exp,na.rm=T),exp_var=matrixStats::colVars(exp,na.rm=T),pearson=unlist(rnai_cors))

example_gene <- "PSMC3 (5702)"

ggplot(basic_stats,aes(x=pearson)) +
  geom_histogram(bins=10,color="grey50",fill="white") +
  geom_vline(xintercept = basic_stats[example_gene,"pearson"],color="blue") +
  theme_classic(base_size=11) +
  xlab("Pearson Correlation") +
  ylab("Gene Count")
ggsave(file.path("figures","pandependency_agreement_RNAi_pandep_exp_rnai_cor_hist.pdf"),height=2.5,width=2.7)

highlight <- data.frame(RNAi=rnai[,example_gene],Expression=exp[,example_gene])
ggplot(highlight,aes(x=Expression,y=RNAi)) +
  geom_point(color="darkblue",alpha=.5,size=.3) +
  geom_smooth(method="lm") +
  theme_bw(base_size=6) +
  xlab("Expression TPM (log2)") +
  ylab("RNAi Gene Effect")
ggsave(file.path("figures","pandependency_agreement_PSMC3_scatter.pdf"),height=1,width=1)

cor(highlight$RNAi,highlight$Expression,use="pairwise.complete")

exp <- scale(exp,center=T,scale=T)
rnai <- scale(rnai,center=T,scale=T)

long_exp <- melt(exp)
rownames(long_exp) <- paste0(long_exp$Var1,"_",long_exp$Var2)
long_exp %<>% dplyr::select(.,value) %>% dplyr::rename(.,Expression=value) %>% rownames_to_column(.,var="ID")
long_rnai <- melt(rnai)
rownames(long_rnai) <- paste0(long_rnai$Var1,"_",long_rnai$Var2)
long_rnai %<>% dplyr::select(.,value) %>% dplyr::rename(.,RNAi=value) %>% rownames_to_column(.,var="ID")

plot_df <- left_join(long_exp,long_rnai,by="ID")

pal <- c("#899DA4","#C93312","#DC863B","#FAEFD1")

ggplot(plot_df) +
  geom_hex(aes(y=Expression,x=RNAi,colour = ..count..),bins=100) +
  scale_color_gradientn(colors=pal,name = "count", trans = "log",
                        breaks = 10^(0:6)) +
  scale_fill_gradientn(colors=pal,name = "count", trans = "log",
                       breaks = 10^(0:6)) +
  xlim(c(-8,8)) +
  ylim(c(-8,8)) +
  xlab("Relative RNAseq Expression") +
  ylab("Relative RNAi Gene Effect") +
  geom_smooth(aes(y=Expression,x=RNAi),method="lm") +
  theme_classic(base_size=11) 
  # theme(legend.position="none")
ggsave(file.path("figures","pandependency_agreement_RNAi_pandep_relative_expression.pdf"),height=2.6,width=3.75)

cor(plot_df$RNAi,plot_df$Expression,use="pairwise.complete")

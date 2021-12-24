
source("src/packages_paths.R")

#Get expression, CRISPR, and RNAi data
dataset <- list("Expression"=fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-18Q4.csv")),
                "CRISPR"=fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")),
                "RNAi" = fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")))
dataset <- lapply(dataset,function(x){x %<>% column_to_rownames(.,var="V1") %>% as.matrix(.)})

cls <- lapply(dataset,rownames)
cls <- Reduce(intersect,cls)

genes <- lapply(dataset,colnames)
genes <- Reduce(intersect,genes)

dataset <- lapply(dataset,function(d){d[cls,genes]})

#Calculate gene means
gene_means <- lapply(dataset,function(d){colMeans(d,na.rm=T)})
gene_means <- bind_cols(gene_means)
gene_means %<>% add_column(.,Gene=genes,.before=1)

#plot mean expression vs mean dependency for RNAi and CRISPR
gm1 <- mutate(gene_means,Dependency=RNAi,tech="RNAi") %>% dplyr::select(.,Gene,Expression,Dependency,tech)
gm2 <- mutate(gene_means,Dependency=CRISPR,tech="CRISPR") %>% dplyr::select(.,Gene,Expression,Dependency,tech)
gene_means_df <- rbind(gm1,gm2)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1 %<>% subset(.,high_confidence)
t1$CDS_ID <- entrez_to_cds(t1$entrez_id,hgnc)
gene_means_libconf <- subset(gene_means_df, Gene %in% t1$CDS_ID)

mypal <- c("CRISPR"="#DC0000FF","RNAi"="#8491B4FF")

gene_means_df$grouping <- ">0"
gene_means_df$grouping[gene_means_df$Dependency <= 0] <- "[-0.5,0]"
gene_means_df$grouping[gene_means_df$Dependency < -.5] <- "[-1,-0.5]"
gene_means_df$grouping[gene_means_df$Dependency < -1] <- "<-1"

gene_means_df$grouping <- factor(gene_means_df$grouping,levels=c("<-1","[-1,-0.5]","[-0.5,0]",">0"))
ggplot(gene_means_df,aes(x=grouping,y=Expression,color=tech)) +
  geom_boxplot() +
  theme_classic(base_size=11) +
  scale_color_manual(values=mypal,name="") +
  theme(legend.position = "top") + 
  xlab("Mean Gene Effect") +
  ylab("Mean Expression") +
  theme(legend.spacing = unit(0, "null")) +
  theme(plot.margin = unit(c(-.03, .01, .01, .01), "null"))
ggsave(file.path("figures","pandependency_agreement_dependency_vs_expression_boxplot.pdf"),height=2.85,width=2.3)

table(subset(gene_means_df,tech=="CRISPR")$grouping)
table(subset(gene_means_df,tech=="RNAi")$grouping)

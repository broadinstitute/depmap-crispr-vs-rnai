library(tidyverse)
library(magrittr)
library(data.table)
library(ggrepel)

se_genes <- c("ASH2L", "FUBP1", "RNF125")

#### Gene effect correlation

lib_gene_effect <- list("Avana"="data/raw/gene-effect-scaled-crispr-avana.csv",
                        "KY"="data/raw/gene-effect-scaled-crispr-ky.csv",
                        "Achilles"="data/raw/gene-effect-scaled-rnai-achilles.csv",
                        "DRIVE"="data/raw/gene-effect-scaled-rnai-drive.csv")

lib_gene_effect <- lapply(lib_gene_effect,function(x){fread(x) %>% column_to_rownames(.,var="V1")})
lib_gene_effect <- lapply(lib_gene_effect,function(x){colnames(x) <- gsub(" .*","",colnames(x));return(x)})

cls <- lapply(lib_gene_effect,function(x){rownames(x)})
cls <- Reduce(intersect,cls)

lib_gene_effect <- lapply(lib_gene_effect, function(x){x[cls,]})

res_df <- list()
for (g in se_genes){
  
  crispr_r <- cor(lib_gene_effect[["Avana"]][,g],lib_gene_effect[["KY"]][,g],use="pairwise.complete")
  rnai_r <- cor(lib_gene_effect[["Achilles"]][,g],lib_gene_effect[["DRIVE"]][,g],use="pairwise.complete")
  
  res_df[[g]] <- data.frame(gene=g,CRISPR=crispr_r,RNAi=rnai_r)
  
}
plot_df <- bind_rows(res_df)

ggplot(plot_df,aes(x=CRISPR,y=RNAi)) +
  geom_point() +
  geom_abline(slope=1,linetype="dashed") +
  xlab("CRISPR correlation (Avana,KY)") +
  ylab("RNAi correlation (Achilles,DRIVE)") +
  geom_label_repel(data=plot_df,
                   aes(x=CRISPR,y=RNAi,label=gene),
                   min.segment.length = 0) +
  theme_bw(base_size=16) 
ggsave("figures/efficacy_specificity_ASH2L_FUBP1_RNF125_correlation.pdf",
              height=4,width=4)

#### Predictive accuracy
se_genes <- c("ASH2L (9070)", "FUBP1 (8880)", "RNF125 (54941)")

crispr_pred <- fread("data/processed/ensemble-prediction-summary-crispr-matched.csv")
crispr_pred %<>% subset(.,gene %in% se_genes) %>% 
  subset(.,best) %>%
  dplyr::select(.,gene,CRISPR=pearson,CRISPR_feature=feature0)
crispr_pred$CRISPR_feature_type <- gsub(".+_","",crispr_pred$CRISPR_feature)
crispr_pred$CRISPR_feature <- gsub(" .*","",crispr_pred$CRISPR_feature)

rnai_pred <- fread("data/processed/ensemble-prediction-summary-rnai-matched.csv")
rnai_pred %<>% subset(.,gene %in% se_genes) %>% 
  subset(.,best) %>%
  dplyr::select(.,gene,RNAi=pearson,RNAi_feature=feature0)
rnai_pred$RNAi_feature_type <- gsub(".+_","",rnai_pred$RNAi_feature)
rnai_pred$RNAi_feature <- gsub(" .*","",rnai_pred$RNAi_feature)

pred <- full_join(crispr_pred,rnai_pred,by="gene")
pred$target <- gsub(" .*","",pred$gene)

ggplot(pred,aes(x=CRISPR,y=RNAi)) +
  geom_point() +
  geom_abline(slope=1,linetype="dashed") +
  xlab("CRISPR predictive accuracy") +
  ylab("RNAi predictive accuracy") +
  geom_label_repel(data=pred,
                   aes(x=CRISPR,y=RNAi,label=target),
                   min.segment.length = 0) +
  theme_bw(base_size=16) 
ggsave("figures/efficacy_specificity_ASH2L_FUBP1_RNF125_predictability.pdf",
       height=4,width=4)

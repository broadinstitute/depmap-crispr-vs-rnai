#setwd("/Users/mburger/tda-pipeline/pipeline")

library(purrr)
library(data.table)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
feat_file <- args[1]
pred_file <- args[2]
obs_file <- args[3]
outfile <- args[4]

#Load files
pr <- as.data.frame(fread(obs_file)) %>% column_to_rownames(.,var="Row.name")
feats <- read_csv(feat_file)
preds <- as.data.frame(read_csv(pred_file)) %>% column_to_rownames(.,var="Row.name")
preds <- preds[rownames(pr),]

print(head(pr))
print(head(preds))

#Calculate correlation between obs and pred
r <- mapply(cor,preds,pr[,colnames(preds)],use="complete")

feats <- feats[match(names(r),feats$gene),]
feats$pearson <- r

#Write summary file
performance <- feats[,c("gene","model","pearson")]
feat_imp <- feats[,grepl("feature",colnames(feats))]
result_df <- cbind(performance,feat_imp)

write.csv(result_df,file=outfile,row.names=F)

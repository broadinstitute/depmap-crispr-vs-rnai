library(purrr)
library(data.table)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
model_summary_list <- args[2:length(args)]

results_list <- list()
for (f_name in model_summary_list){
  feats <- read_csv(f_name)
  results_list[[f_name]] <- feats
}

result_df <- do.call(rbind,results_list)
result_df <- result_df[order(result_df$pearson,decreasing = T),]
result_df %<>% add_column(.,best=!duplicated(result_df$gene),.after="pearson")

write.csv(result_df,file=outfile,row.names=F)


library(purrr)
library(data.table)
library(tidyverse)
library(magrittr)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
sparkle_params <- args[1]
tmp_file_dir <- args[2]
out_file <- args[3]

param_table <- read.csv(sparkle_params,stringsAsFactors = F)
param_table$suffix <- paste0("-",param_table$start,"-",param_table$end)
param_table$prefix <- file_path_sans_ext(out_file)
in_files <- paste0(param_table$prefix,param_table$suffix,".",file_ext(out_file))

LRT_list <- list()
for (i in 1:length(in_files)){
  LRT_list[[as.character(i)]] <- fread(file.path(tmp_file_dir,in_files[i]))
}
LRT_result <- bind_rows(LRT_list)
write_csv(LRT_result,path=out_file)
  

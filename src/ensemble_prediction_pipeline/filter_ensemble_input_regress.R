
library(jsonlite)
library(matrixStats)
library(data.table)
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
model_config <- args[1]
task_filter <- args[2]
gs_dataset <- args[3]
var_thresh <- as.numeric(args[4])
outfile_name <- args[5]
test <- args[6]
feature_index <- args[7]

model_def <- read_json(model_config)

# Set input to gene effect
ensemble_input <- data.frame(fread(gs_dataset),stringsAsFactors = F,check.names = F) %>% column_to_rownames(.,var="Row.name")
print(paste0("Initial dimensions of input: ",paste(dim(ensemble_input),collapse=",")))

# Filter cell lines for inclusion in required feature sets
required_feats <- lapply(model_def,function(x){x[["Required"]]})
required_feats <- unique(unlist(lapply(required_feats,function(x){unlist(x)})))
if (sum(!(required_feats %in% c("none","Confounders"))) > 0){
  dataset_index <- fread(feature_index)
  dataset_index %<>% subset(.,dataset %in% required_feats)
  dataset_list <- lapply(dataset_index$filename,function(x){fread(x)})
  cls_list <- lapply(dataset_list,function(x){x$`Row.name`})
  cls_keep <- Reduce(intersect,cls_list)
  cls_keep <- intersect(cls_keep,rownames(ensemble_input))
  ensemble_input <- ensemble_input[cls_keep,]
}
print(paste0("Required feature filter of input: ",paste(dim(ensemble_input),collapse=",")))

# Filter ensemble input for genes with deps > folds (count), variance > x (variance), or skip (none)
if (task_filter == "variance"){

  dep_var <- colVars(as.matrix(ensemble_input),na.rm=T)
  dep_rank <- frankv(dep_var,order=-1,ties.method="max")
  names(dep_rank) <- colnames(ensemble_input)
  keep_genes <- names(dep_rank)[dep_rank <= var_thresh]
  ensemble_input <- ensemble_input[,keep_genes]

} else if (task_filter == "none") {
  print("Skipping gene filtering of ensemble input")
} else {
  stop(paste0("Error: Unsupported task filter ",task_filter))
}
print(paste0("Gene filter of input: ",paste(dim(ensemble_input),collapse=",")))

###########################
####### TEST SET ##########
###########################
if (test != "NA"){
  print("Filtering for test set only")
  test_genes <- fread(test,sep="\t") %>% dplyr::select(.,Gene)
  test_genes <- as.character(test_genes$Gene)
  ensemble_input <- ensemble_input[,colnames(ensemble_input) %in% test_genes]
}
print(paste0("Test set gene filter of input: ",paste(dim(ensemble_input),collapse=",")))

# Ouput file
ensemble_input %<>% as.data.frame(.) %>% rownames_to_column(.,var="Row.name")
write_csv(ensemble_input,outfile_name)
if ( !((ncol(ensemble_input) >= 3) & (nrow(ensemble_input) >= 50))){
  stop(paste0("Not enough samples to run ensemble on ",group_name," subset"))
}

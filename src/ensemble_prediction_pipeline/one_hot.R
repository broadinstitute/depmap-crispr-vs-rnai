
# context_table <- "/Users/mburger/tda-pipeline/pipeline/inputs/cell-line-lineage.csv"
# rds_file <- "/Users/mburger/tda-pipeline/pipeline/inputs/cell-line-lineage-one-hot.rds"
# text_file <- "/Users/mburger/tda-pipeline/pipeline/inputs/cell-line-lineage-one-hot.csv"

# args <- commandArgs(trailingOnly = TRUE)
# context_table <- args[1]
# rds_file <- args[2]
# text_file <- args[3]

# context_table <- data.frame(fread(context_table),stringsAsFactors = F) %>% column_to_rownames(.,var="Row.name")

get_dummies <- function(context_table){
  
  mat_list <- list()
  for (v in colnames(context_table)){
    all_columns <- unique(context_table[,v])
    all_columns <- all_columns[!is.na(all_columns)]
    all_columns <- all_columns[all_columns != ""]
    all_columns <- gsub(" ", "-", all_columns)
    all_columns <- gsub("/", "-", all_columns)
    all_columns <- gsub(":", "-", all_columns)
    drop_chars <- c(",",";")
    all_columns <- gsub(pattern = paste(drop_chars, collapse = "|"), replacement = "", all_columns)
    
    new_mat <- matrix(rep(F,length(all_columns)*nrow(context_table)),
                      ncol=length(all_columns),
                      nrow=nrow(context_table),
                      dimnames = list(rownames(context_table),all_columns))
    
    for (v2 in all_columns){
      new_mat[,v2] <- context_table[,v] == v2
    }
    
    mat_list[[v]] <- new_mat + 0
    
  }
  
  #saveRDS(mat_list,file=rds_file)
  
  for (dataset in names(mat_list)){
    mat <- mat_list[[dataset]]
    colnames(mat) <- paste0(colnames(mat),"_",dataset)
    mat_list[[dataset]] <- mat
  }

  text_mat <- do.call(cbind,mat_list)
  text_mat <- as.data.frame(text_mat)
  return(text_mat)
  
}
  
get_dummies_grepl <- function(sample_info,groups){
  
  sample_info[,2] <- tolower(sample_info[,2])
  
  onehot_mat <- matrix(rep(F,nrow(sample_info)*length(groups)),nrow=nrow(sample_info),ncol=length(groups))
  colnames(onehot_mat) <- groups
  rownames(onehot_mat) <- sample_info[,1]
  for (g in groups){
    
    search_g <- tolower(g)
    onehot_mat[,g] <- grepl(search_g,sample_info[,2])
    
  }
  
  return(onehot_mat+0)
  
}
  

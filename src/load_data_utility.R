
library(httr)

# Get mapping of local directories to figshare from config json
figshare_mapping <- list("data/raw"="https://figshare.com/",
                         "data/processed"="https://figshare.com/",
                         "tables"="https://figshare.com/")

load_local <- function(file_path,data_type){
  
  if (data_type == "table"){
    return(fread(file_path))
  } else if (data_type == "matrix"){
    tmp_file <- fread(file_path)
    row_col <- colnames(tmp_file)[1]
    tmp_file %<>% column_to_rownames(.,var=row_col)
    tmp_file %<>% as.matrix(.)
    return(tmp_file)
  } else {
    stop("Unknown data type")
  }
  
}

get_figshare_file <- function(local_dir,filename){
  
  base_url <- figshare_mapping[[local_dir]]
  url <- paste0(base_url,"/",filename)
  GET(url, progress(), write_disk(file.path(local_dir,filename), overwrite=TRUE))
  
}

load_data <- function(local_dir,filename,data_type,figshare_mapping_dict=figshare_mapping){
  
  local_file <- file.path(local_dir,filename)
  
  #Download file from Figshare if it doesn't exist locally
  # if (!file.exists(local_file)){
  #   
  #   get_figshare_file(local_dir,filename)
  #   
  # } 
  
  return(load_local(file_path=local_file,data_type=data_type))
}



library(tidyverse)
library(magrittr)
library(data.table)

source("src/R/id_utility.R")

#' Load predictive feature datasets
#' 
#' Retrieves datasets from Taiga and concatenates them according to the specified operation
#'
#' @param feature_datasets vector containing the names of the predictive feature datasets, all data must have samples as rows
#' @param target_data matrix of target data, samples as rows 
#' @param required_datasets takes intersection of cell lines in required datasets and union of others
#' @param file_key list where names are predictive feature dataset names and values are taiga data.name, data.file, norm
#'
#' @return concat predictive matrix of all features
#' @export
load.features <- function(feature_datasets,required_datasets=NA,file_key,shape="wide"){
  
  stopifnot(all(feature_datasets %in% names(file_key)))
  hgnc <- fread(file_key[["HGNC"]]) %>% select(.,symbol,entrez_id)
  hgnc$entrez_id %<>% as.character()
  hgnc %<>% subset(.,!is.na(entrez_id) & (entrez_id != ""))
  
  X <- list()
  for (fd in feature_datasets){
    #data <- load.from.taiga.retry(data.name=file_key[[fd]]$data.name,data.file=file_key[[fd]]$data.file,no.save=T)
    #colnames(data) <- gsub(" .*","",colnames(data))
    data <- fread(file_key[[fd]]) %>% as.data.frame(.) %>% column_to_rownames(.,var="Row.name") %>% 
      set_colnames(.,extract_entrez(colnames(.)))
    
    data <- data[,colnames(data) %in% hgnc$entrez_id]
    
    #check if binary to determine if norm is ind_zscore or none
    d_check <- data[,1:10]
    d_check[is.na(d_check)] <- 0
    binary <- all(unique(as.numeric(unlist(d_check))) %in% c(0,1)) 
    
    if (binary){
      X[[fd]] <- list(data = data,norm = "none")
    } else {
      X[[fd]] <- list(data = data,norm = "ind_zscore")
    }
    
  }
  
  if (length(required_datasets) == 1 && is.na(required_datasets)){
    used_CLs <- lapply(X,function(x){rownames(x$data)})
    used_CLs <- Reduce(union,used_CLs)
  } else {
    used_CLs <- lapply(required_datasets,function(x){rownames(X[[x]]$data)})
    used_CLs <- Reduce(intersect,used_CLs)
  }
  
  X_norm <- list()
  for (x_name in names(X)){
    x_norm <- X[[x_name]]$norm
    x_data <- X[[x_name]]$data
    
    #gene_variance <- colVars(x_data,rows=rownames(x_data) %in% used_CLs,na.rm=T)
    #keep_genes <- subset(colnames(x_data),gene_variance > .01)
    #x_data <- x_data[,keep_genes]
    
    if (identical(x_norm,"ind_zscore")) {
      norm_fac <- pmax(1e-10, apply(x_data, 2, sd, na.rm = T))
      x_data <- scale(x_data, center = TRUE, scale = norm_fac)
    }
    
    if (identical(shape,"wide")){
      colnames(x_data) <- paste0(colnames(x_data),"_",x_name)
    }
    
    x_data <- x_data[rownames(x_data) %in% used_CLs, ,drop=F]
    missing_cls <- used_CLs[!(used_CLs %in% rownames(x_data))]
    if (length(missing_cls) > 0){
      missing_mat <- matrix(rep(NA,length(missing_cls)*ncol(x_data)),
                            nrow=length(missing_cls),
                            ncol=ncol(x_data),
                            dimnames = list(missing_cls,colnames(x_data)))
      x_data <- rbind(x_data,missing_mat)
    }
    
    if (identical(shape,"wide")){
      X_norm[[x_name]] <- x_data[used_CLs, ,drop=F] %>% as.data.frame(.) %>% rownames_to_column(.,var="Row.name")
    } else if (identical(shape,"combined")){
      X_norm[[x_name]] <- x_data[used_CLs, ,drop=F]
    }
  }
  
  #Make all matrices match columns if combining 
  if (identical(shape,"combined")){
    used_genes <- lapply(X_norm,function(x){colnames(x)})
    used_genes <- Reduce(union,used_genes)
    
    for (x_name in names(X_norm)){
      
      x_data <- X_norm[[x_name]]
      missing_genes <- used_genes[!(used_genes %in% colnames(x_data))]
      if (length(missing_genes) > 0){
        missing_mat <- matrix(rep(NA,length(missing_genes)*nrow(x_data)),
                              nrow=nrow(x_data),
                              ncol=length(missing_genes),
                              dimnames = list(rownames(x_data),missing_genes))
        x_data <- cbind(x_data,missing_mat)
        x_data <- x_data[,used_genes,drop=F]
      }
      
      X_norm[[x_name]] <- x_data
    }
    
  } else if (identical(shape,"wide")){
    X_norm %<>% reduce(full_join, by="Row.name")
  }
  
  return(X_norm)
}

#' Filter X predictive matrix for features that contain the target gene name
#' 
#' Current required data types for self filtering are MUThs and GE. Cell lines will be dropped without these.
#' 
#' @param gene target gene to use as search term
#' @param X_all matrix of all features
#' 
#' @return filtered X matrix for variables that have a basename that matches the target gene
filter.self <- function(gene,X_all,required_features=c("GE")){
  
  X_basename <- gsub("_.*","",colnames(X_all))
  keep <- grep(gene,X_basename,fixed=TRUE)
  X <- X_all[,keep,drop=F]
  
  if (!is.na(required_features[1])){
    X <- check.requirements(X,required_features)
  }
  
  return(X)
}

#' Check dataset requirements
#'
#' Filter cell lines based on which ones have the required features and check that they are not all NA. 
#' 
#' @param X predictive matrix (cell lines by features) with feature types as suffix of feature name
#' @param required_features features that must be included 
#' 
#' @return X filtered for cell lines with required features or NULL if all lines fail
check.requirements <- function(X,required_features){
  
  X_vtypes <- gsub(".*_","",colnames(X))
  for (req in required_features){
    if(sum(X_vtypes %in% req) > 0 && !is.null(X)){
      X <- X[rowSums(!is.na(X[,X_vtypes %in% req,drop=F])) > 0,,drop=F]
    } else {
      X <- NULL
    }
  }
  
  return(X)
  
}
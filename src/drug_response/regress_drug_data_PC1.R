
# require(WGCNA)

source("src/packages_paths.R")

master_file <- fread(file.path(data_raw,"depmap-cell-line-annotations-v846.csv"))
depmap_lines <- master_file$DepMap_ID

dose_datasets <- list("ctd2"="drug-screen-viability-filtered-ctd2.csv",
                      "gdsc"="drug-screen-viability-filtered-gdsc.csv",
                      "prism"="drug-screen-viability-filtered-prism.csv")
dose_datasets <- lapply(dose_datasets,function(x){fread(file.path(data_processed,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

for (dset in names(dose_datasets)){
  tmp.dose <- dose_datasets[[dset]]
  
  #Filter for cell lines included in DepMap
  common_cls <- intersect(depmap_lines, rownames(tmp.dose))
  tmp.dose <- tmp.dose[common_cls,]
  
  #Impute missing values using column median
  indx <- which(is.na(tmp.dose), arr.ind = TRUE)
  tmp.filled <- tmp.dose
  tmp.filled[indx] <- matrixStats::colMedians(tmp.filled, na.rm = TRUE)[indx[, 2]]
  
  #Remove PC1
  tmp.noPC1 <- WGCNA::removePrincipalComponents(tmp.filled, 1)
  
  #Replace imputed values with original NA
  tmp.noPC1[is.na(tmp.dose)] <- NA
  
  tmp.noPC1 %<>% as.data.frame() %>% rownames_to_column('Row.name')
  write_csv(tmp.noPC1,file.path(data_processed,paste0("drug-screen-viability-filtered-noPC1-",dset,".csv")))
}


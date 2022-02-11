
source("src/packages_paths.R")

gene_cor_recall_rank <- function(d1,d2,d1_name,d2_name,hgnc){
  
  #Align cell lines as rows
  cls <- intersect(rownames(d1),rownames(d2))
  d1 <- d1[cls,]
  d2 <- d2[cls,]
  d1 <- d1[,!grepl("&",colnames(d1))]
  d2 <- d2[,!grepl("&",colnames(d2))]
  
  #Align genes as columns (switch to stable entrez ids)
  colnames(d1) <- extract_entrez(colnames(d1))
  colnames(d2) <- extract_entrez(colnames(d2))
  all_ids <- intersect(colnames(d1),colnames(d2))
  d2 <- d2[,all_ids]
  d1 <- d1[,all_ids]
  
  #Calculate pearson
  pearson <- cor(d1,d2,use="pairwise.complete")
  
  #Add pearson coefficient to results (diagonal is gene-gene)
  result_table <- data.frame(entrez_id=all_ids,r=diag(pearson),stringsAsFactors = F)
  
  #rank cors and find match
  #For each gene in D2 correlate to all genes in D1 and determine how high the gene itself ranks
  rank_tmp <- pearson
  for (i in 1:ncol(rank_tmp)){
    rank_tmp[,i] <- frankv(rank_tmp[,i],order=-1)
  }
  rank_tmp <- diag(rank_tmp)
  result_table %<>% add_column(!!(paste0("rank_in_",d1_name)):=rank_tmp)
  
  #For each gene in D1 correlate to all genes in D2 and determine how high the gene itself ranks
  rank_tmp <- pearson
  for (i in 1:ncol(rank_tmp)){
    rank_tmp[i,] <- frankv(rank_tmp[i,],order=-1)
  }
  rank_tmp <- diag(rank_tmp)
  result_table %<>% add_column(!!(paste0("rank_in_",d2_name)):=rank_tmp)
  
  #annotate entrez with HUGO symbol
  result_table %<>% add_column(.,symbol=entrez_to_symbol(result_table$entrez_id,hgnc),.before=1)
  
  return(result_table)
  
}

# crispr_ky_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-ky.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
# crispr_avana_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-avana.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
# 
# rnai_achilles_gs <- fread(file.path(data_processed,"gene-effect-scaled-rnai-achilles.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
# rnai_drive_gs <- fread(file.path(data_processed,"gene-effect-scaled-rnai-drive.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)

gs_list_crispr <- list("crispr_ky_gs" = "gene-effect-scaled-crispr-ky.csv",
                       "crispr_avana_gs" = "gene-effect-scaled-crispr-avana.csv")
gs_list_crispr <- lapply(gs_list_crispr,function(x){load_data(local_dir=data_raw,filename=x,data_type="matrix")})

gs_list_rnai <- list("rnai_achilles_gs" = "gene-effect-scaled-rnai-achilles.csv",
                     "rnai_drive_gs" = "gene-effect-scaled-rnai-drive.csv")
gs_list_rnai <- lapply(gs_list_rnai,function(x){load_data(local_dir=data_processed,filename=x,data_type="matrix")})
gs_list <- c(gs_list_crispr,gs_list_rnai)

# hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc <- load_data(local_dir=data_raw,filename="hgnc-complete-set.csv",data_type="table")
hgnc$entrez_id %<>% as.character()

comparisons <- list(
  "CRISPR-CRISPR"=list("d1_name"="CRISPR-Avana",
                       "d2_name"="CRISPR-KY",
                       "d1"=gs_list[["crispr_avana_gs"]],
                       "d2"=gs_list[["crispr_ky_gs"]],
                       "hgnc"=hgnc),
  "RNAi-RNAi"=list("d1_name"="RNAi-Achilles",
                   "d2_name"="RNAi-DRIVE",
                   "d1"=gs_list[["rnai_achilles_gs"]],
                   "d2"=gs_list[["rnai_drive_gs"]],
                   "hgnc"=hgnc))

for (comp in names(comparisons)){
  comp_params <- comparisons[[comp]]
  
  result_table <- gene_cor_recall_rank(d1=comp_params[["d1"]],
                                       d2=comp_params[["d2"]],
                                       d1_name=comp_params[["d1_name"]],
                                       d2_name=comp_params[["d2_name"]],
                                       hgnc=comp_params[["hgnc"]])
  
  outfile <- file.path(data_processed,paste0(comp_params[["d1_name"]],"_vs_",comp_params[["d2_name"]],"_cor_recall_vals.csv")) 
  write_csv(result_table,file=outfile)
  
}




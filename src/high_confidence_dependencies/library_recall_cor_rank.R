
library(tidyverse)
library(magrittr)
library(data.table)

source("/Users/mburger/dynamic-duo/src/R/id_utility.R")

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

# gene_scores <- list(all=list(crispr_avana = load.from.taiga(data.name='avana-public-tentative-19q1-6956', data.version=2, data.file='gene_effect_corrected'),
#                              crispr_ky = load.from.taiga(data.name='ceres-8a62', data.version=1, data.file='ceres_ky_scaled'),
#                              rnai_achilles = load.from.taiga(data.name='demeter2-achilles-5386', data.version=13, data.file='gene_effect'),
#                              rnai_drive = load.from.taiga(data.name='demeter2-drive-0591', data.version=12, data.file='gene_effect')),
#                     overlap=list(ceres=load.from.taiga(data.name='dependency-data-9cb5', data.version=2, data.file='crispr_gs'),
#                                  d2=load.from.taiga(data.name='dependency-data-9cb5', data.version=2, data.file='rnai_gs')))

# crispr_broad_gs <- load.from.taiga(data.name='dependency-data-9cb5', data.version=2, data.file='crispr_gs')
# rnai_combined_gs <- load.from.taiga(data.name='dependency-data-9cb5', data.version=2, data.file='rnai_gs')
crispr_ky_gs <- load.from.taiga(data.name='inputs-ef6d', data.version=1, data.file='ceres_ky_scaled')
crispr_avana_gs <- load.from.taiga(data.name='inputs-ef6d', data.version=1, data.file='gene_effect_corrected')
rnai_achilles_gs <- load.from.taiga(data.name='inputs-7713', data.version=1, data.file='DEMETER2_Achilles_gene_effect')
rnai_drive_gs <- load.from.taiga(data.name='inputs-7713', data.version=1, data.file='DEMETER2_DRIVE_gene_effect')

hgnc <- load.from.taiga(data.name='hgnc-6825', data.version=2, data.file='hgnc_complete_set_090318')
hgnc$entrez_id %<>% as.character()

comparisons <- list(
  # "CRISPR-RNAi"=list("d1_name"="CRISPR-DepMap",
  #                    "d2_name"="RNAi-DepMap",
  #                    "d1"=crispr_broad_gs,
  #                    "d2"=rnai_combined_gs,
  #                    "hgnc"=hgnc),
  "CRISPR-CRISPR"=list("d1_name"="CRISPR-Avana",
                       "d2_name"="CRISPR-KY",
                       "d1"=crispr_avana_gs,
                       "d2"=crispr_ky_gs,
                       "hgnc"=hgnc),
  "RNAi-RNAi"=list("d1_name"="RNAi-Achilles",
                   "d2_name"="RNAi-DRIVE",
                   "d1"=rnai_achilles_gs,
                   "d2"=rnai_drive_gs,
                   "hgnc"=hgnc))

for (comp in names(comparisons)){
  comp_params <- comparisons[[comp]]
  
  result_table <- gene_cor_recall_rank(d1=comp_params[["d1"]],
                                       d2=comp_params[["d2"]],
                                       d1_name=comp_params[["d1_name"]],
                                       d2_name=comp_params[["d2_name"]],
                                       hgnc=comp_params[["hgnc"]])
  
  outfile <- file.path("/Users/mburger/dynamic-duo-biorxiv/figures/library_agreement",paste0(comp_params[["d1_name"]],"_vs_",comp_params[["d2_name"]],"_cor_recall_vals.csv")) 
  write_csv(result_table,path=outfile)
  
}




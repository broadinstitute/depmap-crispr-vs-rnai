
source("src/packages_paths.R")

#Construct related list by adding KEGG to PPI,CORUM,paralog
related <- fread(file.path(data_raw,"gene-set-related-features.csv"))
gs_df <- fread(file.path(data_raw,"gene-set-kegg.csv"))
related <- bind_rows(related,gs_df)

#Enrichment 

# master <- load.from.taiga(data.name='master-table-na-matched-6315', data.version=1, data.file='master-table')

cor_data <- list(CRISPR="codependency_CRISPR_pearson_baseline.rds",
                 RNAi="codependency_RNAi_pearson_baseline.rds",
                 SNF="codependency-CRISPR-RNAi-SNF.rds")
cor_data <- lapply(cor_data,function(x){readRDS(file.path(data_processed,x))})

#Verify that all genes have at least 3 dependent cell lines and matching rownames and colnames
t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,(`CRISPR_numDeps(50-100)` >= 3) | (`RNAi_numDeps(50-100)` >= 3))

cor_cols <- colnames(cor_data[["CRISPR"]])
stopifnot(all(extract_entrez(cor_cols) %in% t2$entrez_id))

cor_data <- lapply(cor_data,function(x){x[cor_cols,cor_cols]})

#Descriptive stats of each similarity matrix
cor_means <- lapply(cor_data,function(x){mean(x,na.rm=T)})
cor_sd <- lapply(cor_data,function(x){sd(x,na.rm=T)})

tmp_related <- subset(related,target %in% cor_cols)
tmp_related %<>% subset(.,partner %in% cor_cols)
tmp_related %<>% subset(.,target %in% cor_cols)
tmp_related %<>% dplyr::select(.,target,partner)
tmp_related %<>% unique(.)

pvals <- list()
for (gene_target in cor_cols){
  
  in_set <- subset(tmp_related,target == gene_target)$partner
  in_set <- intersect(cor_cols,in_set)
  if(length(in_set) > 0){
    out_set <- setdiff(cor_cols,c(in_set,gene_target))
    
    dset_results <- list()
    for (dset in names(cor_data)){
      
      sim_net <- cor_data[[dset]]
      
      #KS
      x <- sim_net[in_set,gene_target]
      y <- sim_net[out_set,gene_target]
      ks_pval <- ks.test(x,y,alternative="less")$p.value
      ks_ABSpval <- ks.test(abs(x),abs(y),alternative="less")$p.value
      
      #Min rank of related gene
      codeps <- sort(abs(sim_net[,gene_target]),decreasing=T)
      codeps <- codeps[codeps < 1]
      codeps <- data.frame(codep=names(codeps),rank=1:length(codeps))
      codeps %<>% subset(.,codep %in% in_set)
      codep_minrank_abs <- min(codeps$rank)
      
      codeps <- sort(sim_net[,gene_target],decreasing=T)
      codeps <- codeps[codeps < 1]
      codeps <- data.frame(codep=names(codeps),rank=1:length(codeps))
      codeps %<>% subset(.,codep %in% in_set)
      codep_minrank_pos <- min(codeps$rank)
      
      dset_results[[dset]] <- data.frame(target=gene_target,
                                         tech=dset,
                                         num_partners=length(in_set),
                                         ks_pvalue=ks_pval,
                                         ks_ABS_pvalue=ks_ABSpval,
                                         related_min_rank_abs=codep_minrank_abs,
                                         related_min_rank_pos=codep_minrank_pos,
                                         stringsAsFactors = F)
    }
    gene_res <- bind_rows(dset_results)
    pvals[[gene_target]] <- gene_res
  }
}


melted_res <- bind_rows(pvals)

write_csv(melted_res,file.path(data_processed,"codependency_network_enrichment.csv"))





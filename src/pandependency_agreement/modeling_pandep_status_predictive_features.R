# require(matrixStats)

source("src/packages_paths.R")

#Get expression, CRISPR, and RNAi data
dataset <- list("Expression"=fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-18Q4.csv")),
                "Protein"=fread(file.path(data_raw,"depmap-omics-proteomics-normalized-abundance.csv")),
                "CRISPR"= fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")),
                "RNAi" = fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")))
dataset <- lapply(dataset,function(x){x %<>% column_to_rownames(.,var="V1") %>% as.matrix(.); return(x)})

cls <- lapply(dataset,rownames)
cls <- Reduce(intersect,cls)
dataset <- lapply(dataset,function(d){d[cls,]})

basic_stats <- function(d_name,d_list){
  d <- d_list[[d_name]]
  u <- colMeans(d,na.rm=T)
  q1 <- matrixStats::colQuantiles(d,prob = c(0.1),na.rm=T)
  q2 <- matrixStats::colQuantiles(d,prob = c(0.9),na.rm=T)
  sd <- matrixStats::colVars(d,na.rm=T)
  missing <- colSums(is.na(d)) / nrow(d)
  symbols <- gsub(" .*","",colnames(d))
  mult_iso_list <- symbols[duplicated(symbols)]
  mult_iso <- symbols %in% mult_iso_list
  df <- data.frame(gene_id=colnames(d),
                   u=u,
                   q1=q1,
                   q2=q2,
                   sd=sd,
                   missing=missing,
                   mult_iso=mult_iso,
                   stringsAsFactors = F
  )
  colnames(df) <- paste0(d_name,"_",colnames(df))
  df %<>% add_column(.,symbol=symbols,.before=1) %>% remove_rownames(.)
  return(df)
  
}

stats_list <- list()
for (d_name in names(dataset)){
  stats_list[[d_name]] <- basic_stats(d_name,dataset)
}

full_df <- full_join(stats_list[["CRISPR"]],stats_list[["RNAi"]],by="symbol")
full_df %<>% full_join(.,stats_list[["Expression"]],by="symbol")
full_df %<>% full_join(.,stats_list[["Protein"]],by="symbol")

combos <- names(dataset)

while (length(combos) > 1){
  s1 <- combos[1]
  s1_d <- dataset[[s1]]
  combos <- combos[combos != s1]
  for (s2 in combos){
    s2_d <- dataset[[s2]]
    cors_vec <- vector(mode="numeric",length(nrow(full_df)))
    for (i in 1:nrow(full_df)){
      g1 <- full_df[i,paste0(s1,"_gene_id")]
      g2 <- full_df[i,paste0(s2,"_gene_id")]
      if (any(is.na(c(g1,g2)))){
        cors_vec[i] <- NA
      } else {
        cors_vec[i] <- cor(s1_d[,g1],s2_d[,g2],use="pairwise.complete")
      }
    }
    col_name <- paste0(s1,"_",s2,"_cor")
    full_df %<>% add_column(.,!!col_name :=cors_vec)
    print(paste0("Finished ",col_name))
  }
  print(paste0("Finished ",s1))
}

saveRDS(full_df,file.path(data_processed,"modeling_pandep_status_predictive_features.rds"))

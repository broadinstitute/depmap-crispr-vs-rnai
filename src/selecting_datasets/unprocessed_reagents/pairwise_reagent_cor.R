source("src/packages_paths.R")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))

#Load all reagent sets
file_dict <- list("RNAi-Achilles"="lfc-unscaled-rnai-achilles.csv" ,
                  "RNAi-DRIVE"="lfc-unscaled-rnai-drive.csv",
                  "CRISPR-Avana"="lfc-unscaled-crispr-avana.csv",
                  "CRISPR-KY"="lfc-unscaled-crispr-ky.csv")
file_dict <- lapply(file_dict,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="V1")})

#Scale by cell line to reduce correlations driven by batch effects
file_dict <- lapply(file_dict,function(x){scale(x)})

#Intersect genes across all 4 datasets
CRISPR_map <- fread(file.path(data_raw,"reagent-to-gene-map-sgrna.csv"))
RNAi_map <- fread(file.path(data_raw,"reagent-to-gene-map-shrna.csv"))

genes <- intersect(CRISPR_map$entrez_id[CRISPR_map$Avana],CRISPR_map$entrez_id[CRISPR_map$KY])
genes <- intersect(genes,RNAi_map$entrez_id[RNAi_map$Achilles_55k])
genes <- intersect(genes,RNAi_map$entrez_id[RNAi_map$DRIVE])

reagent_dict <- list("CRISPR-Avana"=subset(CRISPR_map,Avana & (entrez_id %in% genes)),
                     "CRISPR-KY"=subset(CRISPR_map,KY & (entrez_id %in% genes)),
                     "RNAi-Achilles"=subset(RNAi_map,(Achilles_55k | Achilles_98k) & (entrez_id %in% genes)),
                     "RNAi-DRIVE"=subset(RNAi_map,DRIVE & (entrez_id %in% genes)))


genes <- entrez_to_cds(genes,hgnc)
gene_set_ids <- extract_entrez(genes)

#Pairwise correlation quantiles
dataset_res <- list()
for (dataset in names(file_dict)){
  reagents <- file_dict[[dataset]] %>% as.data.frame(.)
  reagent_set_CLs <- colnames(reagents)
  
  reagent_map <- reagent_dict[[dataset]]
  reagents <- reagents[reagent_map$reagent,]
  
  reagent_map$entrez_id %<>% as.character()
  
  reagent_map <- reagent_map[order(reagent_map$entrez_id),]
  reagent_map %<>% select(.,reagent,entrez_id)
  reagent_map %<>% subset(.,reagent %in% rownames(reagents))
  reagent_map %<>% unique()
  reagents %<>% rownames_to_column(.,var="reagent")
  reagents <- left_join(reagent_map,reagents,by="reagent")
  
  reagents <- split(reagents,factor(reagents$entrez_id))
  
  gene_res <- list()
  for (g in names(reagents)){
    gene_df <- reagents[[g]]
    if (any(duplicated(gene_df$reagent))){
      stop(paste0("Duplicated reagents for ",g))
    }
    gene_mat <- gene_df %>% remove_rownames() %>% column_to_rownames(.,var="reagent") %>% select(.,-entrez_id) %>% as.matrix()
    reagent_vars <- matrixStats::rowVars(gene_mat,na.rm=T)
    
    stat_median_var <- median(reagent_vars,na.rm=T)
    stat_max_var <- max(reagent_vars,na.rm=T)
    
    gene_cors <- cor(t(gene_mat),use="pairwise.complete")
    gene_cors_vec <- as.numeric(gene_cors[lower.tri(gene_cors)])
    
    stat_median_cor <- median(gene_cors_vec,na.rm=T)
    stat_top_cor <- max(gene_cors_vec,na.rm=T)
    
    max_var_index <- which.max(reagent_vars)
    diag(gene_cors) <- NA
    stat_maxvar_cor <- max(gene_cors[rownames(gene_mat)[max_var_index],],na.rm=T)
    
    res_df <- data.frame(gene=g,
                         median_var=stat_median_var,
                         max_var=stat_max_var,
                         median_cor=stat_median_cor,
                         top_cor=stat_top_cor,
                         maxvar_cor=stat_maxvar_cor,
                         stringsAsFactors = F)
    
    gene_res[[g]] <- res_df
  }
  
  res_table <- bind_rows(gene_res)
  res_table$dataset <- dataset
  
  res_table <- res_table[order(res_table$median_var),]
  res_table$median_var_perc <- 1:nrow(res_table)
  res_table$median_var_perc <- res_table$median_var_perc / nrow(res_table)
  
  res_table <- res_table[order(res_table$max_var),]
  res_table$max_var_perc <- 1:nrow(res_table)
  res_table$max_var_perc <- res_table$max_var_perc / nrow(res_table)
  
  dataset_res[[dataset]] <- res_table
  
}

plot_df <- bind_rows(dataset_res)
plot_df$median_var_bin <- "0-25"
plot_df$median_var_bin[plot_df$median_var_perc > .25] <- "25-50"
plot_df$median_var_bin[plot_df$median_var_perc > .5] <- "50-75"
plot_df$median_var_bin[plot_df$median_var_perc > .75] <- "75-100"
plot_df$median_var_bin <- factor(plot_df$median_var_bin,levels=c("0-25","25-50","50-75","75-100"))

plot_df$max_var_bin <- "0-25"
plot_df$max_var_bin[plot_df$max_var_perc > .25] <- "25-50"
plot_df$max_var_bin[plot_df$max_var_perc > .5] <- "50-75"
plot_df$max_var_bin[plot_df$max_var_perc > .75] <- "75-100"
plot_df$max_var_bin <- factor(plot_df$max_var_bin,levels=c("0-25","25-50","50-75","75-100"))

mypal <- c("CRISPR-KY"="#E64B35FF","RNAi-DRIVE"="#4DBBD5FF","RNAi-Achilles"="#3C5488FF","CRISPR-Avana"="#F39B7FFF")

ggplot(plot_df,aes(x=median_var_bin,y=median_cor,fill=dataset)) +
  geom_boxplot(size=.25,outlier.size = .25) +
  xlab("Median Variance Percentile") +
  ylab("Median Pearson Correlation") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")
ggsave(file.path("figures","selecting_datasets_reagents_pairwise_cors_boxplot.pdf"),height=3,width=3.5)



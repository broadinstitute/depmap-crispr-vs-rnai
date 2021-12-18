
source("src/packages_paths.R")

args <- commandArgs(trailingOnly = TRUE)
datasets <- strsplit(args[1],",")[[1]]
analyses <- strsplit(args[2],",")[[1]]
hgnc_file <- args[3]
out_file <- args[4]

#### Input result files for each dependency metric required for dependency classifications
pipeline_results <- list()

for (d in datasets){
  d_results <- list()
  for (a in analyses){
    d_results[[a]] <- fread(file.path(data_processed,paste0(a,"-",d,".csv")))
  }
  pipeline_results[[d]] <- d_results
}

######################  Filter for relevant columns and join #######################

keep_columns <- list(`dependency-counts`=c("numDeps(0-50)","numDeps(50-100)"),
                     `gene-effect-moments`=c("Skewness"),
                     `dependency-probability-variance`=c("dep_var","high_variance"),
                     `pandependency-score`=c("CE_percentile","Common_Essential"),
                     `gene-effect-LRT`=c("LRT")
)

#start the result table with ensembl id, symbol, name
gene_deps <- as.data.frame(fread(hgnc_file)) %>% 
  dplyr::select(.,entrez_id,symbol,name) %>% subset(.,!is.na(entrez_id))
gene_deps$entrez_id <- as.character(gene_deps$entrez_id)

for (a in names(keep_columns)){
  
  if (a %in% analyses){
    
    for (d in datasets){
      y <- as.data.frame(pipeline_results[[d]][[a]])
      y$Row.name <- extract_entrez(y$Row.name)
      y <- y[,c("Row.name",keep_columns[[a]])]
      colnames(y) <- c("Row.name",paste0(d,"_",keep_columns[[a]]))
      
      gene_deps <- full_join(gene_deps,y,by=c("entrez_id"="Row.name"))
      
    }
    
  }
  
}

gene_deps <- gene_deps[rowSums(!is.na(gene_deps)) > 3,]

#### Label genes in each class T/F
gene_deps$entrez_id %<>% as.character(.)

gene_lists <- list()
#ND
gene_lists[["CRISPR_ND"]] <- subset(gene_deps,`crispr-matched_numDeps(50-100)` == 0)$entrez_id
gene_lists[["RNAi_ND"]] <- subset(gene_deps,`rnai-matched_numDeps(50-100)` == 0)$entrez_id

#PD
gene_lists[["CRISPR_PD"]] <- subset(gene_deps,`crispr-matched_Common_Essential`)$entrez_id
gene_lists[["RNAi_PD"]] <- subset(gene_deps,`rnai-matched_Common_Essential`)$entrez_id

#HVD
gene_lists[["CRISPR_HVD"]] <- subset(gene_deps,`crispr-matched_high_variance`)$entrez_id
gene_lists[["RNAi_HVD"]] <- subset(gene_deps,`rnai-matched_high_variance`)$entrez_id

#SSD
gene_lists[["CRISPR_SSD"]] <- subset(gene_deps,`crispr-matched_LRT` > 100)$entrez_id
gene_lists[["RNAi_SSD"]] <- subset(gene_deps,`rnai-matched_LRT` > 100)$entrez_id

#SD
all_crispr <- unique(c(gene_lists[["CRISPR_ND"]],gene_lists[["CRISPR_PD"]],gene_lists[["CRISPR_SSD"]],gene_lists[["CRISPR_HVD"]]))
all_rnai <- unique(c(gene_lists[["RNAi_ND"]],gene_lists[["RNAi_PD"]],gene_lists[["RNAi_SSD"]],gene_lists[["RNAi_HVD"]]))
gene_lists[["CRISPR_WSD"]] <- subset(gene_deps,!(entrez_id %in% all_crispr))$entrez_id
gene_lists[["RNAi_WSD"]] <- subset(gene_deps,!(entrez_id %in% all_rnai))$entrez_id

class_order <- c("ND","WSD","SSD","HVD","PD")
for (dep_class in rev(c(paste0("RNAi_",class_order),paste0("CRISPR_",class_order)))){
  gene_deps %<>% add_column(.,!!dep_class := gene_deps$entrez_id %in% gene_lists[[dep_class]],.after=3)
}

####### Add the disjoint mappings
gene_deps %<>% add_column(.,CRISPR_class="Weakly Selective Dependency",.after=3)
gene_deps$CRISPR_class[gene_deps$CRISPR_ND] <- "Non-dependency"
gene_deps$CRISPR_class[gene_deps$CRISPR_SSD] <- "Strongly Selective Dependency"
gene_deps$CRISPR_class[gene_deps$CRISPR_PD] <- "Pan-dependency"
gene_deps$CRISPR_class[gene_deps$CRISPR_HVD] <- "High-variance Dependency"
gene_deps$CRISPR_class[gene_deps$CRISPR_SSD & (gene_deps$`crispr-matched_Skewness` < 0)] <- "Strongly Selective Dependency"

gene_deps %<>% add_column(.,RNAi_class="Weakly Selective Dependency",.after=3)
gene_deps$RNAi_class[gene_deps$RNAi_ND] <- "Non-dependency"
gene_deps$RNAi_class[gene_deps$RNAi_SSD] <- "Strongly Selective Dependency"
gene_deps$RNAi_class[gene_deps$RNAi_PD] <- "Pan-dependency"
gene_deps$RNAi_class[gene_deps$RNAi_HVD] <- "High-variance Dependency"
gene_deps$RNAi_class[gene_deps$RNAi_SSD & (gene_deps$`rnai-matched_Skewness` < 0)] <- "Strongly Selective Dependency"

colnames(gene_deps) <- gsub("crispr-matched","CRISPR",colnames(gene_deps))
colnames(gene_deps) <- gsub("rnai-matched","RNAi",colnames(gene_deps))

gene_deps %<>% dplyr::select(.,-CRISPR_high_variance,-RNAi_high_variance)
gene_deps %<>% rename(.,CRISPR_dependency_variance=CRISPR_dep_var,RNAi_dependency_variance=RNAi_dep_var)

gene_deps %<>% dplyr::select(.,-CRISPR_Common_Essential,-RNAi_Common_Essential)
gene_deps %<>% rename(.,CRISPR_pandependency_score=CRISPR_CE_percentile,RNAi_pandependency_score=RNAi_CE_percentile)

write_csv(gene_deps,out_file)

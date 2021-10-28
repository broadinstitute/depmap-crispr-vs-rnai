
# require(pROC)
# require(ggsci)
# require(scales)
# require(grid)

source("src/id_utility.R")

#################################### All datasets ####################################

hgnc <- fread("data/raw/hgnc-complete-set.csv")

#Get reagent mean data for overlapping CLs and genes
file_dict <- readRDS("data/processed/mean_reagent_lfc.rds")
names(file_dict) <- paste0(names(file_dict),"-Mean")
file_dict <- lapply(file_dict,function(x){t(x)})

# proc_dict <- readRDS("/Users/mburger/dynamic-duo-biorxiv/figures/benchmarking/processed/processed_unscaled_gene_effects.rds")
proc_dict <- list("CRISPR-Avana"=fread("data/raw/gene-effect-unscaled-crispr-avana.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                  "CRISPR-KY"=fread("data/raw/gene-effect-unscaled-crispr-ky.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                  "RNAi-DRIVE"=fread("data/raw/gene-effect-unscaled-rnai-drive.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                  "RNAi-Achilles"=fread("data/raw/gene-effect-unscaled-rnai-achilles.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.))

names(proc_dict) <- paste0(names(proc_dict),"-Corrected")
proc_dict <- lapply(proc_dict,function(x){colnames(x) <- extract_entrez(colnames(x)); return(x)})

file_dict <- c(file_dict,proc_dict)

cls <- lapply(file_dict,function(x){rownames(x)})
cls <- Reduce(intersect,cls)
genes <- lapply(file_dict,function(x){colnames(x)})
genes <- Reduce(intersect,genes)

file_dict <- lapply(file_dict,function(x){x[cls,genes]})

#unbiased essential genes
gene_set <- fread("data/raw/control-essential-genes-unbiased.csv")
gene_set %<>% subset(.,unbiased_essential)
gene_set$entrez_id %<>% as.character(.)
gene_set %<>% subset(., entrez_id %in% genes)
pos_set_ids <- gene_set$entrez_id

#Get nonexpressed genes
tpm <- fread("data/raw/depmap-omics-expression-rnaseq-tpm.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
cls <- intersect(cls,rownames(tpm))
genes_exp <- intersect(genes,extract_entrez(colnames(tpm)))
colnames(tpm) <- extract_entrez(colnames(tpm))
tpm <- tpm[cls,genes_exp]
tpm <- tpm < .2

#Calculate AUC
cl_results <- list()
for (cl in cls){
  
  df <- data.frame(cell_line=cl,dataset=names(file_dict),auc=NA,lower=NA,upper=NA)
  rownames(df) <- df$dataset
  
  nonexp <- colnames(tpm)[tpm[cl,]]
  # nonexp <- neg_set_ids
  keep_genes <- union(nonexp,pos_set_ids)
  
  for (d in names(file_dict)){
    tmp_mat <- file_dict[[d]]
    keep_ids <- intersect(keep_genes,colnames(tmp_mat))
    zscore_vec <- tmp_mat[cl,keep_ids]
    class_vec <- keep_ids %in% pos_set_ids
    r <- pROC::roc(class_vec, zscore_vec)
    conf <- pROC::ci(r)
    df[d,"auc"] <- as.numeric(r$auc)
    df[d,"lower"] <- conf[1]
    df[d,"upper"] <- conf[3]
  }
  
  df$mean <- mean(df$auc)
  cl_results[[cl]] <- df
  
}

cl_results <- do.call(rbind,cl_results)
cl_results %<>% remove_rownames(.)

#Plot dot plot ordered by cell line mean
cl_results <- cl_results[order(cl_results$mean),]
cl_results$cell_line %<>% as.character(cl_results$cell_line)
cl_results$cell_line %<>% factor(cl_results$cell_line,levels=cl_results$cell_line)

# cl_results$dataset <- gsub("RNAi-","",cl_results$dataset)
# cl_results$dataset <- gsub("CRISPR-","",cl_results$dataset)
cl_results$dataset <- factor(cl_results$dataset,levels=c("RNAi-Achilles-Mean","RNAi-Achilles-Corrected",
                                                         "RNAi-DRIVE-Mean","RNAi-DRIVE-Corrected",
                                                         "CRISPR-KY-Mean","CRISPR-KY-Corrected",
                                                         "CRISPR-Avana-Mean","CRISPR-Avana-Corrected"))

#################################### Chosen datasets ####################################

cl_test <- cl_results

mypal = ggsci::pal_npg()(9)
#mypal <- mypal[c(1:2,4:5,8,6)]
#names(mypal) <- c("KY","DRIVE","Achilles","Avana","CERES","DEMETER2")
dsets_names <- c("CRISPR-KY-Mean","CRISPR-KY-Corrected",
                 "RNAi-DRIVE-Mean","RNAi-DRIVE-Corrected",
                 "RNAi-Achilles-Mean","RNAi-Achilles-Corrected",
                 "CRISPR-Avana-Mean","CRISPR-Avana-Corrected",
                 "R-Combined-D2")
dsets_cols <- c("gray",mypal[1],
                "gray",mypal[2],
                "gray",mypal[4],
                "gray",mypal[8],
                mypal[6])
names(dsets_cols) <- dsets_names

#### Dataset label for facet
cl_test$experiment <- "Achilles"
cl_test$experiment[grepl("DRIVE",cl_test$dataset)] <- "DRIVE"
cl_test$experiment[grepl("Avana",cl_test$dataset)] <- "Avana"
cl_test$experiment[grepl("KY",cl_test$dataset)] <- "KY"

#draw horizontal mean lines by separate aes data with experiment label indicated
dsets <- unique(as.character(cl_test$dataset))
u <- lapply(dsets,function(x){mean(subset(cl_test,dataset == x)$auc)})
u %<>% unlist(.)
names(u) <- dsets
u_df <- data.frame(means=u,dataset=names(u),stringsAsFactors = F)
u_df$experiment <- "Achilles"
u_df$experiment[grepl("DRIVE",u_df$dataset)] <- "DRIVE"
u_df$experiment[grepl("Avana",u_df$dataset)] <- "Avana"
u_df$experiment[grepl("KY",u_df$dataset)] <- "KY"
u_df$experiment <- factor(u_df$experiment,levels=c("DRIVE","Achilles","KY","Avana"))

cl_test$experiment <- factor(cl_test$experiment,levels=c("DRIVE","Achilles","KY","Avana"))
cl_test$CL_label <- paste0(cl_test$cell_line,"-",cl_test$experiment)

cl_means <- aggregate(dplyr::select(cl_test,auc), list(cl_test$CL_label), mean)
cl_means <- cl_means[order(cl_means$auc),]

cl_test$CL_label <- factor(cl_test$CL_label,levels=cl_means$Group.1)

ggplot(cl_test, aes(x=CL_label, color=dataset,fill=dataset)) +
  # geom_ribbon(data=cl_test,aes(x=CL_label,ymin=lower,ymax=upper),stat="identity",color="black") +
  geom_point(aes(y=auc, color=dataset),size=.5) +
  scale_color_manual(values=dsets_cols) +
  scale_fill_manual(values=dsets_cols) +
  ylab("AUC") +
  xlab("Cell Line") +
  theme_classic(base_size=11) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(1,3,1,1), "lines")) +
  geom_hline(data=u_df,aes(yintercept=means,color=dataset)) +
  facet_grid(.~experiment, scales = "free", space = "free")

ggsave("figures/selecting_datasets_processed_vs_unproc_AUC.pdf",height=3,width=7.6)







#require(reshape)

source("src/packages_paths.R")

#Mean LFC for reagents targeting the same gene -- 7595 genes (entrez_id) by 62 cell lines (DepMap ID)
gene_score_dict <- readRDS(file.path(data_processed,"mean_reagent_lfc.rds"))
genes <- rownames(gene_score_dict[[1]])
cls <- colnames(gene_score_dict[[1]])

# Core essentials
ceg <- fread(file.path(data_raw,"control-essential-genes-core.csv"),sep=",")
ceg <- ceg$gene
gene_set_ids <- extract_entrez(ceg)
ceg_ids <- intersect(gene_set_ids,genes)

#Unbiased essentials
gene_set <- fread(file.path(data_raw,"control-essential-genes-unbiased.csv"))
gene_set %<>% subset(.,unbiased_essential)
gene_set$entrez_id %<>% as.character(.)
gene_set %<>% subset(., entrez_id %in% genes)
ueg_ids <- gene_set$entrez_id

#Get non-essential genes
nonessential.genes <- fread(file.path(data_raw,"control-nonessential-genes.csv"),sep=",")
neg_ids <- intersect(extract_entrez(nonessential.genes$gene),genes) 

#Get non-expressed genes
tpm <- fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-19Q1.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
cls_exp <- intersect(cls,rownames(tpm))
genes_exp <- intersect(genes,extract_entrez(colnames(tpm)))
colnames(tpm) <- extract_entrez(colnames(tpm))
tpm <- tpm[cls,genes_exp]
tpm <- t(tpm)
tpm <- tpm < .2
#Randomly select 50 non-expressed genes per cell line
for (i in 1:ncol(tpm)){
  sampled_entrez <- sample(rownames(tpm)[tpm[,i]],10)
  tpm[,i] <- rownames(tpm) %in% sampled_entrez
}


ceg_file_dict <- lapply(gene_score_dict,function(x){x[ceg_ids,]})
ueg_file_dict <- lapply(gene_score_dict,function(x){x[ueg_ids,]})
neg_file_dict <- lapply(gene_score_dict,function(x){x[neg_ids,]})

plot_vals <- list()
for (dataset in names(gene_score_dict)){
  
  #Gold-standards
  ceg_mat <- ceg_file_dict[[dataset]] 
  ceg_mat %<>% reshape::melt(.)
  ceg_mat$dataset <- dataset
  ceg_mat %<>% dplyr::select(.,value,dataset)
  ceg_mat$group <- "CEG"
  
  neg_mat <- neg_file_dict[[dataset]] 
  neg_mat %<>% reshape::melt(.)
  neg_mat$dataset <- dataset
  neg_mat %<>% dplyr::select(.,value,dataset)
  neg_mat$group <- "NEG"
  
  goldstandard_mat <- rbind(ceg_mat,neg_mat)
  
  #Larger sets
  ueg_mat <- ueg_file_dict[[dataset]] 
  ueg_mat %<>% reshape::melt(.)
  ueg_mat$dataset <- dataset
  ueg_mat %<>% dplyr::select(.,value,dataset)
  ueg_mat$group <- "UEG"
  
  exp_mat <- gene_score_dict[[dataset]]
  exp_mat <- exp_mat[rownames(tpm),colnames(tpm)]
  exp_mat <- data.frame(value=exp_mat[tpm],dataset=dataset,group="Non-expressed",stringsAsFactors = F)
  
  newstandard_mat <- rbind(ueg_mat,exp_mat)
  
  plot_vals[[dataset]] <- rbind(goldstandard_mat,newstandard_mat)
  
}

plot_vals <- bind_rows(plot_vals)
plot_vals %<>% subset(.,!is.na(value))
plot_vals$tech <- "CRISPR"
plot_vals$tech[grepl("RNAi",plot_vals$dataset)] <- "RNAi"
plot_vals$tech %<>% factor(.,levels=c("RNAi","CRISPR"))
plot_vals$dataset <- gsub("CRISPR-","",plot_vals$dataset)
plot_vals$dataset <- gsub("RNAi-","",plot_vals$dataset)
plot_vals$dataset <- factor(plot_vals$dataset,levels = c("Achilles","DRIVE","Avana","KY"))
plot_vals$group <- factor(plot_vals$group,levels = c("CEG","UEG","Non-expressed","NEG"))
plot_vals$box <- paste0(plot_vals$dataset,"-",plot_vals$group)

mypal <- c("UEG"="#DBBD68","NEG"="#766751","Non-expressed"="#D1C5B4","CEG"="#CC6633")

g <- ggplot(plot_vals, aes(x=dataset, y=value, fill=group)) +
  geom_boxplot(size=.25,outlier.size = .25,coef=4) + 
  xlab("") +
  ylab("Mean LFC") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank())
g <- g + facet_grid(. ~ tech, scales = "free", space = "free") +
  theme(legend.position = "none",legend.margin=margin(t=-.9,unit="cm"))

ggsave(plot=g, file.path("figures","efficacy_specificity_mean_reagent_LFC.pdf"),height=3,width=3.5)



#require(matrixStats)
#require(reshape)

source("src/packages_paths.R")

hgnc <- fread("data/raw/hgnc-complete-set.csv")
hgnc$entrez_id %<>% as.character(.)

file_dict <- readRDS(file.path(data_processed,"processed_unscaled_gene_effects.rds"))
cls <- lapply(file_dict,function(x){rownames(x)})
cls <- Reduce(intersect,cls)
gene_ids <- lapply(file_dict,function(x){extract_entrez(colnames(x))})
gene_ids <- Reduce(intersect,gene_ids)

control_scaling <- function(gene_effect_mat,pos_set_ids,neg_set_ids){
  
  colnames(gene_effect_mat) <- extract_entrez(colnames(gene_effect_mat))
  
  neg_mat <- gene_effect_mat[,colnames(gene_effect_mat) %in% neg_set_ids]
  center_vec <- matrixStats::rowMedians(neg_mat,na.rm=T)
  gene_effect_mat <- gene_effect_mat - center_vec
  
  pos_mat <- gene_effect_mat[,colnames(gene_effect_mat) %in% pos_set_ids]
  scaling_vec <-  matrixStats::rowMedians(pos_mat,na.rm=T)
  gene_effect_mat <- (gene_effect_mat / scaling_vec) * -1
  
  return(gene_effect_mat)
}

#Get non-essential genes
nonessential.genes <- fread(file.path(data_raw,"control-nonessential-genes.csv"),sep=",")
neg_set_ids <- extract_entrez(nonessential.genes$gene) 

#Get core essential genes
ceg <- fread(file.path(data_raw,"control-essential-genes-core.csv"),sep=",")
ceg <- ceg$gene
pos_set_ids <- extract_entrez(ceg)

file_dict <- lapply(file_dict,function(x){control_scaling(x,pos_set_ids,neg_set_ids)})

#unbiased essential genes
gene_set <- fread(file.path(data_raw,"control-essential-genes-unbiased.csv"))
gene_set %<>% subset(.,unbiased_essential)
gene_set$entrez_id %<>% as.character(.)
gene_set %<>% subset(., entrez_id %in% gene_ids)
gene_set_ids <- gene_set$entrez_id

ueg_file_dict <- lapply(file_dict,function(x){x[,colnames(x) %in% gene_set_ids]})

#Get non-expressed genes
tpm <- fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-19Q1.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
cls_exp <- intersect(cls,rownames(tpm))
genes_exp <- intersect(gene_ids,extract_entrez(colnames(tpm)))
colnames(tpm) <- extract_entrez(colnames(tpm))
tpm <- tpm[cls_exp,genes_exp]
tpm <- tpm < .2

plot_vals <- list()
for (dataset in names(file_dict)){
  
  #Larger sets
  ueg_mat <- ueg_file_dict[[dataset]] 
  ueg_mat %<>% reshape::melt(.)
  ueg_mat$dataset <- dataset
  ueg_mat %<>% select(.,value,dataset)
  ueg_mat$group <- "UEG"
  
  exp_mat <- file_dict[[dataset]]
  colnames(exp_mat) <- colnames(exp_mat)
  exp_mat <- exp_mat[rownames(tpm),colnames(tpm)]
  exp_mat <- data.frame(value=exp_mat[tpm],dataset=dataset,group="Non-expressed",stringsAsFactors = F)
  
  newstandard_mat <- rbind(ueg_mat,exp_mat)
  
  plot_vals[[dataset]] <- newstandard_mat
  
}
plot_vals <- bind_rows(plot_vals)
plot_vals %<>% subset(.,!is.na(value))

plot_vals$tech <- "CRISPR"
plot_vals$tech[grepl("RNAi",plot_vals$dataset)] <- "RNAi"
plot_vals$tech %<>% factor(.,levels=c("RNAi","CRISPR"))
plot_vals$dataset <- gsub("CRISPR-","",plot_vals$dataset)
plot_vals$dataset <- gsub("RNAi-","",plot_vals$dataset)
plot_vals$dataset <- factor(plot_vals$dataset,levels = c("Achilles","DRIVE","Avana","KY"))
plot_vals$group <- factor(plot_vals$group,levels = c("UEG","Non-expressed"))
plot_vals$box <- paste0(plot_vals$dataset,"-",plot_vals$group)


mypal <- c("UEG"="#DBBD68","NEG"="#766751","Non-expressed"="#D1C5B4","CEG"="#CC6633")


g <- ggplot(plot_vals, aes(x=dataset, y=value, fill=group)) +
  geom_hline(yintercept = 0,linetype="dashed",color="#766751") +
  geom_hline(yintercept = -1,linetype="dashed",color="#CC6633") +
  geom_violin(width=.9,lwd=.5,draw_quantiles=c(.5)) + 
  xlab("") +
  ylab("Mean Gene Effect (Scaled)") +
  theme_bw(base_size=11) +
  scale_color_manual(values=mypal) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank())
g <- g + facet_grid(. ~ tech, scales = "free", space = "free") +
  theme(legend.position = "none",legend.margin=margin(t=-.9,unit="cm"))

ggsave(plot=g,file.path("figures","efficacy_specificity_processed_scaled_benchmark_sets.pdf"),height=2.7,width=3.4)


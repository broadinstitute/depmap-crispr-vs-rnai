
#################################### All datasets ####################################

source("src/packages_paths.R")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))

#Get reagent mean unprocessed and CERES or D2 processed data for overlapping CLs and genes
file_dict <- readRDS(file.path(data_processed,"mean_reagent_lfc.rds"))
names(file_dict) <- paste0(names(file_dict),"-Mean")
file_dict <- lapply(file_dict,function(x){t(x)})

proc_dict <- list("CRISPR-Avana"="gene-effect-unscaled-crispr-avana.csv",
                  "CRISPR-KY"="gene-effect-unscaled-crispr-ky.csv",
                  "RNAi-DRIVE"="gene-effect-unscaled-rnai-drive.csv",
                  "RNAi-Achilles"="gene-effect-unscaled-rnai-achilles.csv")
proc_dict <- lapply(proc_dict,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)})

names(proc_dict) <- paste0(names(proc_dict),"-Corrected")
proc_dict <- lapply(proc_dict,function(x){colnames(x) <- extract_entrez(colnames(x)); return(x)})

file_dict <- c(file_dict,proc_dict)

cls <- lapply(file_dict,function(x){rownames(x)})
cls <- Reduce(intersect,cls)
genes <- lapply(file_dict,function(x){colnames(x)})
genes <- Reduce(intersect,genes)

file_dict <- lapply(file_dict,function(x){x[cls,genes]})

control_scaling <- function(gene_effect_mat,pos_set_ids,neg_set_ids){
  
  # colnames(gene_effect_mat) <- extract_entrez(colnames(gene_effect_mat))
  
  neg_mat <- gene_effect_mat[,colnames(gene_effect_mat) %in% neg_set_ids]
  center_vec <- matrixStats::rowMedians(neg_mat,na.rm=T)
  gene_effect_mat <- gene_effect_mat - center_vec
  
  pos_mat <- gene_effect_mat[,colnames(gene_effect_mat) %in% pos_set_ids]
  scaling_vec <- matrixStats::rowMedians(pos_mat,na.rm=T)
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
saveRDS(file_dict,file=file.path(data_processed,"processed_unprocessed_scaled_gene_effects.rds"))


#Pearson correlation
pairs <- expand.grid(x = names(file_dict), y = names(file_dict))
pairs %<>% subset(., x != y) %>% remove_rownames(.)
pairs <- t(pairs)

for (i in 1:ncol(pairs)){
  pairs[,i] <- sort(pairs[,i])
}
pairs <- t(pairs)
pairs <- unique(pairs)

res_df <- list()
for (p in 1:nrow(pairs)){
  
  d1_name = pairs[p,1]
  d2_name = pairs[p,2]
  p_name <- paste0(d1_name,"_",d2_name)
  
  
  d1_gs <- file_dict[[d1_name]]
  d2_gs <- file_dict[[d2_name]]
  
  r <- vector(mode="numeric",length(genes))
  for (i in 1:length(genes)){
    r[i] <- cor(d1_gs[,i],d2_gs[,i],use="pairwise.complete.obs")
  }
  
  df <- data.frame(entrez_id=genes,D1=d1_name,D2=d2_name,pearson_r=r,D1_var=matrixStats::colVars(d1_gs,na.rm=T),D2_var=matrixStats::colVars(d2_gs,na.rm=T),stringsAsFactors = F)
  df$D1_perc <- (rank(df$D1_var) / nrow(df)) * 100
  df$D2_perc <- (rank(df$D2_var) / nrow(df)) * 100
  df$D_mean <- rowMeans(select(df,D1_perc,D2_perc))
  df$mean_perc <- (rank(df$D_mean) / nrow(df)) * 100
  df$pair <- p_name
  
  res_df[[p_name]] <- df
}

res_df <- do.call(rbind,res_df)
res_df %<>% remove_rownames(.)

plot_vals <- res_df

plot_vals$perc_bin <- "75-100"
plot_vals$perc_bin[plot_vals$mean_perc <= 75] <- "50-75"
plot_vals$perc_bin[plot_vals$mean_perc <= 50] <- "0-50"
plot_vals$perc_bin <- factor(plot_vals$perc_bin,levels=c("0-50","50-75","75-100"))

plot_vals$group <- "RNAi-CRISPR"
plot_vals$group[grepl("RNAi",plot_vals$D1) & grepl("RNAi",plot_vals$D2)] <- "RNAi-RNAi"
plot_vals$group[grepl("CRISPR",plot_vals$D1) & grepl("CRISPR",plot_vals$D2)] <- "CRISPR-CRISPR"

plot_vals$type <- "Mean-Processed"
plot_vals$type[grepl("Mean",plot_vals$D1) & grepl("Mean",plot_vals$D2)] <- "Mean-Mean"
plot_vals$type[grepl("Corrected",plot_vals$D1) & grepl("Corrected",plot_vals$D2)] <- "Processed-Processed"

plot_vals %<>% subset(.,type != "Mean-Processed")

plot_vals$experiment <- plot_vals$pair
plot_vals$experiment <- gsub("RNAi-","",plot_vals$experiment)
plot_vals$experiment <- gsub("CRISPR-","",plot_vals$experiment)
plot_vals$experiment <- gsub("-Mean","",plot_vals$experiment)
plot_vals$experiment <- gsub("-Corrected","",plot_vals$experiment)

plot_vals$status <- "Unprocessed"
plot_vals$status[plot_vals$type == "Processed-Processed"] <- "Processed"

plot_vals$var_bin <- "0-75"
plot_vals$var_bin[plot_vals$perc_bin == "75-100"] <- "75-100"

mypal <- c("Processed"="#00A087FF","Unprocessed"="#7E6148FF" )

plot_vals %<>% subset(., group == "RNAi-CRISPR")
plot_vals$experiment <- factor(plot_vals$experiment,levels=c("KY_Achilles","Avana_Achilles","Avana_DRIVE","KY_DRIVE"))


#Plot figure
ggplot(plot_vals, aes(x=experiment, y=pearson_r, fill=status)) +
  geom_boxplot() + 
  xlab("") +
  ylab("Pearson Correlation") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  # theme(legend.position = "bottom",legend.margin=margin(t=-.3,unit="cm")) +
  # theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none") +
  facet_grid(.~var_bin, scales = "free", space = "free") 
ggsave(file.path("figures","selecting_datasets_proc_vs_unproc_dataset_cor_boxplot.pdf"),height=2.5,width=4)





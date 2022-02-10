
source("src/packages_paths.R")

#################################### All datasets ####################################

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))

file_dict <- readRDS(file.path(data_processed,"processed_unprocessed_scaled_gene_effects.rds"))

cls <- lapply(file_dict,function(x){rownames(x)})
cls <- Reduce(intersect,cls)
genes <- lapply(file_dict,function(x){colnames(x)})
genes <- Reduce(intersect,genes)
file_dict <- lapply(file_dict,function(x){x[cls,genes]})

#Generating pairs
pairs <- expand.grid(x = names(file_dict), y = names(file_dict))
pairs %<>% subset(., x != y) %>% remove_rownames(.)
pairs$group <- "RNAi-CRISPR"
pairs$group[grepl("RNAi",pairs$x) & grepl("RNAi",pairs$y)] <- "RNAi-RNAi"
pairs$group[grepl("CRISPR",pairs$x) & grepl("CRISPR",pairs$y)] <- "CRISPR-CRISPR"
pairs %<>% subset(.,group == "RNAi-CRISPR")

pairs$type <- "Mean-Processed"
pairs$type[grepl("Mean",pairs$x) & grepl("Mean",pairs$y)] <- "Mean-Mean"
pairs$type[grepl("Corrected",pairs$x) & grepl("Corrected",pairs$y)] <- "Processed-Processed"
pairs %<>% subset(.,type != "Mean-Processed")
pairs <- t(pairs)

for (i in 1:ncol(pairs)){
  pairs[1:2,i] <- sort(pairs[1:2,i])
}
pairs <- t(pairs)
pairs <- unique(pairs)


#Pearson correlation
res_df <- list()
for (p in 1:nrow(pairs)){
  
  d1_name = pairs[p,1]
  d2_name = pairs[p,2]
  p_name <- paste0(d1_name,"_",d2_name)
  
  
  d1_gs <- file_dict[[d1_name]]
  d2_gs <- file_dict[[d2_name]]
  
  r <- cor(d1_gs,d2_gs,use="pairwise.complete.obs")
  
  #diagonal is gene-gene
  result_table <- data.frame(entrez_id=genes,D1=d1_name,D2=d2_name,r=diag(r),stringsAsFactors = F)
  
  #rank cors and find match
  #For each gene in RNAi correlate to all genes in CRISPR and determine how high the gene itself ranks
  rank_d1 <- r
  for (i in 1:ncol(rank_d1)){
    rank_d1[,i] <- frankv(rank_d1[,i],order=-1)
  }
  rank_d1 <- diag(rank_d1)
  result_table$D1_rank <- rank_d1
  
  #For each gene in CRISPR correlate to all genes in RNAi and determine how high the gene itself ranks
  rank_d2 <- r
  for (i in 1:ncol(rank_d2)){
    rank_d2[i,] <- frankv(rank_d2[i,],order=-1)
  }
  rank_d2 <- diag(rank_d2)
  result_table$D2_rank <- rank_d2
  
  result_table$min_rank <- matrixStats::rowMins(as.matrix(select(result_table,D1_rank,D2_rank)))
  
  
  tmp_df <- data.frame(pair=p_name,top_cors=sum(result_table$min_rank == 1,na.rm=T),stringsAsFactors = F)
  
  res_df[[p_name]] <- tmp_df
}

res_df <- do.call(rbind,res_df)
res_df %<>% remove_rownames(.)

plot_vals <- res_df
plot_vals$d_label <- ""
plot_vals$group <- ""
for (i in 1:nrow(plot_vals)){
  tmp_pair <- strsplit(plot_vals$pair[i],"_")[[1]]
  d1 <- tmp_pair[1]
  d2 <- tmp_pair[2]
  
  type_d1 <- strsplit(d1,"-")[[1]]
  type_d2 <- strsplit(d2,"-")[[1]]
  
  plot_vals$d_label[i] <- paste0(type_d1[2],"-",type_d2[2])
  plot_vals$group[i] <- paste0(type_d1[3],"-",type_d2[3])
  
}

plot_vals$status <- "Unprocessed"
plot_vals$status[plot_vals$group == "Corrected-Corrected"] <- "Processed"

mypal <- c("Processed"="#00A087FF","Unprocessed"="#7E6148FF" )

plot_vals$d_label <- factor(plot_vals$d_label,levels=c("KY-Achilles","Avana-Achilles","Avana-DRIVE","KY-DRIVE"))

#Plot figure
ggplot(plot_vals, aes(x=d_label, y=top_cors, fill=status)) +
  geom_bar(stat="identity",position="dodge",col="black") + 
  xlab("") +
  ylab("Number of Genes\nTop Correlation Rank") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") 
ggsave(file.path("figures","selecting_datasets_proc_vs_unprc_topcor_recall.pdf"),height=2.5,width=2.7)





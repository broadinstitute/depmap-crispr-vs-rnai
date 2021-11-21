source("src/packages_paths.R")

hgnc <- fread("data/raw/hgnc-complete-set.csv")

plot_params <- readRDS(file.path(data_processed,"mean_reagent_combRNAi_boxplot_params.rds"))

file_dict <- list("CRISPR"="gene-effect-unscaled-crispr-matched.csv",
                  "RNAi"="gene-effect-unscaled-rnai-matched.csv")
file_dict <- lapply(file_dict,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)})


file_dict <- lapply(file_dict,function(x){x[rownames(x) %in% plot_params$CLs,]})
file_dict <- lapply(file_dict,function(x){colnames(x) <- extract_entrez(colnames(x));return(x)})

file_dict <- lapply(file_dict,function(x){t(scale(t(x)))})

ceg_file_dict <- lapply(file_dict,function(x){x[,colnames(x) %in% plot_params$CEG]})
ueg_file_dict <- lapply(file_dict,function(x){x[,colnames(x) %in% plot_params$UEG]})
neg_file_dict <- lapply(file_dict,function(x){x[,colnames(x) %in% plot_params$NEG]})

#Get non-expressed genes
tpm <- t(plot_params$`Non-expressed`)

#plot all remaining reagents by 4 groups in boxplot
plot_vals <- list()
for (dataset in names(file_dict)){
  
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
  
  exp_mat <- file_dict[[dataset]]
  
  exp_cls <- intersect(rownames(tpm),rownames(exp_mat))
  exp_genes <- intersect(colnames(tpm),colnames(exp_mat))
  exp_mat <- exp_mat[exp_cls,exp_genes]
  tpm <- tpm[exp_cls,exp_genes]
  exp_mat <- data.frame(value=exp_mat[tpm],dataset=dataset,group="Non-expressed",stringsAsFactors = F)
  
  newstandard_mat <- rbind(ueg_mat,exp_mat)
  
  
  
  plot_vals[[dataset]] <- rbind(goldstandard_mat,newstandard_mat)
  
}
plot_vals <- bind_rows(plot_vals)
plot_vals %<>% subset(.,!is.na(value))
plot_vals$tech <- "CRISPR"
plot_vals$tech[grepl("RNAi",plot_vals$dataset)] <- "RNAi"
plot_vals$tech %<>% factor(.,levels=c("RNAi","CRISPR"))
plot_vals$dataset <- factor(plot_vals$dataset,levels = c("RNAi","CRISPR"))
plot_vals$group <- factor(plot_vals$group,levels = c("CEG","UEG","Non-expressed","NEG"))
plot_vals$box <- paste0(plot_vals$dataset,"-",plot_vals$group)


mypal <- c("UEG"="#DBBD68","NEG"="#766751","Non-expressed"="#D1C5B4","CEG"="#CC6633")

g <- ggplot(plot_vals, aes(x=dataset, y=value, fill=group)) +
  geom_boxplot(size=.25,outlier.size = .25,coef=4) +
  xlab("") +
  ylab("Gene Effect Estimates Z-score") +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank())
g<- g + theme(legend.position = "none",legend.margin=margin(t=-.9,unit="cm"))

ggsave(plot=g,file.path("figures","efficacy_specificity_processed_zscored_benchmark_sets_combRNAi.pdf"),height=2.75,width=2.125)



source("src/packages_paths.R")

library(ggridges)
library(ggbeeswarm)

#### Reagent level data
file_dict <- list("RNAi-Achilles"="lfc-unscaled-rnai-achilles.csv",
                  "RNAi-DRIVE"="lfc-unscaled-rnai-drive.csv",
                  "CRISPR-Avana"="lfc-unscaled-crispr-avana.csv",
                  "CRISPR-KY"="lfc-unscaled-crispr-ky.csv")
file_dict <- lapply(file_dict,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="V1")})

CRISPR_map <- fread(file.path(data_raw,"reagent-to-gene-map-sgrna.csv"))
CRISPR_map$entrez_id %<>% as.character(.)
RNAi_map <- fread(file.path(data_raw,"reagent-to-gene-map-shrna.csv"))
RNAi_map$entrez_id %<>% as.character(.)

genes <- intersect(CRISPR_map$entrez_id[CRISPR_map$Avana],CRISPR_map$entrez_id[CRISPR_map$KY])
genes <- intersect(genes,RNAi_map$entrez_id[RNAi_map$Achilles_55k])
genes <- intersect(genes,RNAi_map$entrez_id[RNAi_map$DRIVE])

reagent_dict <- list("CRISPR-Avana"=subset(CRISPR_map,Avana & (entrez_id %in% genes)),
                     "CRISPR-KY"=subset(CRISPR_map,KY & (entrez_id %in% genes)),
                     "RNAi-Achilles"=subset(RNAi_map,(Achilles_55k | Achilles_98k) & (entrez_id %in% genes)),
                     "RNAi-DRIVE"=subset(RNAi_map,DRIVE & (entrez_id %in% genes)))

pos_cntrl <- fread(file.path(data_raw,"control-essential-genes-core.csv"),sep=",")
pos_cntrl <- pos_cntrl$gene
pos_cntrl <- extract_entrez(pos_cntrl)
pos_cntrl <- intersect(pos_cntrl,genes)

neg_cntrl <- fread(file.path(data_raw,"control-nonessential-genes.csv"),sep=",")
neg_cntrl <- neg_cntrl$gene
neg_cntrl <- extract_entrez(neg_cntrl)
neg_cntrl <- intersect(neg_cntrl,genes)

ssmd <- function(p,n){
  result <- (mean(p,na.rm=T) - mean(n,na.rm=T)) / (sqrt(sd(p,na.rm=T)^2 + sd(n,na.rm=T)^2))
  return(result)
}

SSMD_list <- list()
for (dset in names(file_dict)){
  gs_dataset = file_dict[[dset]]
  reagent_list = reagent_dict[[dset]]

  pos_cntrl_reagents <- subset(reagent_list,entrez_id %in% pos_cntrl)$reagent
  neg_cntrl_reagents <- subset(reagent_list,entrez_id %in% neg_cntrl)$reagent

  ssmd_list <- list()
  for (cl in colnames(gs_dataset)){

    pos <- as.numeric(gs_dataset[rownames(gs_dataset) %in% pos_cntrl_reagents,cl])
    neg <- as.numeric(gs_dataset[rownames(gs_dataset) %in% neg_cntrl_reagents,cl])
    ssmd_list[[cl]] <- ssmd(pos,neg)

  }

  SSMD_list[[dset]] <- data.frame(`DepMap_ID`=names(ssmd_list),cell_line_SSMD=unlist(ssmd_list),dataset=dset,stringsAsFactors = F)

}

SSMD_df <- bind_rows(SSMD_list)

cl_mean_SSMD <- aggregate(SSMD_df, list(SSMD_df$DepMap_ID), FUN=mean, na.rm=TRUE)
cl_mean_SSMD <- cl_mean_SSMD[order(cl_mean_SSMD$cell_line_SSMD),]
SSMD_df$DepMap_ID <- factor(SSMD_df$DepMap_ID,levels=cl_mean_SSMD$Group.1)

dataset_pal <- c("CRISPR-Avana"="#F39B7FFF","CRISPR-KY"="#E64B35FF","RNAi-Achilles"="#3C5488FF","RNAi-DRIVE"="#4DBBD5FF")

dsets <- names(file_dict)
u <- lapply(dsets,function(x){mean(subset(SSMD_df,dataset == x)$cell_line_SSMD)})
u %<>% unlist(.)
names(u) <- dsets

#Simplified plot

ggplot(SSMD_df,aes(x=dataset,y=cell_line_SSMD,color=dataset)) +
  geom_quasirandom() +
  # geom_beeswarm() +
  scale_color_manual(values = dataset_pal) +
  theme_bw(base_size=11) +
  ylab("SSMD") +
  xlab("") +
  theme(legend.position = "none")
ggsave(file.path("figures","selecting_datasets_reagents_SSMD_per_dataset.pdf"),width=3.5,height=2.5)


#### Example cell line distributions
dset_cl_values <- list()
for (dset in names(file_dict)){
  gs_dataset = file_dict[[dset]]
  reagent_list = reagent_dict[[dset]]

  pos_cntrl_reagents <- subset(reagent_list,entrez_id %in% pos_cntrl)$reagent
  neg_cntrl_reagents <- subset(reagent_list,entrez_id %in% neg_cntrl)$reagent

  cl <- "ACH-000856"

  pos <- as.numeric(gs_dataset[rownames(gs_dataset) %in% pos_cntrl_reagents,cl])
  neg <- as.numeric(gs_dataset[rownames(gs_dataset) %in% neg_cntrl_reagents,cl])

  pos_df <- data.frame(control_set="CEG",
                          LFC=pos,
                          dataset=dset,
                          stringsAsFactors = F)

  neg_df <- data.frame(control_set="NEG",
                       LFC=neg,
                       dataset=dset,
                       stringsAsFactors = F)

  dset_cl_values[[dset]] <- rbind(pos_df,neg_df)

}

dset_cl_values <- bind_rows(dset_cl_values)

mypal <- c("UEG"="#DBBD68","NEG"="#766751","Non-expressed"="#D1C5B4","CEG"="#CC6633")

ggplot(data=dset_cl_values,aes(x=LFC,y=dataset,fill=control_set,color=control_set))+
  geom_density_ridges(alpha=.7,scale = .95, rel_min_height = .01) +
  scale_shape_identity() +
  # scale_y_discrete(expand = c(.01, 0)) +
  scale_fill_manual(values=mypal) +
  scale_color_manual(values=mypal) +
  theme_ridges() +
  theme(legend.box = 'horizontal') +
  theme(legend.position = "none") +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        plot.title=element_text(size=16,hjust = 0.5)) +
  labs(x="Reagent log2 Fold-change",y="",caption="") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(margin = margin(.25,0,.25,0, "cm"),face = "bold")) +
  theme(plot.margin = margin(l=0,r=0,unit="cm"))
ggsave(file.path("figures","selecting_datasets_reagents_SSMD_CL_example.pdf"),width=4,height=2.5)

#SSMD of the single cell line for each experiment
for (dset in names(reagent_dict)){
  pos_vals <- subset(dset_cl_values,(control_set == "CEG") & (dataset == dset))$LFC
  neg_vals <- subset(dset_cl_values,(control_set == "NEG") & (dataset == dset))$LFC
  print(paste0(dset,": ",ssmd(pos_vals,neg_vals)))
}

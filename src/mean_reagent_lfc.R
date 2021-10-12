# require(ggsci)
# require(scales)

source("src/id_utility.R")

hgnc <- fread("data/raw/hgnc-complete-set.csv")

#Load all unscaled reagent sets
file_dict <- list("RNAi-Achilles"=fread("data/raw/unscaled-rnai-achilles-lfc.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                  "RNAi-DRIVE"=fread("data/raw/unscaled-rnai-drive-lfc.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                  "CRISPR-Avana"=fread("data/raw/unscaled-crispr-broad-lfc.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                  "CRISPR-KY"=fread("data/raw/unscaled-crispr-sanger-lfc.csv") %>% column_to_rownames(.,var="V1") %>% as.matrix(.))

#filter for cell lines that overlap all 4 datasets
cls <- lapply(file_dict,function(x){colnames(x)})
cls <- Reduce(intersect,cls)
file_dict <- lapply(file_dict,function(x){x[,cls]})

#Intersect genes across all 4 datasets
CRISPR_map <- fread("data/raw/reagent-to-gene-map-sgrna.csv")
CRISPR_map$entrez_id %<>% as.character(.)
RNAi_map <- fread("data/raw/reagent-to-gene-map-shrna.csv")
RNAi_map$entrez_id %<>% as.character(.)

genes <- intersect(CRISPR_map$entrez_id[CRISPR_map$Avana],CRISPR_map$entrez_id[CRISPR_map$KY])
genes <- intersect(genes,RNAi_map$entrez_id[RNAi_map$Achilles_55k])
genes <- intersect(genes,RNAi_map$entrez_id[RNAi_map$DRIVE])

#Get reagents that correpsond to overlapping genes
reagent_dict <- list("RNAi-Achilles"=subset(RNAi_map,(Achilles_55k | Achilles_98k) & (entrez_id %in% genes)),
                     "RNAi-DRIVE"=subset(RNAi_map,DRIVE & (entrez_id %in% genes)),
                     "CRISPR-Avana"=subset(CRISPR_map,Avana & (entrez_id %in% genes)),
                     "CRISPR-KY"=subset(CRISPR_map,KY & (entrez_id %in% genes)))

#Filter LFC datasets for reagent list
for (dname in names(file_dict)){
  reagent_map = reagent_dict[[dname]]
  file_dict[[dname]] <- file_dict[[dname]][reagent_map$reagent,]
}

#Aggregate function
mean_narm <- function(x){mean(x,na.rm=T)}

gene_score <- function(tmp_reagents,reagent_map,aggregate_func){
  tmp_reagents <- aggregate(tmp_reagents, list(reagent_map$entrez_id), aggregate_func)
  tmp_reagents %<>% column_to_rownames(.,var="Group.1")
  return(tmp_reagents)
}

gene_score_dict <- lapply(names(file_dict),function(x){gene_score(file_dict[[x]],reagent_dict[[x]],mean_narm)})
names(gene_score_dict) <- names(file_dict)
saveRDS(gene_score_dict,file="data/processed/mean_reagent_lfc.rds")


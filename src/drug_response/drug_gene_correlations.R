
source("src/packages_paths.R")

#Matched genetic screens
genetic_data <- list("crispr"="gene-effect-scaled-crispr-matched.csv",
                     "rnai"="gene-effect-scaled-rnai-matched.csv"
)
genetic_data <- lapply(genetic_data,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

#Dose-level drug data
drug_data <- list("ctd2"="drug-screen-viability-filtered-noPC1-ctd2.csv",
                  "gdsc"="drug-screen-viability-filtered-noPC1-gdsc.csv",
                  "prism"="drug-screen-viability-filtered-noPC1-prism.csv"
)
drug_data <- lapply(drug_data,function(x){fread(file.path(data_processed,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})


#### Gene targets are included in both RNAi datasets and have at least 3 deps in either CRISPR or RNAi
drive <- fread(file.path(data_raw,"gene-effect-unscaled-target-dict.csv")) %>% subset(.,dataset == "RNAi-DRIVE")
drive_genes <- drive$gene

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,(`CRISPR_numDeps(50-100)` >= 3) | (`RNAi_numDeps(50-100)` >= 3))

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv")) 
hgnc$entrez_id %<>% as.character(.)
t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc)
keep_genes <- intersect(t2$CDS_ID,drive_genes)

genetic_data <- lapply(genetic_data,function(x){x[,colnames(x) %in% keep_genes]})

#Long format compound, gene target pairs
target_annotations <- fread(file.path(data_processed,"drug-screen-target-annotations-rephub-shared.csv"))
target_annotations %<>% subset(.,target_id %in% colnames(genetic_data[["crispr"]]))

#### Drug dataset specific
across_drug_datasets <- list()
for (dset in names(drug_data)){
  drug_response <- drug_data[[dset]]
  
  #Filter for shared cell lines
  cls <- intersect(rownames(genetic_data[["crispr"]]),rownames(drug_response))
  crispr <- genetic_data[["crispr"]][cls,]
  rnai <- genetic_data[["rnai"]][cls,]
  drug_response <- drug_response[cls,]
  
  #Filter for perturbations with < 20% missing values
  drug_response <- drug_response[,colSums(is.na(drug_response)) < .2*nrow(drug_response)]
  
  print(paste0(dset,": CL=",length(cls)," genes=",ncol(crispr)," drugs=",length(unique(gsub("::.*","",colnames(drug_response))))))
  
  #Compute all pairwise Pearson correlations between drug and genetic perturbation datasets
  cors_data <- list("CRISPR"=cor(crispr,drug_response,use = 'pairwise.complete'),
                    "RNAi"=cor(rnai,drug_response,use = 'pairwise.complete'))
  
  summary <- list()
  for (tech in names(cors_data)){
    genetic_cors <- cors_data[[tech]]
    
    #Convert every column to a rank -- this gives rank of each gene (rows) per drug perturbation
    genetic_ranks <- apply(as.data.frame(genetic_cors),2,frankv,order=-1,na.last=T)
    rownames(genetic_ranks) <- rownames(genetic_cors)
    
    drug_res <- list()
    #Get the pearson correlation and rank of each gene target per column (drug perturbation)
    for (drug in colnames(genetic_ranks)){
      tmp_targets = subset(target_annotations,column_name == drug)$target_id
      tmp_targets <- tmp_targets[tmp_targets %in% keep_genes]
      
      if (length(tmp_targets >= 1)){
        tmp_cors_mat = genetic_cors[tmp_targets,drug,drop=F]
        tmp_cors_df = reshape2::melt(tmp_cors_mat)
        colnames(tmp_cors_df) <- c("target","RepHub_ID","pearson_r")
        
        tmp_rank_mat = genetic_ranks[tmp_targets,drug,drop=F]
        tmp_rank_df = reshape2::melt(tmp_rank_mat)
        colnames(tmp_rank_df) <- c("target","RepHub_ID","rank")
        
        tmp_cors_df$rank <- tmp_rank_df$rank
        
        drug_res[[drug]] <- tmp_cors_df
      }
      
    }
    drug_res <- bind_rows(drug_res)
    drug_res$broad_id <- gsub("::.*","",drug_res$RepHub_ID)
    drug_res$dose <- gsub(".+::","",drug_res$RepHub_ID)
    
    drug_res$genetic_perturbation <- tech
    drug_res %<>% remove_rownames(.)
    drug_res %<>% rename(.,drug_perturbation=RepHub_ID)
    
    #Label the best correlation for each drug to any of its targets (gives best target+dose)
    drug_res <- drug_res[order(drug_res$pearson_r,decreasing=T),]
    drug_res$best_drug_cor <- !duplicated(drug_res$broad_id)

    #Label the best rank for each drug to any of its targets/doses (gives best target+dose)
    drug_res <- drug_res[order(drug_res$rank,decreasing=F),]
    drug_res$best_drug_rank <- !duplicated(drug_res$broad_id)
    
    drug_res$target_drug_pair <- paste0(drug_res$target,"::",drug_res$broad_id)
    
    #Label the best correlation for each drug to each of its targets (gives best dose per target)
    drug_res <- drug_res[order(drug_res$pearson_r,decreasing=T),]
    drug_res$best_drug_target_cor <- !duplicated(drug_res$target_drug_pair)
    
    #Label the best rank for each drug to each of its targets (gives best dose per target)
    drug_res <- drug_res[order(drug_res$rank,decreasing=F),]
    drug_res$best_drug_target_rank <- !duplicated(drug_res$target_drug_pair)
    
    summary[[tech]] <- drug_res
  }
  
  summary %<>% bind_rows(.)
  summary$drug_dataset <- dset
  
  across_drug_datasets[[dset]] <- summary
}

across_drug_datasets %<>% bind_rows(.)
write_csv(across_drug_datasets,file=file.path(data_processed,"drug-screen-genetic-targets-correlations.csv"))


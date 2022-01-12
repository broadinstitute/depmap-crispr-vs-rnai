
source("src/packages_paths.R")

rephub <- fread(file.path(data_raw,"drug-screen-target-annotations-rephub.csv"))

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

hgnc %<>% subset(.,!is.na(entrez_id)) %>% subset(.,locus_group == "protein-coding gene")
hgnc %<>% dplyr::select(.,entrez_id,symbol)
hgnc$target_id <- paste0(hgnc$symbol," (",hgnc$entrez_id,")")

hgnc %<>% dplyr::select(.,target=symbol,target_id)

rephub %<>% left_join(.,hgnc,by="target")
rephub %<>% subset(.,!is.na(target_id))
rephub %<>% dplyr::select(.,broad_id,target_id)

dose_info <- list(prism="drug-screen-viability-info-prism.csv",
                  ctd2="drug-screen-viability-info-ctd2.csv",
                  gdsc="drug-screen-viability-info-gdsc.csv")
dose_info <- lapply(dose_info,function(x){fread(file.path(data_raw,x))})

matchRelated <- list()
for (dset in names(dose_info)){
  tmp_dose = dose_info[[dset]]
  tmp_dose %<>% dplyr::select(.,broad_id,column_name)
  tmp_dose %<>% subset(.,broad_id %in% rephub$broad_id)
  tmp_dose %<>% left_join(.,rephub,by="broad_id")
  matchRelated[[dset]] <- tmp_dose
}

matchRelated <- bind_rows(matchRelated)
matchRelated %<>% dplyr::select(.,column_name,target_id)
matchRelated$source <- "RepHub"
matchRelated <- unique(matchRelated)

drive <- fread(file.path(data_raw,"gene-effect-unscaled-target-dict.csv")) %>% subset(.,dataset == "RNAi-DRIVE")
drive_genes <- drive$gene

crispr <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.) 
all_targets <- intersect(drive_genes,colnames(crispr))

matchRelated %<>% subset(.,target_id %in% all_targets)

write_csv(matchRelated,file.path(data_processed,"drug-screen-target-annotations-rephub-shared.csv"))

drug_datasets <- list("ctd2"="drug-screen-viability-ctd2.csv",
                      "gdsc"="drug-screen-viability-gdsc.csv",
                      "prism"="drug-screen-viability-prism.csv")
drug_datasets <- lapply(drug_datasets,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})


drug_datasets <- lapply(drug_datasets,function(x){x<-x[,colnames(x) %in% matchRelated$column_name]; return(x)})

for (dset in names(drug_datasets)){
  tmp_data <- drug_datasets[[dset]]
  tmp_data %<>% as.data.frame(.) %>% rownames_to_column(.,var="Row.name")
  write_csv(tmp_data,file.path(data_processed,paste0("drug-screen-viability-filtered-",dset,".csv")))   
}


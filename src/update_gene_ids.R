
data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")
source(file.path("src","id_utility.R"))

# gene ID updates which were performed on RNAi datasets as part of target discovery pipeline
hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
hgnc %<>% subset(.,locus_group == "protein-coding gene")

ach_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-achilles.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
colnames(ach_gs) <- extract_entrez(colnames(ach_gs))
ach_gs <- ach_gs[,colnames(ach_gs) %in% hgnc$entrez_id]
colnames(ach_gs) <- entrez_to_cds(colnames(ach_gs),hgnc)
write.csv(ach_gs,file.path(data_processed,"gene-effect-scaled-rnai-achilles.csv"))

drive_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-drive.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
colnames(drive_gs) <- extract_entrez(colnames(drive_gs))
drive_gs <- drive_gs[,colnames(drive_gs) %in% hgnc$entrez_id]
colnames(drive_gs) <- entrez_to_cds(colnames(drive_gs),hgnc)
write.csv(drive_gs,file.path(data_processed,"gene-effect-scaled-rnai-drive.csv"))


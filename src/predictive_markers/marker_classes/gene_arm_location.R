source("src/packages_paths.R")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

gene_location <- dplyr::select(hgnc,entrez_id,location) %>%
  subset(.,!is.na(entrez_id))
gene_location$arm <- gsub("p.*","",gene_location$location)
gene_location$arm <- gsub("q.*","",gene_location$arm)
gene_location$tmp <- ""
gene_location$tmp[grepl("p",gene_location$location)] <- "p"
gene_location$tmp[grepl("q",gene_location$location)] <- "q"
gene_location$arm_level <- paste0(gene_location$arm,gene_location$tmp)
gene_location %<>% dplyr::select(.,entrez_id,arm_level)
gene_location$gene <- entrez_to_cds(gene_location$entrez_id,hgnc)
gene_location %<>% dplyr::select(.,gene,arm_level)

write_csv(gene_location,file.path("data/processed","hgnc_gene_arm_location.csv"))


setwd("/Users/mburger/dynamic-duo/figures/genedep_profiles/common_essential")
source("/Users/mburger/dynamic-duo/src/R/id_utility.R")

hgnc <- load.from.taiga(data.name='hgnc-6825', data.version=2, data.file='hgnc_complete_set_090318')
hgnc$entrez_id %<>% as.character(.)

## RNAi controls
ceg <- load.from.taiga(data.name='avana-internal-18q3-6b2c', data.version=3, data.file='essential_genes')
ceg <- ceg$gene
neg <- load.from.taiga(data.name='avana-public-tentative-19q1-6956', data.version=2, data.file='nonessential_genes')
neg <- neg$gene

## CRISPR controls
ceg2 <- fread("/Users/mburger/dynamic-duo/data/processed/essential_gene_lists/CEGv2.txt")
ceg2 <- subset(hgnc,hgnc_id %in% ceg2$HGNC_ID)$entrez_id
ceg2 <- entrez_to_cds(ceg2,hgnc)

exp <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=12, data.file='CCLE_depMap_18Q4_TPM_ProteinCoding')
exp <- exp > .2
exp_frac <- colSums(exp,na.rm=T) / colSums(!is.na(exp))
nonexp <- names(exp_frac)[exp_frac < .5]

## D2 scores
achilles <- load.from.taiga(data.name='task-results-ecd9', data.version=1, data.file='Achilles-common-essential')
achilles_pos <- subset(achilles,`Row.name` %in% ceg)
achilles_pos$group <- "CEG"
achilles_pos$dataset <- "Achilles"
achilles_neg <- subset(achilles,`Row.name` %in% neg)
achilles_neg$group <- "NEG"
achilles_neg$dataset <- "Achilles"
achilles_df <- bind_rows(achilles_pos,achilles_neg)

drive <- load.from.taiga(data.name='task-results-ecd9', data.version=1, data.file='DRIVE-common-essential')
drive_pos <- subset(drive,`Row.name` %in% ceg)
drive_pos$group <- "CEG"
drive_pos$dataset <- "DRIVE"
drive_neg <- subset(drive,`Row.name` %in% neg)
drive_neg$group <- "NEG"
drive_neg$dataset <- "DRIVE"
drive_df <- bind_rows(drive_pos,drive_neg)

rnai_df <- bind_rows(achilles_df,drive_df)

## CERES scores
avana <- load.from.taiga(data.name='task-results-59b8', data.version=1, data.file='Avana-common-essential')
avana_pos <- subset(avana,`Row.name` %in% ceg2)
avana_pos$group <- "CEG2"
avana_pos$dataset <- "Avana"
avana_neg <- subset(avana,`Row.name` %in% nonexp)
avana_neg$group <- "Non-expressed"
avana_neg$dataset <- "Avana"
avana_df <- bind_rows(avana_pos,avana_neg)

ky <- load.from.taiga(data.name='task-results-59b8', data.version=1, data.file='KY-common-essential')
ky_pos <- subset(ky,`Row.name` %in% ceg2)
ky_pos$group <- "CEG2"
ky_pos$dataset <- "KY"
ky_neg <- subset(ky,`Row.name` %in% nonexp)
ky_neg$group <- "Non-expressed"
ky_neg$dataset <- "KY"
ky_df <- bind_rows(ky_pos,ky_neg)

crispr_df <- bind_rows(avana_df,ky_df)

mypal <- c("CEG"="#DBBD68","NEG"="#766751","CEG2"="#9986A5","Non-expressed"="#D1C5B4")

ggplot(rnai_df,aes(x=CE_percentile,fill=group)) +
  geom_density() +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  xlab("Ranking Within 90th Percentile Least Dependent Cell Line") +
  theme(legend.position = c(0.85, 0.7)) +
  facet_grid(. ~ dataset, scales = "free", space = "free")
ggsave(paste0("/Users/mburger/dynamic-duo/figures/genedep_profiles/common_essential/benchmark_90th_percentile_ranks_RNAi.pdf"),height=2.5,width=4.5)

ggplot(crispr_df,aes(x=CE_percentile,fill=group)) +
  geom_density() +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  xlab("") +
  theme(legend.position = c(0.85, 0.7)) +
  facet_grid(. ~ dataset, scales = "free", space = "free")
ggsave(paste0("/Users/mburger/dynamic-duo/figures/genedep_profiles/common_essential/benchmark_90th_percentile_ranks_CRISPR.pdf"),height=2.5,width=4.5)


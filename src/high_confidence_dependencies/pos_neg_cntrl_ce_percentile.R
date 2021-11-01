

data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")
source(file.path("src","id_utility.R"))

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

## RNAi controls
ceg <- fread(file.path(data_raw,"control-essential-genes-core.csv"),sep=",")
ceg <- ceg$gene
neg <- fread(file.path(data_raw,"control-nonessential-genes.csv"),sep=",")
neg <- neg$gene

## CRISPR controls
ceg2 <- fread(file.path(data_raw,"control-essential-genes-CEGv2.csv"),sep=",")
ceg2 <- ceg2$gene

exp <- fread(file.path(data_raw,"depmap-omics-expression-rnaseq-tpm-19Q1.csv"))
exp <- exp > .2
exp_frac <- colSums(exp,na.rm=T) / colSums(!is.na(exp))
nonexp <- names(exp_frac)[exp_frac < .5]

## 90th percentile rank results
pd_score <- fread(file.path(data_processed,"multilib_ce_percentile_results.csv"))

## D2 scores
achilles <- subset(pd_score,dataset == "rnai_achilles")
achilles_pos <- subset(achilles,CDS_ID %in% ceg)
achilles_pos$group <- "CEG"
achilles_pos$dataset <- "Achilles"
achilles_neg <- subset(achilles,CDS_ID %in% neg)
achilles_neg$group <- "NEG"
achilles_neg$dataset <- "Achilles"
achilles_df <- bind_rows(achilles_pos,achilles_neg)

drive <- subset(pd_score,dataset == "rnai_drive")
drive_pos <- subset(drive,CDS_ID %in% ceg)
drive_pos$group <- "CEG"
drive_pos$dataset <- "DRIVE"
drive_neg <- subset(drive,CDS_ID %in% neg)
drive_neg$group <- "NEG"
drive_neg$dataset <- "DRIVE"
drive_df <- bind_rows(drive_pos,drive_neg)

rnai_df <- bind_rows(achilles_df,drive_df)

## CERES scores
avana <- subset(pd_score,dataset == "crispr_avana")
avana_pos <-subset(avana,CDS_ID %in% ceg2)
avana_pos$group <- "CEG2"
avana_pos$dataset <- "Avana"
avana_neg <- subset(avana,CDS_ID %in% nonexp)
avana_neg$group <- "Non-expressed"
avana_neg$dataset <- "Avana"
avana_df <- bind_rows(avana_pos,avana_neg)

ky<- subset(pd_score,dataset == "crispr_ky")
ky_pos <-subset(ky,CDS_ID %in% ceg2)
ky_pos$group <- "CEG2"
ky_pos$dataset <- "KY"
ky_neg <- subset(ky,CDS_ID %in% nonexp)
ky_neg$group <- "Non-expressed"
ky_neg$dataset <- "KY"
ky_df <- bind_rows(ky_pos,ky_neg)

crispr_df <- bind_rows(avana_df,ky_df)

mypal <- c("CEG"="#DBBD68","NEG"="#766751","CEG2"="#9986A5","Non-expressed"="#D1C5B4")

ggplot(rnai_df,aes(x=rank,fill=group)) +
  geom_density() +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  xlab("Ranking Within 90th Percentile Least Dependent Cell Line") +
  theme(legend.position = c(0.85, 0.7)) +
  facet_grid(. ~ dataset, scales = "free", space = "free")
ggsave(file.path("figures","pandependency_90th_percentile_ranks_RNAi_benchmark.pdf"),height=2.5,width=4.5)

ggplot(crispr_df,aes(x=rank,fill=group)) +
  geom_density() +
  theme_bw(base_size=11) +
  scale_fill_manual(values=mypal) +
  theme(legend.title=element_blank()) +
  xlab("") +
  theme(legend.position = c(0.85, 0.7)) +
  facet_grid(. ~ dataset, scales = "free", space = "free")
ggsave(file.path("figures","pandependency_90th_percentile_ranks_CRISPR_benchmark.pdf"),height=2.5,width=4.5)


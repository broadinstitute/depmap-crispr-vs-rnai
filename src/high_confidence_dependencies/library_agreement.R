
data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")
source(file.path("src","id_utility.R"))

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
#hgnc %<>% subset(.,locus_group == "protein-coding gene")

#Individual library results
avana_pandep_scores <- fread(file.path(data_processed,"pandependency-score-crispr-avana.csv"),sep=",")
avana_pandep_scores %<>% dplyr::rename(.,Avana_Common_Essential=Common_Essential)
avana_pandep_scores$entrez_id <- extract_entrez(avana_pandep_scores$Row.name)
avana_pandep_scores %<>% dplyr::select(.,entrez_id,Avana_Common_Essential)
avana_target <- avana_pandep_scores$entrez_id

ky_pandep_scores <- fread(file.path(data_processed,"pandependency-score-crispr-ky.csv"),sep=",")
ky_pandep_scores %<>% dplyr::rename(.,KY_Common_Essential=Common_Essential)
ky_pandep_scores$entrez_id <- extract_entrez(ky_pandep_scores$Row.name)
ky_pandep_scores %<>% dplyr::select(.,entrez_id,KY_Common_Essential)
ky_target <- ky_pandep_scores$entrez_id

achilles_pandep_scores <- fread(file.path(data_processed,"pandependency-score-rnai-achilles.csv"),sep=",")
achilles_pandep_scores %<>% dplyr::rename(.,Achilles_Common_Essential=Common_Essential)
achilles_pandep_scores$entrez_id <- extract_entrez(achilles_pandep_scores$Row.name)
achilles_pandep_scores %<>% dplyr::select(.,entrez_id,Achilles_Common_Essential)
achilles_target <- achilles_pandep_scores$entrez_id

drive_pandep_scores <- fread(file.path(data_processed,"pandependency-score-rnai-drive.csv"),sep=",")
drive_pandep_scores %<>% dplyr::rename(.,DRIVE_Common_Essential=Common_Essential)
drive_pandep_scores$entrez_id <- extract_entrez(drive_pandep_scores$Row.name)
drive_pandep_scores %<>% dplyr::select(.,entrez_id,DRIVE_Common_Essential)
drive_target <- drive_pandep_scores$entrez_id

# Columns of T/F for pandep in each library
crispr_pandep_df <- full_join(avana_pandep_scores,ky_pandep_scores,by="entrez_id")
crispr_pandep_df$symbol <- entrez_to_symbol(crispr_pandep_df$entrez_id,hgnc)
crispr_pandep_df$CRISPR_avana_target <- crispr_pandep_df$entrez_id %in% avana_target
crispr_pandep_df$CRISPR_ky_target <- crispr_pandep_df$entrez_id %in% ky_target
crispr_pandep_df %<>% dplyr::select(.,entrez_id,symbol,CRISPR_avana_target,CRISPR_ky_target,CRISPR_avana_pandep=Avana_Common_Essential,CRISPR_ky_pandep=KY_Common_Essential)

rnai_pandep_df <- full_join(achilles_pandep_scores,drive_pandep_scores,by="entrez_id")
rnai_pandep_df$symbol <- entrez_to_symbol(rnai_pandep_df$entrez_id,hgnc)
rnai_pandep_df$RNAi_achilles_target <- rnai_pandep_df$entrez_id %in% achilles_target
rnai_pandep_df$RNAi_drive_target <- rnai_pandep_df$entrez_id %in% drive_target
rnai_pandep_df %<>% dplyr::select(.,entrez_id,symbol,RNAi_achilles_target,RNAi_drive_target,RNAi_achilles_pandep=Achilles_Common_Essential,RNAi_drive_pandep=DRIVE_Common_Essential)

# Column of T/F for whether it is top correlate between 2 libraries
rnai_cor <- fread(file.path(data_processed,"RNAi-Achilles_vs_RNAi-DRIVE_cor_recall_vals.csv"))
rnai_cor$entrez_id %<>% as.character(.)
rnai_cor$RNAi_cor_rank <- matrixStats::rowMins(as.matrix(select(rnai_cor,`rank_in_RNAi-Achilles`,`rank_in_RNAi-DRIVE`)))
rnai_cor %<>% select(.,entrez_id,RNAi_cor_rank)

crispr_cor <- fread(file.path(data_processed,"CRISPR-Avana_vs_CRISPR-KY_cor_recall_vals.csv"))
crispr_cor$entrez_id %<>% as.character(.)
crispr_cor$CRISPR_cor_rank <- matrixStats::rowMins(as.matrix(select(crispr_cor,`rank_in_CRISPR-Avana`,`rank_in_CRISPR-KY`)))
crispr_cor %<>% select(.,entrez_id,CRISPR_cor_rank)

# Columns of T/F for no dep in each library
ach_pr <- fread(file.path(data_processed,"dependency-probability-rnai-achilles.csv"))
ach_pr %<>% column_to_rownames(.,var="Row.name")
colnames(ach_pr) <- extract_entrez(colnames(ach_pr))
ach_pr <- ach_pr[,colnames(ach_pr) %in% hgnc$entrez_id]
ach_dep <- ach_pr > .5
ach_depCount <- colSums(ach_dep,na.rm=T)
ach_depCount_df <- data.frame(entrez_id=names(ach_depCount),RNAi_achilles_nodep=ach_depCount == 0,stringsAsFactors = F)

drive_pr <- fread(file.path(data_processed,"dependency-probability-rnai-drive.csv"))
drive_pr %<>% column_to_rownames(.,var="Row.name")
colnames(drive_pr) <- extract_entrez(colnames(drive_pr))
drive_pr <- drive_pr[,colnames(drive_pr) %in% hgnc$entrez_id]
drive_dep <- drive_pr > .5
drive_depCount <- colSums(drive_dep,na.rm=T)
drive_depCount_df <- data.frame(entrez_id=names(drive_depCount),RNAi_drive_nodep=drive_depCount == 0,stringsAsFactors = F)
rnai_depCount_df <- full_join(ach_depCount_df,drive_depCount_df,by="entrez_id")

avana_pr <- fread(file.path(data_processed,"dependency-probability-crispr-avana.csv"))
avana_pr %<>% column_to_rownames(.,var="Row.name")
colnames(avana_pr) <- extract_entrez(colnames(avana_pr))
avana_pr <- avana_pr[,colnames(avana_pr) %in% hgnc$entrez_id]
avana_dep <- avana_pr > .5
avana_depCount <- colSums(avana_dep,na.rm=T)
avana_depCount_df <- data.frame(entrez_id=names(avana_depCount),CRISPR_avana_nodep=avana_depCount == 0,stringsAsFactors = F)

ky_pr <- fread(file.path(data_processed,"dependency-probability-crispr-ky.csv"))
ky_pr %<>% column_to_rownames(.,var="Row.name")
colnames(ky_pr) <- extract_entrez(colnames(ky_pr))
ky_pr <- ky_pr[,colnames(ky_pr) %in% hgnc$entrez_id]
ky_dep <- ky_pr > .5
ky_depCount <- colSums(ky_dep,na.rm=T)
ky_depCount_df <- data.frame(entrez_id=names(ky_depCount),CRISPR_ky_nodep=ky_depCount == 0,stringsAsFactors = F)
crispr_depCount_df <- full_join(avana_depCount_df,ky_depCount_df,by="entrez_id")

# Column for number of deps in shared cell lines between 2 libraries
rnai_genes <- intersect(colnames(ach_dep),colnames(drive_dep))
rnai_cls <- intersect(rownames(ach_dep),rownames(drive_dep))
ach_shared <- ach_dep[rnai_cls,rnai_genes]
drive_shared <- drive_dep[rnai_cls,rnai_genes]
rnai_dep <- ach_shared & drive_shared
rnai_hits <- colSums(rnai_dep,na.rm=T)
rnai_hits_df <- data.frame(entrez_id=names(rnai_hits),RNAi_overlap_deps=rnai_hits,stringsAsFactors = F)

rnai_overlap <- colSums(!is.na(rnai_dep))
if (all(names(rnai_overlap) == rnai_hits_df$entrez_id)){
  rnai_hits_df$RNAi_shared_CLs <- rnai_overlap
} else {
  stop("Entrez ID don't match")
}


crispr_genes <- intersect(colnames(avana_dep),colnames(ky_dep))
crispr_cls <- intersect(rownames(avana_dep),rownames(ky_dep))
avana_shared <- avana_dep[crispr_cls,crispr_genes]
ky_shared <- ky_dep[crispr_cls,crispr_genes]
crispr_dep <- avana_shared & ky_shared
crispr_hits <- colSums(crispr_dep,na.rm=T)
crispr_hits_df <- data.frame(entrez_id=names(crispr_hits),CRISPR_overlap_deps=crispr_hits,stringsAsFactors = F)

crispr_overlap <- colSums(!is.na(crispr_dep))
if (all(names(crispr_overlap) == crispr_hits_df$entrez_id)){
  crispr_hits_df$CRISPR_shared_CLs <- crispr_overlap
} else {
  stop("Entrez ID don't match")
}

t3 <- full_join(crispr_pandep_df,crispr_cor,by="entrez_id")
t3 %<>% full_join(.,crispr_depCount_df,by="entrez_id")
t3 %<>% full_join(.,crispr_hits_df,by="entrez_id") ## x
t3$CRISPR_shared_CLs[is.na(t3$CRISPR_shared_CLs)] <- 0 ## x

t3$CRISPR_pandep <- rowSums(t3[,c("CRISPR_avana_pandep","CRISPR_ky_pandep")],na.rm=T)
t3$CRISPR_topcor <- t3$CRISPR_cor_rank == 1
t3$CRISPR_nodep <- rowSums(t3[,c("CRISPR_avana_nodep","CRISPR_ky_nodep")],na.rm=T)
t3$CRISPR_agreement <- (t3$CRISPR_pandep == 2) | t3$CRISPR_topcor | (t3$CRISPR_nodep == 2)
t3$CRISPR_no_shared_depCLs <- t3$CRISPR_overlap_deps == 0 ## x
t3$CRISPR_not_shared_target <- !(t3$CRISPR_avana_target & t3$CRISPR_ky_target)
t3$CRISPR_agreement[t3$CRISPR_not_shared_target] <- NA

t3_crispr <- t3
t3_crispr %<>% select(.,-symbol)

t3 <- full_join(rnai_pandep_df,rnai_cor,by="entrez_id")
t3 %<>% full_join(.,rnai_depCount_df,by="entrez_id")
t3 %<>% full_join(.,rnai_hits_df,by="entrez_id")
t3$RNAi_shared_CLs[is.na(t3$RNAi_shared_CLs)] <- 0

t3$RNAi_pandep <- rowSums(t3[,c("RNAi_achilles_pandep","RNAi_drive_pandep")],na.rm=T)
t3$RNAi_topcor <- t3$RNAi_cor_rank == 1
t3$RNAi_nodep <- rowSums(t3[,c("RNAi_achilles_nodep","RNAi_drive_nodep")],na.rm=T)
t3$RNAi_agreement <- (t3$RNAi_pandep == 2) | t3$RNAi_topcor | (t3$RNAi_nodep == 2)
t3$RNAi_no_shared_depCLs <- t3$RNAi_overlap_deps == 0
t3$RNAi_not_shared_target <- !(t3$RNAi_achilles_target & t3$RNAi_drive_target)
t3$RNAi_agreement[t3$RNAi_not_shared_target] <- NA

t3_rnai <- t3
t3_rnai %<>% select(.,-symbol)

t3 <- full_join(t3_crispr,t3_rnai,by="entrez_id")
t3 %<>% add_column(.,symbol=entrez_to_symbol(t3$entrez_id,hgnc),.before=1)

write_csv(t3,file.path(data_processed,"library_agreement.csv"))

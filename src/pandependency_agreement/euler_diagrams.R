# require(eulerr)

source("src/packages_paths.R")

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

eg_summary <- fread(file.path(data_raw,"control-essential-genes-unbiased.csv"))
eg_summary$entrez_id %<>% as.character(.)

#Call top 15% of ExAC metric medians a hit
eg_summary$exac_rank <- frankv(eg_summary$human_median,order=-1)
eg_summary$exac_rank[is.na(eg_summary$human_median)] <- NA
eg_summary$exac_hit <- eg_summary$exac_rank < (.15 * sum(!is.na(eg_summary$exac_rank)))

#Call top 15% of gene trap significance median a hit
eg_summary$trap_rank <- frankv(eg_summary$gene_trap_percentile_mean,order=-1)
eg_summary$trap_rank[is.na(eg_summary$gene_trap_percentile_mean)] <- NA
eg_summary$trap_hit <- eg_summary$trap_rank < (.15 * sum(!is.na(eg_summary$trap_rank)))

eg_summary %<>% subset(.,entrez_id %in% t1$entrez_id)

t2$mouse_ko <- t2$entrez_id %in% subset(eg_summary,mouse_essential)$entrez_id
t2$exac_hit <- t2$entrez_id %in% subset(eg_summary,exac_hit)$entrez_id
t2$trap_hit <- t2$entrez_id %in% subset(eg_summary,trap_hit)$entrez_id

t2$evidence <- rowSums(t2[,c("mouse_ko","exac_hit","trap_hit")],na.rm=T)

t2$support <- t2$evidence >= 2

t2$pandep_group <- "unknown"
t2$pandep_group[(!t2$RNAi_PD) & (!t2$CRISPR_PD)] <- "baseline"
t2$pandep_group[(t2$RNAi_PD) & (t2$CRISPR_PD)] <- "shared"
t2$pandep_group[(!t2$RNAi_PD) & (t2$CRISPR_PD)] <- "CRISPR-specific"
t2 %<>% subset(.,pandep_group %in% c("baseline","shared","CRISPR-specific"))

t2$evidence %<>% as.factor(.)

t2$pandep_group <- gsub("CRISPR-specific","CRISPR-spec.",t2$pandep_group)

#CRISPR support
crispr_mat <- dplyr::select(t2,CRISPR=CRISPR_PD,Mus=mouse_ko,ExAC=exac_hit,GT=trap_hit)
crispr_mat <- as.matrix(crispr_mat)
fit_crispr <- eulerr::euler(crispr_mat)
pdf(file.path("figures","pandependency_agreement_CRISPR_euler_venn.pdf"), height=3.9,width=4)
plot(fit_crispr)
dev.off()
nrow(subset(t2,(CRISPR_PD & support))) / nrow(subset(t2,CRISPR_PD))

#RNAi support
rnai_mat <- dplyr::select(t2,RNAi=RNAi_PD,Mus=mouse_ko,ExAC=exac_hit,GT=trap_hit)
rnai_mat <- as.matrix(rnai_mat)
fit_rnai <- eulerr::euler(rnai_mat)
pdf(file.path("figures","pandependency_agreement_RNAi_euler_venn.pdf"), height=3.9,width=4)
plot(fit_rnai)
dev.off()
nrow(subset(t2,(RNAi_PD & support))) / nrow(subset(t2,RNAi_PD))

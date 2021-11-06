
data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")
source(file.path("src","id_utility.R"))

cors_list <- list("crispr"=fread(file.path(data_processed,"CRISPR-Avana_vs_CRISPR-KY_cor_recall_vals.csv")),
                  "rnai"=fread(file.path(data_processed,"RNAi-Achilles_vs_RNAi-DRIVE_cor_recall_vals.csv")))
cors_list[["crispr"]]$entrez_id %<>% as.character(.)
cors_list[["rnai"]]$entrez_id %<>% as.character(.)

gs_list <- list("crispr_ky" = fread(file.path(data_raw,"gene-effect-scaled-crispr-ky.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                "crispr_avana" = fread(file.path(data_raw,"gene-effect-scaled-crispr-avana.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                "rnai_achilles" = fread(file.path(data_processed,"gene-effect-scaled-rnai-achilles.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.),
                "rnai_drive" = fread(file.path(data_processed,"gene-effect-scaled-rnai-drive.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.))

# gene ID updates which were performed on RNAi datasets as part of target discovery pipeline
# hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
# hgnc$entrez_id %<>% as.character(.)
# hgnc %<>% subset(.,locus_group == "protein-coding gene")
# 
# ach_gs <- gs_list[["rnai_achilles"]]
# colnames(ach_gs) <- extract_entrez(colnames(ach_gs))
# ach_gs <- ach_gs[,colnames(ach_gs) %in% hgnc$entrez_id]
# colnames(ach_gs) <- entrez_to_cds(colnames(ach_gs),hgnc)
# gs_list[["rnai_achilles"]] <- ach_gs
# 
# drive_gs <- gs_list[["rnai_drive"]]
# colnames(drive_gs) <- extract_entrez(colnames(drive_gs))
# drive_gs <- drive_gs[,colnames(drive_gs) %in% hgnc$entrez_id]
# colnames(drive_gs) <- entrez_to_cds(colnames(drive_gs),hgnc)
# gs_list[["rnai_drive"]] <- drive_gs

pr_list <- list("crispr_ky" = fread(file.path(data_processed,"dependency-probability-crispr-ky.csv")),
                "crispr_avana"=fread(file.path(data_processed,"dependency-probability-crispr-avana.csv")),
                "rnai_achilles" = fread(file.path(data_processed,"dependency-probability-rnai-achilles.csv")),
                "rnai_drive"= fread(file.path(data_processed,"dependency-probability-rnai-drive.csv")))
pr_list <- lapply(pr_list,function(x){x %>% column_to_rownames(.,var="Row.name")})

tech_dsets <- list("crispr"=c("crispr_ky","crispr_avana"),
                   "rnai"=c("rnai_achilles","rnai_drive"))

# gene <- Reduce(intersect,lapply(gs_list,function(x){colnames(x)}))
# gs_list <- lapply(gs_list,function(x){x[,gene]})
# pr_list <- lapply(pr_list,function(x){x[,gene]})

res_df <- list()
for (dset in names(gs_list)){
  gs=gs_list[[dset]]
  pr=pr_list[[dset]]
  
  gs <- gs[,colnames(pr)]
  stopifnot(all(colnames(pr)==colnames(gs)))
  
  v=matrixStats::colVars(gs,na.rm=T)
  v_norm=v/max(v)
  
  df <- data.frame(gene=colnames(gs),
                   var=v,var_norm=v_norm,
                   frac=(colSums(pr>.5,na.rm=T)/colSums(!is.na(pr))),
                   stringsAsFactors = F)
  colnames(df)[2:4] <- paste0(dset,"_",colnames(df)[2:4]) 
  res_df[[dset]] <- df
  
}

crispr_df <- inner_join(res_df[["crispr_ky"]],res_df[["crispr_avana"]],by="gene")
crispr_df$entrez_id <- extract_entrez(crispr_df$gene)
crispr_df %<>% right_join(.,cors_list[["crispr"]],by="entrez_id")
crispr_df$mean_var_norm <- rowMeans(crispr_df[,c("crispr_avana_var_norm","crispr_ky_var_norm")])
crispr_df$mean_depFrac <- rowMeans(crispr_df[,c("crispr_avana_frac","crispr_ky_frac")])
crispr_df %<>% mutate(.,dataset="CRISPR") %>% dplyr::select(.,gene,mean_var_norm,mean_depFrac,r,dataset)

rnai_df <- inner_join(res_df[["rnai_achilles"]],res_df[["rnai_drive"]],by="gene")
rnai_df$entrez_id <- extract_entrez(rnai_df$gene)
rnai_df %<>% right_join(.,cors_list[["rnai"]],by="entrez_id")
rnai_df$mean_var_norm <- rowMeans(rnai_df[,c("rnai_achilles_var_norm","rnai_drive_var_norm")])
rnai_df$mean_depFrac <- rowMeans(rnai_df[,c("rnai_achilles_frac","rnai_drive_frac")])
rnai_df %<>% mutate(.,dataset="RNAi") %>% dplyr::select(.,gene,mean_var_norm,mean_depFrac,r,dataset)

plot_df <- bind_rows(crispr_df,rnai_df)
tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

ggplot(plot_df,aes(x=mean_depFrac,y=r,color=dataset)) +
  geom_smooth() +
  scale_color_manual(values=tech_pal) +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  ylab("Gene Effect Correlation") +
  xlab("Dependent Cell Line Fraction")
ggsave(file.path("figures","high_confidence_depFrac_vs_correlation.pdf"),width=2,height=2)

ggplot(plot_df,aes(x=mean_depFrac,y=mean_var_norm,color=dataset)) +
  geom_smooth() +
  scale_color_manual(values=tech_pal) +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  ylab("Gene effect variance") +
  xlab("Dependent Cell Line Fraction")
ggsave(file.path("figures","high_confidence_depFrac_vs_variance.pdf"),width=2,height=2)




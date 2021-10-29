
cors_list <- list("crispr"=fread("/Users/mburger/dynamic-duo-biorxiv/figures/library_agreement/CRISPR-Avana_vs_CRISPR-KY_cor_recall_vals.csv"),
                  "rnai"=fread("/Users/mburger/dynamic-duo-biorxiv/figures/library_agreement/RNAi-Achilles_vs_RNAi-DRIVE_cor_recall_vals.csv"))
cors_list[["crispr"]]$entrez_id %<>% as.character(.)
cors_list[["rnai"]]$entrez_id %<>% as.character(.)

gs_list <- list("crispr_ky" = load.from.taiga(data.name='inputs-ef6d', data.version=1, data.file='ceres_ky_scaled'),
                "crispr_avana" = load.from.taiga(data.name='inputs-ef6d', data.version=1, data.file='gene_effect_corrected'),
                "rnai_achilles" = load.from.taiga(data.name='inputs-7713', data.version=1, data.file='DEMETER2_Achilles_gene_effect'),
                "rnai_drive" = load.from.taiga(data.name='inputs-7713', data.version=1, data.file='DEMETER2_DRIVE_gene_effect'))

pr_list <- list("crispr_ky" = fread("/Users/mburger/Data/dynamic-duo/Avana_KY_individual/KY-PR.csv"),
                "crispr_avana"=fread("/Users/mburger/Data/dynamic-duo/Avana_KY_individual/Avana-PR.csv"),
                "rnai_achilles" = fread("/Users/mburger/Data/dynamic-duo/Achilles_DRIVE_individual/Achilles-PR.csv"),
                "rnai_drive"= fread("/Users/mburger/Data/dynamic-duo/Achilles_DRIVE_individual/DRIVE-PR.csv"))
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
  
  v=colVars(gs,na.rm=T)
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
ggplot(plot_df,aes(x=r,y=mean_var_norm,color=dataset)) +
  geom_smooth() +
  scale_color_manual(values=tech_pal) +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  ylab("Mean Normalized Variance") +
  xlab("Pearson Correlation") +
  coord_flip()
ggsave("/Users/mburger/dynamic-duo-biorxiv/figures/library_agreement/variance_vs_correlation.pdf",width=2,height=2)

ggplot(plot_df,aes(x=mean_depFrac,y=r,color=dataset)) +
  geom_smooth() +
  scale_color_manual(values=tech_pal) +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  ylab("Pearson Correlation") +
  xlab("Mean Dep. Cell Line Fraction")
ggsave("/Users/mburger/dynamic-duo-biorxiv/figures/library_agreement/depFrac_vs_correlation.pdf",width=2,height=2)

ggplot(plot_df,aes(x=mean_depFrac,y=mean_var_norm,color=dataset)) +
  geom_smooth() +
  scale_color_manual(values=tech_pal) +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  ylab("Gene effect variance") +
  xlab("Mean Dep. Cell Line Fraction")
ggsave("/Users/mburger/dynamic-duo-biorxiv/figures/library_agreement/depFrac_vs_variance.pdf",width=2,height=2)





# require(pROC)

source("src/packages_paths.R")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)

#Get genes that both CRISPR datasets agree are pan-dependent or agree are not pan-dependent
t1 <- fread(file.path(data_processed,"library_agreement.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,(CRISPR_avana_target & CRISPR_ky_target))
t1 %<>% subset(., (CRISPR_pandep == 2) | (CRISPR_pandep == 0))

keep_genes <- entrez_to_cds(t1$entrez_id,hgnc)
conf_ce <- subset(t1,CRISPR_pandep == 2) %>% select(., symbol, entrez_id)
conf_ce %<>% add_column(.,CDS_ID=entrez_to_cds(conf_ce$entrez_id,hgnc),.before=1)

#Get RNAi gene effect mean 
rnai_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix()
rnai_means <- colMeans(rnai_gs,na.rm=T)
rnai_means <- data.frame(Gene=names(rnai_means),entrez_id=extract_entrez(names(rnai_means)),value=rnai_means,stringsAsFactors = F)
rnai_means %<>% subset(.,Gene %in% keep_genes)

#Determine how accurately you can predict CRISPR pandep status using each technology
eg_summary <- fread(file.path(data_raw,"control-essential-genes-unbiased.csv"))
eg_summary$entrez_id %<>% as.character(.)

eg_summary %<>% dplyr::select(.,entrez_id,ExAC=human_median,`Gene Trap`=gene_trap_percentile_mean)
rnai_means %<>% dplyr::select(.,entrez_id,RNAi=value) %>% remove_rownames(.)
ranks <- left_join(rnai_means,eg_summary,by="entrez_id")

ranks$class <- ranks$entrez_id %in% conf_ce$entrez_id

r_list <- list()
# for (dset in c("RNAi","ExAC","Gene Trap","DRIVE","Multilib")){
for (dset in c("RNAi","ExAC","Gene Trap")){
  r <- pROC::roc(ranks$class, ranks[,dset])
  r_list[[dset]] <- data.frame(dataset=dset,
                               specificity=r$specificities,
                               sensitivity=r$sensitivities,
                               thresholds=r$thresholds,
                               auc=as.numeric(r$auc),
                               stringsAsFactors = F)
}

df_r <- bind_rows(r_list)

mypal <- c("Gene Trap"="#F16667","RNAi"="#996699","ExAC"="#653614")

ggplot(subset(df_r,dataset %in% c("RNAi","ExAC","Gene Trap")), aes(x=specificity, y=sensitivity,color=dataset)) +
  geom_line(size=1) +
  geom_abline(intercept = 1,slope=1,col="gray",linetype="dashed") +
  scale_color_manual(values=mypal) +
  scale_x_reverse() +
  ylab("sensitivity") +
  xlab("specificity") +
  theme_classic(base_size=11) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.75,.25))
ggsave(file.path("figures","pandependency_agreement_ROC_predicting_CRISPR_class.pdf"),height=2.5,width=2.1)

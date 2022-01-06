
source("src/packages_paths.R")

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix()

mut_hot <- fread(file.path(data_raw,"depmap-omics-mutation-hotspot.csv")) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)
mut_hot %<>% as.data.frame(.) %>% rownames_to_column(.,var="DepMap_ID")
mut_hot %<>% dplyr::select(.,DepMap_ID,`BRAF (673)`)
colnames(mut_hot)[2] <- "BRAF hotspot"

crispr_gs %<>% as.data.frame(.) %>% rownames_to_column(.,var="DepMap_ID")
crispr_gs %<>% dplyr::select(.,DepMap_ID,`BRAF (673)`)

crispr_gs %<>% left_join(.,mut_hot,by="DepMap_ID")
crispr_gs$`BRAF hotspot` <- crispr_gs$`BRAF hotspot` == 1

quick_pal <- c("FALSE"="grey","TRUE"="#4B68B1")
ggplot(crispr_gs,aes(x=`BRAF (673)`,fill=`BRAF hotspot`)) +
  geom_density(alpha=1) +
  geom_vline(xintercept = 0,linetype = "dashed",color="#625441",size=1) +
  geom_vline(xintercept = -1,linetype = "dashed",color="#C15327",size=1) +
  theme_classic(base_size=11) +
  theme(legend.position = c(.3,.8)) +
  scale_fill_manual(values=quick_pal) +
  xlab("CRISPR gene effect")
ggsave(file.path("figures","gene_deps_profiles_BRAF_hotspot_density.pdf"),width=3,height=2.5)


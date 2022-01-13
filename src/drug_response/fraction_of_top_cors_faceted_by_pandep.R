
source("src/packages_paths.R")

#Get the best correlated dose for each drug-gene target pair
univariate_res <- fread(file.path(data_processed,"drug-screen-genetic-targets-correlations.csv"))
plot_df <- subset(univariate_res,best_drug_target_cor)

#Annotate the gene targets that are CRISPR pan-dependencies
t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
pandeps <- subset(t2,CRISPR_PD)$entrez_id

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
pandeps <- entrez_to_cds(pandeps,hgnc)
plot_df$pandep <- plot_df$target %in% pandeps

#total pairs, iterate across drug datasets and pandep status
totals <- list()
for (dset in c("ctd2","prism","gdsc")){
  tmp_df <- subset(plot_df,(drug_dataset == dset) & (genetic_perturbation == "CRISPR"))
  pandep_pairs <- nrow(subset(tmp_df,pandep))
  other_pairs <- nrow(subset(tmp_df,!pandep))
  
  hits_df <- subset(plot_df,(rank <= 5) & (drug_dataset == dset))
  crispr_pandep_pairs <- nrow(subset(hits_df ,pandep & (genetic_perturbation == "CRISPR")))
  crispr_other_pairs <- nrow(subset(hits_df ,(!pandep) & (genetic_perturbation == "CRISPR")))
  
  rnai_pandep_pairs <- nrow(subset(hits_df ,pandep & (genetic_perturbation == "RNAi")))
  rnai_other_pairs <- nrow(subset(hits_df ,(!pandep) & (genetic_perturbation == "RNAi")))
  
  tmp_res <- data.frame(perturbation=c("CRISPR","CRISPR","RNAi","RNAi"),
                        pandep=c(T,F,T,F),
                        total=c(pandep_pairs,other_pairs,pandep_pairs,other_pairs),
                        hits=c(crispr_pandep_pairs,crispr_other_pairs,rnai_pandep_pairs,rnai_other_pairs),
                        stringsAsFactors = F)
  tmp_res$frac <- tmp_res$hits / tmp_res$total
  tmp_res$drug_dataset <- dset
  
  totals[[dset]] <- tmp_res
}

hit_frac <- bind_rows(totals)

hit_frac$pandep %<>% as.character(.)
hit_frac$pandep <- factor(hit_frac$pandep,levels=c("TRUE","FALSE"))
tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")
ggplot(hit_frac,aes(x=drug_dataset,y=frac,fill=perturbation)) +
  geom_bar(stat="identity",position="dodge",width=.8) +
  scale_fill_manual(values=tech_pal) +
  facet_grid(. ~ pandep, scales = "free") +
  ylab("Frac. of drug-gene target pairs\n(gene is top 5 correlate of drug)") +
  xlab("") +
  theme_bw(base_size=11) +
  theme(legend.position = "none")
ggsave(file.path("figures","drug_response_fraction_of_top_cors_faceted_by_pandep.pdf"),width=3.5,height=2.3)

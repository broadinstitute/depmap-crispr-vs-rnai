
source("src/packages_paths.R")

standard_doses <- list("ctd2"=16,
                      "prism"=8)

univariate_res <- fread(file.path(data_processed,"drug-screen-genetic-targets-correlations.csv"))

for (dset in unique(univariate_res$drug_dataset)){
  moa <- subset(univariate_res,drug_dataset == dset)
  
  #Get drug-target pairs where the specific target is a top 5 correlate to the drug_pert in either CRISPR or RNAi
  success_pairs <- unique(subset(moa,rank <= 5)$target_drug_pair)
  print(paste0(dset,": ",length(success_pairs)," drug-gene target pairs"))
  
  moa %<>% subset(.,target_drug_pair %in% success_pairs)
  
  #For the filtered pairs, plot the pearson for the best dose of each drug-target pair
  best_cor <- subset(moa,best_drug_target_cor)

  tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")
  ggplot(best_cor,aes(x=genetic_perturbation,y=pearson_r,color=genetic_perturbation)) +
    geom_boxplot() +
    scale_color_manual(values=tech_pal) +
    theme_bw(base_size=11) +
    theme(legend.position = "none") +
    xlab("")
  ggsave(file.path("figures",paste0("drug_response_",dset,"_best_cor_boxplot.pdf")),width=1.5,height=2.5)
}
  
for (dset in names(standard_doses)){

  #Per dose plots 
  moa <- subset(univariate_res,drug_dataset == dset)
  
  moa$dose <- as.numeric(moa$dose)
  
  #Remove drugs with abnormal numbers of doses
  blacklist <- unique(subset(moa,dose > standard_doses[[dset]])$broad_id)
  moa %<>% subset(.,!(broad_id %in% blacklist))
  
  print(paste0(dset,": ",length(unique(moa$broad_id))," drugs using standard dose range"))
  
  moa$dose <- as.character(moa$dose)
  #take best target per drug at each dose
  best_target <- list()
  for (tech in unique(moa$genetic_perturbation)){
    for (dose_name in unique(moa$dose)){
      tmp_res <- subset(moa,(genetic_perturbation == tech) & (dose == dose_name))
      tmp_res <- tmp_res[order(tmp_res$pearson_r,decreasing=T),]
      tmp_res$best_target_cor <- !duplicated(tmp_res$broad_id)
      
      tmp_res <- tmp_res[order(tmp_res$rank,decreasing=F),]
      tmp_res$best_target_rank <- !duplicated(tmp_res$broad_id)
      best_target[[paste0(dose_name,"-",tech)]] <- tmp_res
    }
  }
  best_target <- bind_rows(best_target)
  
  best_target$dose <- factor(best_target$dose,levels=as.character(1:standard_doses[[dset]]))
  
  tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")
  ggplot(subset(best_target, best_target_rank & (rank <= 5)),aes(x=dose,fill=genetic_perturbation)) +
    geom_bar(stat="count",position="dodge",width=.8) +
    scale_fill_manual(values=tech_pal) +
    theme_bw(base_size=11) +
    xlab("Dose") +
    theme(legend.position = "none")
  ggsave(file.path("figures",paste0("drug_response_",dset,"_drugs_with_top5_targets.pdf")),width=2.5,height=2.5)
  
  
}


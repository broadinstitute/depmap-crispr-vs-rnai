
library(stringr)

source("src/packages_paths.R")

standard_doses <- list("ctd2"=16,
                       "prism"=8)

univariate_res <- fread(file.path(data_processed,"drug-screen-genetic-targets-correlations.csv"))

for (dset in names(standard_doses)){
  
  #Per dose plots 
  summary <- subset(univariate_res,drug_dataset == dset)
  
  summary$dose <- as.numeric(summary$dose)
  blacklist <- unique(subset(summary,dose > standard_doses[[dset]])$broad_id)
  summary %<>% subset(.,!(broad_id %in% blacklist))
  summary %<>% subset(.,best_drug_target_cor)
  
  #Keep drug-gene pairs that are top 5 using either CRISPR or RNAi
  keep_pairs <- unique(subset(summary,rank <= 5)$target_drug_pair)
  print(paste0(dset,": ",length(keep_pairs)," drug-gene target pairs"))
  
  crispr_summary <- subset(summary,genetic_perturbation == "CRISPR") %>% 
    subset(.,target_drug_pair %in% keep_pairs) %>% 
    dplyr::select(.,target_drug_pair,CRISPR_dose=dose)
  rnai_summary <- subset(summary,genetic_perturbation == "RNAi") %>% 
    subset(.,target_drug_pair %in% keep_pairs) %>% 
    dplyr::select(.,target_drug_pair,RNAi_dose=dose)
  wide_df <- inner_join(crispr_summary,rnai_summary,by="target_drug_pair")
  
  #Count the doses 
  level_vec <- sort(unique(c(wide_df$CRISPR_dose,wide_df$RNAi_dose)))
  count_mat <- matrix(rep(0,(length(level_vec))^2),nrow=length(level_vec),ncol=length(level_vec))
  colnames(count_mat) <- as.character(level_vec)
  rownames(count_mat) <- as.character(level_vec)
  
  for (rnai_dose in level_vec){
    
    tmp_df <- subset(wide_df,RNAi_dose == rnai_dose)
    
    for (crispr_dose in level_vec){
      
      count_mat[as.character(rnai_dose),as.character(crispr_dose)] <- sum(tmp_df$CRISPR_dose == crispr_dose)
      
    }
  }
  
  longData<-melt(count_mat)
  
  longData$Var1 <- factor(longData$Var1,levels=as.character(level_vec))
  longData$Var2 <- factor(longData$Var2,levels=as.character(level_vec))
  
  pal_func <- function(n){
    
    gradpal <- colorRampPalette(c("grey95", "red"))(n-1)
    final_pal <- c("white",gradpal)
    names(final_pal) <- 0:(max(count_mat))
    return(final_pal)
    return(gradpal)
  }
  
  longData$value <- factor(longData$value)
  
  p <- ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    geom_text(aes(label=value)) +
    scale_fill_manual(values = pal_func(max(count_mat)+1)) +
    labs(x="CRISPR Top Correlated Dose", y="RNAi Top Correlated Dose") +
    theme_classic(base_size=11) 

  ggsave(plot=p+theme(legend.position = "none"),file.path("figures",paste0("drug_response_",dset,"_top_dose_heatmap.pdf")),width=2.5,height=2.5)
  
  legend <- cowplot::get_legend(p)
  grid::grid.newpage()
  
  pdf(file.path("figures",paste0("drug_response_",dset,"_top_dose_heatmap_legend.pdf")))
  grid::grid.draw(legend) 
  dev.off()
}

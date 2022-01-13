
source("src/packages_paths.R")

#Matched genetic screens
genetic_data <- list("crispr"="gene-effect-scaled-crispr-matched.csv",
                     "rnai"="gene-effect-scaled-rnai-matched.csv"
)
genetic_data <- lapply(genetic_data,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

#Dose-level drug data
drug_data <- list(
  # "gdsc"="drug-screen-viability-filtered-noPC1-gdsc.csv",
  # "prism"="drug-screen-viability-filtered-noPC1-prism.csv",
  "ctd2"="drug-screen-viability-filtered-noPC1-ctd2.csv"
)
drug_data <- lapply(drug_data,function(x){fread(file.path(data_processed,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

dose_info <- list(
  # "prism"="drug-screen-viability-info-prism.csv",
  # "gdsc"="drug-screen-viability-info-gdsc.csv",
  "ctd2"="drug-screen-viability-info-ctd2.csv"
)
dose_info <- lapply(dose_info,function(x){fread(file.path(data_raw,x))})

#Filter for shared cell lines
cls <- intersect(rownames(genetic_data[["crispr"]]),rownames(drug_data[["ctd2"]]))
crispr <- genetic_data[["crispr"]][cls,]
rnai <- genetic_data[["rnai"]][cls,]
drugs <- drug_data[["ctd2"]][cls,]

dose_info <- dose_info[["ctd2"]]
dose_info$factor_dose %<>% as.character(.)

rephub_annot <- fread(file.path(data_raw,"drug-screen-target-annotations-rephub.csv"))

#Getting doses from main scatter plot plot_df

example_cors <- list("BCL2L1"=list(drug_id = "BRD-K82746043",gene = "BCL2L1 (598)",tmp_dose = "13"),
                     "CHEK1"=list(drug_id = "BRD-K86525559",gene = "CHEK1 (1111)",tmp_dose = "6"),
                     "WEE1"=list(drug_id = "BRD-K54256913",gene = "WEE1 (7465)",tmp_dose = "11"),
                     "CDK2"=list(drug_id = "BRD-K43389698",gene = "CDK2 (1017)",tmp_dose = "3"))

for (p in names(example_cors)){
  
  drug_id=example_cors[[p]][["drug_id"]]
  gene=example_cors[[p]][["gene"]]
  tmp_dose=example_cors[[p]][["tmp_dose"]]
  
  drug_name = unique(subset(rephub_annot,broad_id == drug_id)$name)
  gene_name = gsub(" .*","",gene)
  dose_conc = unique(subset(dose_info,(broad_id == drug_id) & (factor_dose == tmp_dose))$dose)
  crispr_xy <- data.frame(GENE=crispr[,gene],DRUG=drugs[,paste0(drug_id,"::",tmp_dose)],pert="CRISPR")
  rnai_xy <- data.frame(GENE=rnai[,gene],DRUG=drugs[,paste0(drug_id,"::",tmp_dose)],pert="RNAi")
  xy <- bind_rows(crispr_xy,rnai_xy)
  
  tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")
  ggplot(xy,aes(x=DRUG,y=GENE,color=pert)) +
    geom_density2d(size=.5) +
    geom_smooth(method="lm") +
    theme_bw(base_size=11) +
    scale_color_manual(values=tech_pal) +
    theme(legend.position = "none") +
    xlab(paste0(drug_name, " (",dose_conc,")")) +
    ylab(paste0(gene_name," gene effect"))
  ggsave(file.path("figures",paste0("drug_response_CTD2_",drug_name,"_",gene_name,".pdf")),width=1.5,height=1.5)
  
  
}



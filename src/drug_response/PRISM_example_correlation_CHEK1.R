
source("src/packages_paths.R")

#Matched genetic screens
genetic_data <- list("crispr"="gene-effect-scaled-crispr-matched.csv",
                     "rnai"="gene-effect-scaled-rnai-matched.csv"
)
genetic_data <- lapply(genetic_data,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

#Dose-level drug data
drug_data <- list(
  # "gdsc"="drug-screen-viability-filtered-noPC1-gdsc.csv",
  "prism"="drug-screen-viability-filtered-noPC1-prism.csv"
  # "ctd2"="drug-screen-viability-filtered-noPC1-ctd2.csv"
)
drug_data <- lapply(drug_data,function(x){fread(file.path(data_processed,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

dose_info <- list(
  "prism"="drug-screen-viability-info-prism.csv"
  # "gdsc"="drug-screen-viability-info-gdsc.csv",
  # "ctd2"="drug-screen-viability-info-ctd2.csv"
)
dose_info <- lapply(dose_info,function(x){fread(file.path(data_raw,x))})

#broad_id: BRD-K75009076-001-02-1
#name: SCH-900776
#dose: 0.607250000

dose_info <- dose_info[["prism"]]
drug_pert <- subset(dose_info,(broad_id == "BRD-K75009076") & (dose == 0.60725))$column_name

#Filter compound data for launched oncology compounds
drug_data <- drug_data[["prism"]]

crispr <- genetic_data[["crispr"]]
rnai <- genetic_data[["rnai"]]

#Filter for shared cell lines
cls <- intersect(rownames(drug_data),rownames(crispr))
crispr <- crispr[cls,]
rnai <- rnai[cls,]
drug_data <- drug_data[cls,]

df1 <- data.frame(gene=crispr[,"CHEK1 (1111)"],drug=drug_data[,drug_pert],tech="CRISPR",dose="0.60725",stringsAsFactors = F)
df2 <- data.frame(gene=rnai[,"CHEK1 (1111)"],drug=drug_data[,drug_pert],tech="RNAi",dose="0.60725",stringsAsFactors = F)

plot_df <- bind_rows(df1,df2)

tech_pal <- c("CRISPR"="#0099B4FF","RNAi"="#925E9FFF")

ggplot(plot_df,aes(x=drug,y=gene,color=tech)) +
  geom_density_2d() +
  geom_smooth(method="lm") +
  scale_color_manual(values=tech_pal) +
  theme_classic(base_size=11) +
  theme(legend.position = c(.2,.8)) +
  ylab("CHEK1 Gene Effect") +
  xlab("CHK inhibitor viability\nSCH-900776 (0.60725)")
ggsave(file.path("figures","drug_response_PRISM_CHKinhibitor_CHEK1.pdf"),width=2,height=2.5)

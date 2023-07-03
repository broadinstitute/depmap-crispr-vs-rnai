library(tidyverse)
library(magrittr)
library(data.table)

ce_genes <- c("RBX1", "USP39","RNF125")

# CRISPR
sgrna_map <- fread("data/raw/reagent-to-gene-map-sgrna.csv")
sgrna_map %<>% subset(.,symbol %in% ce_genes)

# Sanger
ky_lfc <- fread("data/raw/lfc-unscaled-crispr-ky.csv") %>% column_to_rownames(.,var="V1")
ky_lfc <- scale(ky_lfc)
ky_lfc <- ky_lfc[rownames(ky_lfc) %in% sgrna_map$reagent,]
ky_u <- rowMeans(ky_lfc)
ky_mean <- data.frame(reagent=names(ky_u),reagent_mean_zscore_lfc=ky_u) %>% remove_rownames(.)
ky_mean %<>% left_join(.,sgrna_map,by="reagent")
ky_mean %<>% dplyr::select(.,reagent,reagent_mean_zscore_lfc,entrez_id,symbol)
ky_mean$library <- "KY"
ky_mean$type <- "CRISPR"

# Broad
avana_lfc <- fread("data/raw/lfc-unscaled-crispr-avana.csv") %>% column_to_rownames(.,var="V1")
avana_lfc <- scale(avana_lfc)
avana_lfc <- avana_lfc[rownames(avana_lfc) %in% sgrna_map$reagent,]
avana_u <- rowMeans(avana_lfc)
avana_mean <- data.frame(reagent=names(avana_u),reagent_mean_zscore_lfc=avana_u) %>% remove_rownames(.)
avana_mean %<>% left_join(.,sgrna_map,by="reagent")
avana_mean %<>% dplyr::select(.,reagent,reagent_mean_zscore_lfc,entrez_id,symbol)
avana_mean$library <- "Avana"
avana_mean$type <- "CRISPR"

# RNAi
shrna_map <- fread("data/raw/reagent-to-gene-map-shrna.csv")
shrna_map %<>% subset(.,symbol %in% ce_genes)

# DRIVE
drive_lfc <- fread("data/raw/lfc-unscaled-rnai-drive.csv") %>% column_to_rownames(.,var="V1")
drive_lfc <- scale(drive_lfc)
drive_lfc <- drive_lfc[rownames(drive_lfc) %in% shrna_map$reagent,]
drive_u <- rowMeans(drive_lfc,na.rm=T)
drive_mean <- data.frame(reagent=names(drive_u),reagent_mean_zscore_lfc=drive_u) %>% remove_rownames(.)
drive_mean %<>% left_join(.,shrna_map,by="reagent")
drive_mean %<>% dplyr::select(.,reagent,reagent_mean_zscore_lfc,entrez_id,symbol)
drive_mean$library <- "DRIVE"
drive_mean$type <- "RNAi"

# Achilles
achilles_lfc <- fread("data/raw/lfc-unscaled-rnai-achilles.csv") %>% column_to_rownames(.,var="V1")
achilles_lfc <- scale(achilles_lfc)
achilles_lfc <- achilles_lfc[rownames(achilles_lfc) %in% shrna_map$reagent,]
achilles_u <- rowMeans(achilles_lfc,na.rm=T)
achilles_mean <- data.frame(reagent=names(achilles_u),reagent_mean_zscore_lfc=achilles_u) %>% remove_rownames(.)
achilles_mean %<>% left_join(.,shrna_map,by="reagent")
achilles_mean %<>% dplyr::select(.,reagent,reagent_mean_zscore_lfc,entrez_id,symbol)
achilles_mean$library <- "Achilles"
achilles_mean$type <- "RNAi"

plot_df <- bind_rows(ky_mean,avana_mean,drive_mean,achilles_mean)

lib_pal <- c("KY"="#E64B35FF","Avana"="#F39B7FFF","Achilles"="#3C5488FF","DRIVE"="#4DBBD5FF")

plot_df_2 <- subset(plot_df,symbol == "RNF125")
plot_df <- subset(plot_df,symbol != "RNF125" )

ggplot(plot_df,aes(x=type,y=reagent_mean_zscore_lfc,color=library)) +
  geom_point(position=position_jitterdodge(jitter.width = .15),size=2,alpha=.7) +
  geom_boxplot(outlier.shape=NA,fill=NA) +
  facet_grid(. ~ symbol, scales = "free", space = "free") +
  theme_bw(base_size=16) +
  xlab("") +
  scale_color_manual(values=lib_pal) +
  ylab("Reagent mean log fold-change") +
  theme(legend.position = "bottom")
ggsave("figures/efficacy_specificity_RBX1_USP39_reagents.pdf",
       height=5,width=6)

plot2_crispr <- subset(plot_df_2,type == "CRISPR")
plot2_rnai <- subset(plot_df_2,type == "RNAi")

shared_reagents <- plot2_rnai$reagent[duplicated(plot2_rnai$reagent)]
plot2_rnai_shared <- subset(plot2_rnai,reagent %in% shared_reagents)
plot2_rnai_shared$type <- "RNAi (shared)"
plot2_rnai_unique <- subset(plot2_rnai,!(reagent %in% shared_reagents))
plot2_rnai_unique$type <- "RNAi (unique)"

plot2 <- bind_rows(plot2_crispr,plot2_rnai_shared,plot2_rnai_unique)

ggplot(plot2,aes(x=type,y=reagent_mean_zscore_lfc,color=library)) +
  geom_point(position=position_jitterdodge(jitter.width = .15),size=2,alpha=.7) +
  geom_boxplot(outlier.shape=NA,fill=NA) +
  facet_grid(. ~ symbol, scales = "free", space = "free") +
  theme_bw(base_size=16) +
  xlab("") +
  scale_color_manual(values=lib_pal) +
  ylab("Reagent mean log fold-change") +
  theme(legend.position = "bottom")
ggsave("figures/efficacy_specificity_RNF125_reagents.pdf",
       height=5,width=4)

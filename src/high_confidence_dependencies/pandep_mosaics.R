
source("src/packages_paths.R")

library(ggmosaic)

#### CRISPR pan-dependency

t3 <- fread(file.path(data_processed,"library_agreement.csv"))
t3$entrez_id %<>% as.character(.)

t3 %<>% subset(.,CRISPR_avana_target & CRISPR_ky_target)

crispr_pandep <- subset(t3,CRISPR_avana_pandep & CRISPR_ky_pandep)$entrez_id
avana_only_pandep <- subset(t3,CRISPR_avana_pandep & !CRISPR_ky_pandep)$entrez_id
ky_only_pandep <- subset(t3,!CRISPR_avana_pandep & CRISPR_ky_pandep)$entrez_id
crispr_not_pandep <- subset(t3,!CRISPR_avana_pandep & !CRISPR_ky_pandep)$entrez_id

a <- length(crispr_pandep)
b <- length(avana_only_pandep)
c <- length(ky_only_pandep)
d <- length(crispr_not_pandep)

plot_df <- select(t3,entrez_id,CRISPR_avana_pandep,CRISPR_ky_pandep) %>% column_to_rownames(.,var = "entrez_id")

for (dset in colnames(plot_df)){
  
  x_name <- rep("No",nrow(plot_df))
  x_name[plot_df[,dset]] <- "Yes"
  plot_df[,dset] <- x_name
  
}

plot_df$CRISPR_ky_pandep <- factor(plot_df$CRISPR_ky_pandep,levels=c("Yes","No"))
crispr_mosiac <- ggplot(data = plot_df) +
  geom_mosaic(aes(x = product(CRISPR_avana_pandep, CRISPR_ky_pandep)), na.rm=TRUE) +
  annotate("text", x = c(.2,.8,.2,.8), y = c(.8,.8,.2,.2), label =	as.character(c(a,b,c,d)),size=5) +
  xlab("KY") +
  ylab("Avana") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
ggsave(plot=crispr_mosiac,file.path("figures","high_confidence_crispr_pandep_mosaic.pdf"),width=1.25,height=1.25)

#### RNAi pan-dependency

t3 <- fread(file.path(data_processed,"library_agreement.csv"))
t3$entrez_id %<>% as.character(.)

t3 %<>% subset(.,RNAi_achilles_target & RNAi_drive_target)

rnai_pandep <- subset(t3,RNAi_achilles_pandep & RNAi_drive_pandep)$entrez_id
achilles_only_pandep <- subset(t3,RNAi_achilles_pandep & !RNAi_drive_pandep)$entrez_id
drive_only_pandep <- subset(t3,!RNAi_achilles_pandep & RNAi_drive_pandep)$entrez_id
rnai_not_pandep <- subset(t3,!RNAi_achilles_pandep & !RNAi_drive_pandep)$entrez_id

a <- length(rnai_pandep)
b <- length(achilles_only_pandep)
c <- length(drive_only_pandep)
d <- length(rnai_not_pandep)

plot_df <- select(t3,entrez_id,RNAi_achilles_pandep,RNAi_drive_pandep) %>% column_to_rownames(.,var = "entrez_id")

for (dset in colnames(plot_df)){
  
  x_name <- rep("No",nrow(plot_df))
  x_name[plot_df[,dset]] <- "Yes"
  plot_df[,dset] <- x_name
  
}

plot_df$RNAi_drive_pandep <- factor(plot_df$RNAi_drive_pandep,levels=c("Yes","No"))
rnai_mosiac <- ggplot(data = plot_df) +
  geom_mosaic(aes(x = product(RNAi_achilles_pandep, RNAi_drive_pandep)), na.rm=TRUE) +
  annotate("text", x = c(.2,.8,.2,.8), y = c(.8,.8,.2,.2), label =	as.character(c(a,b,c,d)),size=5) +
  xlab("DRIVE") +
  ylab("Achilles") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
ggsave(plot=rnai_mosiac,file.path("figures","high_confidence_rnai_pandep_mosaic.pdf"),width=1.25,height=1.25)


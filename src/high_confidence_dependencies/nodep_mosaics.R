
library(ggmosaic)

data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")

#### CRISPR non-dependency

t3 <- fread(file.path(data_processed,"library_agreement.csv"))
t3$entrez_id %<>% as.character(.)

t3 %<>% subset(.,CRISPR_avana_target & CRISPR_ky_target)

crispr_nodep <- subset(t3,CRISPR_avana_nodep & CRISPR_ky_nodep)$entrez_id
avana_only_nodep <- subset(t3,CRISPR_avana_nodep & !CRISPR_ky_nodep)$entrez_id
ky_only_nodep <- subset(t3,!CRISPR_avana_nodep & CRISPR_ky_nodep)$entrez_id
crispr_not_nodep <- subset(t3,!CRISPR_avana_nodep & !CRISPR_ky_nodep)$entrez_id

a <- length(crispr_nodep)
b <- length(avana_only_nodep)
c <- length(ky_only_nodep)
d <- length(crispr_not_nodep)

plot_df <- select(t3,entrez_id,CRISPR_avana_nodep,CRISPR_ky_nodep) %>% column_to_rownames(.,var = "entrez_id")

for (dset in colnames(plot_df)){
  
  x_name <- rep("No",nrow(plot_df))
  x_name[plot_df[,dset]] <- "Yes"
  plot_df[,dset] <- x_name
  
}

plot_df$CRISPR_ky_nodep <- factor(plot_df$CRISPR_ky_nodep,levels=c("Yes","No"))
crispr_mosiac <- ggplot(data = plot_df) +
  geom_mosaic(aes(x = product(CRISPR_avana_nodep, CRISPR_ky_nodep)), na.rm=TRUE) +
  annotate("text", x = c(.2,.8,.2,.8), y = c(.8,.8,.2,.2), label =	as.character(c(a,b,c,d)),size=5) +
  xlab("KY") +
  ylab("Avana") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
ggsave(plot=crispr_mosiac,file.path("figures","high_confidence_crispr_nodep_mosaic.pdf"),width=1.25,height=1.25)

#### RNAi non-dependency

t3 <- fread(file.path(data_processed,"library_agreement.csv"))
t3$entrez_id %<>% as.character(.)

t3 %<>% subset(.,RNAi_achilles_target & RNAi_drive_target)

rnai_nodep <- subset(t3,RNAi_achilles_nodep & RNAi_drive_nodep)$entrez_id
achilles_only_nodep <- subset(t3,RNAi_achilles_nodep & !RNAi_drive_nodep)$entrez_id
drive_only_nodep <- subset(t3,!RNAi_achilles_nodep & RNAi_drive_nodep)$entrez_id
rnai_not_nodep <- subset(t3,!RNAi_achilles_nodep & !RNAi_drive_nodep)$entrez_id

a <- length(rnai_nodep)
b <- length(achilles_only_nodep)
c <- length(drive_only_nodep)
d <- length(rnai_not_nodep)

plot_df <- select(t3,entrez_id,RNAi_achilles_nodep,RNAi_drive_nodep) %>% column_to_rownames(.,var = "entrez_id")

for (dset in colnames(plot_df)){
  
  x_name <- rep("No",nrow(plot_df))
  x_name[plot_df[,dset]] <- "Yes"
  plot_df[,dset] <- x_name
  
}

plot_df$RNAi_drive_nodep <- factor(plot_df$RNAi_drive_nodep,levels=c("Yes","No"))
rnai_mosiac <- ggplot(data = plot_df) +
  geom_mosaic(aes(x = product(RNAi_achilles_nodep, RNAi_drive_nodep)), na.rm=TRUE) +
  annotate("text", x = c(.2,.8,.2,.8), y = c(.8,.8,.2,.2), label =	as.character(c(a,b,c,d)),size=5) +
  xlab("DRIVE") +
  ylab("Achilles") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
ggsave(plot=rnai_mosiac,file.path("figures","high_confidence_rnai_nodep_mosaic.pdf"),width=1.25,height=1.25)


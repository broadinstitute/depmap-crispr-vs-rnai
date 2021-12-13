

# library(purrr)
# library(ggbeeswarm)
# library(ggjoy)
# library(viridis)
# library(ggsci)
# library(scales)
# require(cowplot)
# library(reshape2)
# library(grid)

source("src/packages_paths.R")

#Load CRISPR and RNAi gene effects
ceres <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix()
d2 <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix()

colnames(ceres) <- extract_entrez(colnames(ceres))
colnames(d2) <- extract_entrez(colnames(d2))

master <- fread(file.path("tables","Supplemental-Table-2.csv"))
master$entrez_id %<>% as.character(.)
  
gene_lists <- list()
#ND
gene_lists[["CRISPR_ND"]] <- subset(master,CRISPR_ND)$entrez_id
gene_lists[["RNAi_ND"]] <- subset(master,RNAi_ND)$entrez_id
  
#PD
gene_lists[["CRISPR_PD"]] <- subset(master,CRISPR_PD)$entrez_id
gene_lists[["RNAi_PD"]] <- subset(master,RNAi_PD)$entrez_id
  
#HVD
gene_lists[["CRISPR_HVD"]] <- subset(master,CRISPR_HVD)$entrez_id
gene_lists[["RNAi_HVD"]] <- subset(master,RNAi_HVD)$entrez_id
  
#SSD
gene_lists[["CRISPR_SSD"]] <- subset(master,CRISPR_SSD)$entrez_id
gene_lists[["RNAi_SSD"]] <- subset(master,RNAi_SSD)$entrez_id
  
#SD
all_crispr <- unique(c(gene_lists[["CRISPR_ND"]],gene_lists[["CRISPR_PD"]],gene_lists[["CRISPR_SSD"]],gene_lists[["CRISPR_HVD"]]))
all_rnai <- unique(c(gene_lists[["RNAi_ND"]],gene_lists[["RNAi_PD"]],gene_lists[["RNAi_SSD"]],gene_lists[["RNAi_HVD"]]))
gene_lists[["CRISPR_WSD"]] <- subset(master,!(entrez_id %in% all_crispr))$entrez_id
gene_lists[["RNAi_WSD"]] <- subset(master,!(entrez_id %in% all_rnai))$entrez_id

get_data_long <- function(data_mats,ceres,d2,sample_size){
  for (i in seq_along(data_mats)){
    indiv_list <- names(data_mats)[i]
    if (grepl("CRISPR",indiv_list)){
      filtered_mat <- ceres[,colnames(ceres) %in% data_mats[[i]],drop=F]
    } else {
      filtered_mat <- d2[,colnames(ceres) %in% data_mats[[i]],drop=F]
    }
    
    tmp_mat <- gather(as.data.frame(filtered_mat),key="Entrez ID","Gene Score") %>% mutate(.,Group=indiv_list)
    if (!is.na(sample_size)){
      tmp_mat <- tmp_mat[sample(nrow(tmp_mat),sample_size),]
    }
    data_mats[[indiv_list]] <- tmp_mat
    
  }
  
  long_data <- bind_rows(data_mats)
  long_data$tech <- "CRISPR"
  long_data$tech[grepl("RNAi_",long_data$Group)] <- "RNAi"
  long_data$profile <- "Weakly Selective"
  long_data$profile[grepl("_HVD",long_data$Group)] <- "High-variance"
  long_data$profile[grepl("_PD",long_data$Group)] <- "Pan-dependency"
  long_data$profile[grepl("_SSD",long_data$Group)] <- "Strongly Selective"
  long_data$profile[grepl("_ND",long_data$Group)] <- "Non-dependency"
  long_data$profile <- factor(long_data$profile,levels=c("Pan-dependency","High-variance","Strongly Selective","Weakly Selective","Non-dependency"))
  return(long_data)
}

#Comparison to use for density plots
#Full distributions
long_data <- get_data_long(gene_lists,ceres,d2,sample_size=NA)

mypal = c("Non-dependency"="#00204DFF","Strongly Selective"="#145A32","Pan-dependency"="#C71000B2","Weakly Selective"="Gray","High-variance"="#FF6F00B2")

####################################################

ridge_rug <- ggplot(data=long_data,aes(x=`Gene Score`,y=profile,fill=profile,color=profile))+
  geom_density_ridges(alpha=.7,scale = .95, rel_min_height = .01) +
  scale_shape_identity() +
  # geom_point(data=eg_data,mapping=aes(shape ="|",color=profile),size=3) +
  theme_ridges(center = TRUE) +
  scale_y_discrete(expand = c(.01, 0)) +
  scale_x_continuous(labels=c("-2","-1","0",""),limits=c(-2,1)) +
  scale_fill_manual(values=mypal,name = "Densities")+
  scale_color_manual(values=mypal,name = "Points", labels = c("Pan-Dependency"="SF3B1","SSD"="PAX8","No Dependency"="A1BG","Selective"="NOTCH3","High-variance"="ADAR"))+
  theme(legend.box = 'horizontal') +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        plot.title=element_text(size=16,hjust = 0.5)) +
  labs(x="Gene Score",y="",caption="") +
  facet_grid(. ~ tech, scales = "free", space = "free") +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(margin = margin(.25,0,.25,0, "cm"),face = "bold")) +
  theme(plot.margin = margin(l=0,r=0,unit="cm"))

ridge_rug + theme(legend.position="none")
ggsave(file.path("figures","gene_deps_profiles_density_joyplot.pdf"),height=4.1,width=4.6)


####################################################

#Absolute bars

genedeps_counts <- data.frame(Group=c("Non-dependency","Pan-dependency","Weakly Selective","Strongly Selective","High-variance"),stringsAsFactors = F)
genedeps_counts$`CRISPR` <- 0
genedeps_counts$`RNAi` <- 0
rownames(genedeps_counts) <- genedeps_counts$Group
genedeps_counts %<>% dplyr::select(.,-Group)
genedeps_counts <- t(genedeps_counts)
genedeps_counts <- data.frame(genedeps_counts,check.names = F)

genesets <- gene_lists
genedeps_counts["CRISPR","Non-dependency"] <- length(genesets[["CRISPR_ND"]])
genedeps_counts["CRISPR","Pan-dependency"] <- length(genesets[["CRISPR_PD"]])
genedeps_counts["CRISPR","Weakly Selective"] <- length(genesets[["CRISPR_WSD"]])
genedeps_counts["CRISPR","Strongly Selective"] <- length(genesets[["CRISPR_SSD"]])
genedeps_counts["CRISPR","High-variance"] <- length(genesets[["CRISPR_HVD"]])

genedeps_counts["RNAi","Non-dependency"] <- length(genesets[["RNAi_ND"]])
genedeps_counts["RNAi","Pan-dependency"] <- length(genesets[["RNAi_PD"]])
genedeps_counts["RNAi","Weakly Selective"] <- length(genesets[["RNAi_WSD"]])
genedeps_counts["RNAi","Strongly Selective"] <- length(genesets[["RNAi_SSD"]])
genedeps_counts["RNAi","High-variance"] <- length(genesets[["RNAi_HVD"]])

genedeps_counts <- t(genedeps_counts)
genedeps_counts <- data.frame(genedeps_counts,check.names = F)
genedeps_labels <- genedeps_counts
genedeps_labels$ID <- rownames(genedeps_labels)
genedeps_labels <- melt(genedeps_labels)

bars_pal <- c("CRISPR"="#0099B4FF","Overlap"="#868686FF","RNAi"="#925E9FFF")

genedeps_labels$ID <- factor(genedeps_labels$ID,levels=c("Pan-dependency","High-variance","Strongly Selective","Weakly Selective","Non-dependency"))

#No log 
double_nodep <- subset(genedeps_labels,ID %in% c("Weakly Selective","Non-dependency"))
double_nodep <- rbind(double_nodep,subset(genedeps_labels,ID %in% c("Non-dependency")))
double_nodep$ID <- as.character(double_nodep$ID)
double_nodep$ID[5:6] <- c("Placeholder","Placeholder")
double_nodep$value[5:6] <- c(0,0)
double_nodep$ID <- factor(double_nodep$ID,levels=c("Placeholder","Weakly Selective","Non-dependency"))
high_count <- ggplot(double_nodep, aes(x=ID, y=value, fill=variable)) +
  geom_bar(stat="identity",width = 0.5, position = position_dodge(),color="black",size=1) +
  coord_flip() +
  theme_minimal(base_size=12) +
  scale_fill_manual(values=bars_pal,name='',labels=c("CRISPR","RNAi")) +
  theme(legend.position = "none",legend.spacing.x = unit(.1, 'cm')) +
  # theme(axis.text.y=element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Number of Genes")

low_count <- ggplot(subset(genedeps_labels,ID %in% c("Pan-dependency","High-variance","Strongly Selective")), aes(x=ID, y=value, fill=variable)) +
  geom_bar(stat="identity",width = 0.5, position = position_dodge(),color="black",size=1) +
  coord_flip() +
  theme_minimal(base_size=12) +
  scale_fill_manual(values=bars_pal,name='',labels=c("CRISPR","RNAi")) +
  theme(legend.position = "none",legend.spacing.x = unit(.1, 'cm')) +
  # theme(axis.text.y=element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Number of Genes")

grobs <- rbind(ggplotGrob(high_count), ggplotGrob(low_count), size = "first")
grid.newpage()
pdf(file.path("figures","gene_deps_profiles_abs_bars.pdf"), height=3.9,width=4)
grid.draw(grobs)
dev.off()


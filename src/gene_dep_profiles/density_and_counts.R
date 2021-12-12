

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

#Add high variance to table of dependency metrics
#Table needs numDeps(50-100), Common_Essential, HVD, LRT
hv_res <- fread("/Users/mburger/dynamic-duo-biorxiv/figures/strongly_selective/highVar_classification.csv")
hv_res$entrez_id %<>% as.character(.)

master <- load.from.taiga(data.name='master-table-na-matched-6315', data.version=1, data.file='master-table')
master$entrez_id %<>% as.character(.)

master$CRISPR_HVD <- master$entrez_id %in% subset(hv_res,CRISPR_highVar)$entrez_id
master$RNAi_HVD <- master$entrez_id %in% subset(hv_res,RNAi_highVar)$entrez_id

#Replace CDS_ID with entrez ID
colnames(ceres) <- extract_entrez(colnames(ceres))
colnames(d2) <- extract_entrez(colnames(d2))

master$entrez_id %<>% as.character(.)
  
gene_lists <- list()
#ND
gene_lists[["crispr_nd"]] <- subset(master,`CRISPR_numDeps(50-100)` == 0)$entrez_id
gene_lists[["rnai_nd"]] <- subset(master,`RNAi_numDeps(50-100)` == 0)$entrez_id
  
#PD
gene_lists[["crispr_pd"]] <- subset(master,CRISPR_Common_Essential)$entrez_id
gene_lists[["rnai_pd"]] <- subset(master,RNAi_Common_Essential)$entrez_id
  
#HVD
gene_lists[["crispr_hvd"]] <- subset(master,CRISPR_HVD)$entrez_id
gene_lists[["rnai_hvd"]] <- subset(master,RNAi_HVD)$entrez_id
  
#SSD
gene_lists[["crispr_ssd"]] <- subset(master_table,CRISPR_LRT > 100)$entrez_id
gene_lists[["rnai_ssd"]] <- subset(master_table,RNAi_LRT > 100)$entrez_id
  
#SD
all_crispr <- unique(c(gene_lists[["crispr_nd"]],gene_lists[["crispr_pd"]],gene_lists[["crispr_ssd"]]))
all_rnai <- unique(c(gene_lists[["rnai_nd"]],gene_lists[["rnai_pd"]],gene_lists[["rnai_ssd"]]))
gene_lists[["crispr_sd"]] <- subset(master_table,!(entrez_id %in% all_crispr))$entrez_id
gene_lists[["rnai_sd"]] <- subset(master_table,!(entrez_id %in% all_rnai))$entrez_id

#### Example genes
# eg_lists <- list()
# #ND
# eg_lists[["crispr_nd"]] <- "1" #A1BG
# eg_lists[["rnai_nd"]] <- "1" 
#   
# #PD
# eg_lists[["crispr_pd"]] <- "23451" #SF3B1
# eg_lists[["rnai_pd"]] <- "23451" 
#   
# #SSD
# eg_lists[["crispr_ssd"]] <- "7849" #PAX8
# eg_lists[["rnai_ssd"]] <- "7849"
#   
# #SD 
# eg_lists[["crispr_sd"]] <- "4854" #NOTCH3
# eg_lists[["rnai_sd"]] <- "4854"
#   
# #HVD
# eg_lists[["crispr_hvd"]] <- "103" #ADAR
# eg_lists[["rnai_hvd"]] <- "103"
  

get_data_long <- function(data_mats,ceres,d2,sample_size){
  for (i in seq_along(data_mats)){
    indiv_list <- names(data_mats)[i]
    if (grepl("crispr",indiv_list)){
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
  long_data$tech[grepl("rnai_",long_data$Group)] <- "RNAi"
  long_data$profile <- "Weakly Selective"
  long_data$profile[grepl("_hvd",long_data$Group)] <- "High-variance"
  long_data$profile[grepl("_pd",long_data$Group)] <- "Pan-dependency"
  long_data$profile[grepl("_ssd",long_data$Group)] <- "Strongly Selective"
  long_data$profile[grepl("_nd",long_data$Group)] <- "Non-dependency"
  long_data$profile <- factor(long_data$profile,levels=c("Pan-dependency","High-variance","Strongly Selective","Weakly Selective","Non-dependency"))
  return(long_data)
}

#Comparison to use for density plots
#Full distributions
long_data <- get_data_long(gene_lists,ceres,d2,sample_size=NA)

#Example points
# eg_data <- get_data_long(eg_lists,ceres,d2,sample_size=NA)

mypal = c("Non-dependency"="#00204DFF","Strongly Selective"="#145A32","Pan-dependency"="#C71000B2","Weakly Selective"="Gray","High-variance"="#FF6F00B2")


####################################################
#Exploring facet version

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

genedeps_counts <- data.frame(Group=c("No Dependency","Pan-Dependency","Selective","SSD","High-variance"),stringsAsFactors = F)
genedeps_counts$`CRISPR` <- 0
genedeps_counts$`RNAi` <- 0
rownames(genedeps_counts) <- genedeps_counts$Group
genedeps_counts %<>% dplyr::select(.,-Group)
genedeps_counts <- t(genedeps_counts)
genedeps_counts <- data.frame(genedeps_counts,check.names = F)

genesets <- gene_lists[[comp]]
genedeps_counts["CRISPR","No Dependency"] <- length(genesets[["crispr_nd"]])
genedeps_counts["CRISPR","Pan-Dependency"] <- length(genesets[["crispr_pd"]])
genedeps_counts["CRISPR","Selective"] <- length(genesets[["crispr_sd"]])
genedeps_counts["CRISPR","SSD"] <- length(genesets[["crispr_ssd"]])
genedeps_counts["CRISPR","High-variance"] <- length(genesets[["crispr_hvd"]])

genedeps_counts["RNAi","No Dependency"] <- length(genesets[["rnai_nd"]])
genedeps_counts["RNAi","Pan-Dependency"] <- length(genesets[["rnai_pd"]])
genedeps_counts["RNAi","Selective"] <- length(genesets[["rnai_sd"]])
genedeps_counts["RNAi","SSD"] <- length(genesets[["rnai_ssd"]])
genedeps_counts["RNAi","High-variance"] <- length(genesets[["rnai_hvd"]])

genedeps_counts <- t(genedeps_counts)
genedeps_counts <- data.frame(genedeps_counts,check.names = F)
genedeps_labels <- genedeps_counts
genedeps_labels$ID <- rownames(genedeps_labels)
genedeps_labels <- melt(genedeps_labels)

bars_pal <- c("CRISPR"="#0099B4FF","Overlap"="#868686FF","RNAi"="#925E9FFF")

genedeps_labels$ID <- factor(genedeps_labels$ID,levels=c("Pan-Dependency","High-variance","SSD","Selective","No Dependency"))
genedeps_labels$trans_value <- log2(genedeps_labels$value)

l <- c(100,400,1600,6400)
b <- log2(l)

bar_p <- ggplot(genedeps_labels, aes(x=ID, y=trans_value, fill=variable)) +
  geom_bar(stat="identity",width = 0.7, position = position_dodge(),color="black",size=1) +
  coord_flip() +
  theme_minimal(base_size=12) +
  scale_fill_manual(values=bars_pal,name='',labels=c("CRISPR","RNAi")) +
  scale_y_continuous(breaks=b,labels=l,limits = c(log2(50),max(genedeps_labels$trans_value)),oob = rescale_none) + #limits = c(b[1],b[length(b)])
  # scale_y_continuous(breaks=b,labels=l,oob = rescale_none) +
  theme(legend.position = "none",legend.spacing.x = unit(.1, 'cm')) +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=11,angle = 45),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Number of Genes")

ggsave(plot=bar_p,"/Users/mburger/dynamic-duo-biorxiv/figures/gene_deps_summary/dependency_profile_abs_bars.pdf",height=3.9,width=2.5)

#No log 
double_nodep <- subset(genedeps_labels,ID %in% c("Selective","No Dependency"))
double_nodep <- rbind(double_nodep,subset(genedeps_labels,ID %in% c("No Dependency")))
double_nodep$ID <- as.character(double_nodep$ID)
double_nodep$ID[5:6] <- c("Placeholder","Placeholder")
double_nodep$value[5:6] <- c(0,0)
double_nodep$ID <- factor(double_nodep$ID,levels=c("Placeholder","Selective","No Dependency"))
high_count <- ggplot(double_nodep, aes(x=ID, y=value, fill=variable)) +
  geom_bar(stat="identity",width = 0.5, position = position_dodge(),color="black",size=1) +
  coord_flip() +
  theme_minimal(base_size=12) +
  scale_fill_manual(values=bars_pal,name='',labels=c("CRISPR","RNAi")) +
  # scale_y_continuous(breaks=b,labels=l,limits = c(log2(50),max(genedeps_labels$trans_value)),oob = rescale_none) + #limits = c(b[1],b[length(b)])
  # scale_y_continuous(breaks=b,labels=l,oob = rescale_none) +
  theme(legend.position = "none",legend.spacing.x = unit(.1, 'cm')) +
  theme(axis.text.y=element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Number of Genes")

low_count <- ggplot(subset(genedeps_labels,ID %in% c("Pan-Dependency","High-variance","SSD")), aes(x=ID, y=value, fill=variable)) +
  geom_bar(stat="identity",width = 0.5, position = position_dodge(),color="black",size=1) +
  coord_flip() +
  theme_minimal(base_size=12) +
  scale_fill_manual(values=bars_pal,name='',labels=c("CRISPR","RNAi")) +
  # scale_y_continuous(breaks=b,labels=l,limits = c(log2(50),max(genedeps_labels$trans_value)),oob = rescale_none) + #limits = c(b[1],b[length(b)])
  # scale_y_continuous(breaks=b,labels=l,oob = rescale_none) +
  theme(legend.position = "none",legend.spacing.x = unit(.1, 'cm')) +
  theme(axis.text.y=element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Number of Genes")

grobs <- rbind(ggplotGrob(high_count), ggplotGrob(low_count), size = "first")
grid.newpage()
pdf("/Users/mburger/dynamic-duo-biorxiv/figures/gene_deps_summary/dependency_profile_abs_bars_nolog_v2.pdf", height=3.9,width=2.5)
grid.draw(grobs)
dev.off()


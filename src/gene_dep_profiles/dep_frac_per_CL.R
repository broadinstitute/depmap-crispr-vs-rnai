
# require(cowplot)

source("src/packages_paths.R")

master_table <- fread(file.path("tables","Supplemental-Table-2.csv"))
master_table$entrez_id %<>% as.character(.)
master_table$CDS_ID <- paste0(master_table$symbol," (",master_table$entrez_id,")")

master_table$crispr_dep_frac <- master_table$`CRISPR_numDeps(50-100)` / (master_table$`CRISPR_numDeps(0-50)` + master_table$`CRISPR_numDeps(50-100)`)
master_table$rnai_dep_frac <- master_table$`RNAi_numDeps(50-100)` / (master_table$`RNAi_numDeps(0-50)` + master_table$`RNAi_numDeps(50-100)`)

gene_sets <- list()

gene_sets[["crispr_nd"]] <- subset(master_table,CRISPR_class == "Non-dependency")$CDS_ID
gene_sets[["rnai_nd"]] <- subset(master_table,RNAi_class == "Non-dependency")$CDS_ID

gene_sets[["crispr_pd"]] <- subset(master_table,CRISPR_class == "Pan-dependency")$CDS_ID
gene_sets[["rnai_pd"]] <- subset(master_table,RNAi_class == "Pan-dependency")$CDS_ID

gene_sets[["crispr_ssd"]] <- subset(master_table,CRISPR_class == "Strongly Selective Dependency")$CDS_ID
gene_sets[["rnai_ssd"]] <- subset(master_table,RNAi_class == "Strongly Selective Dependency")$CDS_ID

gene_sets[["crispr_wsd"]] <- subset(master_table,CRISPR_class == "Weakly Selective Dependency")$CDS_ID
gene_sets[["rnai_wsd"]] <- subset(master_table,RNAi_class == "Weakly Selective Dependency")$CDS_ID

gene_sets[["crispr_hvd"]] <- subset(master_table,CRISPR_class == "High-variance Dependency")$CDS_ID
gene_sets[["rnai_hvd"]] <- subset(master_table,RNAi_class == "High-variance Dependency")$CDS_ID

#import dependency probabilities 
crispr_pr <- fread(file.path(data_processed,"dependency-probability-crispr-matched.csv")) %>% column_to_rownames(.,var="Row.name")

rnai_pr <- fread(file.path(data_processed,"dependency-probability-rnai-matched.csv")) %>% column_to_rownames(.,var="Row.name")
all(rownames(crispr_pr) == rownames(rnai_pr))

pr_datasets <- list("crispr"=crispr_pr,
                    "rnai"=rnai_pr)

#### Average dependency per cell line stats
mean(rowSums(crispr_pr > .5,na.rm=T) / rowSums(!is.na(crispr_pr)))
mean(rowSums(rnai_pr > .5,na.rm=T) / rowSums(!is.na(rnai_pr)))

total_screened <- rowSums(!is.na(crispr_pr))
results_df <- data.frame(CL=names(total_screened),screened=total_screened)

all_counts <- list()
for (gs in names(gene_sets)){
  tech <- strsplit(gs,"_")[[1]][1]
  tmp_pr <- pr_datasets[[tech]]
  tmp_pr <- tmp_pr[,colnames(tmp_pr) %in% gene_sets[[gs]]]
  dep_count <- rowSums(tmp_pr > .5,na.rm=T)
  res_df <- as.data.frame(dep_count)
  colnames(res_df) <- gs
  all_counts[[gs]] <- res_df
}

all_counts <- bind_cols(all_counts)
crispr_counts <- all_counts[,grepl("crispr",colnames(all_counts))]
rnai_counts <- all_counts[,grepl("rnai",colnames(all_counts))]

crispr_frac <- crispr_counts / total_screened
crispr_frac %<>% dplyr::select(.,crispr_pd,crispr_ssd,crispr_wsd,crispr_hvd)
colnames(crispr_frac) <- gsub("crispr_","",colnames(crispr_frac))
crispr_frac %<>% rownames_to_column(.,var="CL") %>% as.data.table(.)
crispr_frac <- melt.data.table(crispr_frac,id.vars="CL")
crispr_frac$tech <- "CRISPR"

rnai_frac <- rnai_counts / total_screened
rnai_frac %<>% dplyr::select(.,rnai_pd,rnai_ssd,rnai_wsd,rnai_hvd)
colnames(rnai_frac) <- gsub("rnai_","",colnames(rnai_frac))
rnai_frac %<>% rownames_to_column(.,var="CL") %>% as.data.table(.)
rnai_frac <- melt.data.table(rnai_frac,id.vars="CL")
rnai_frac$tech <- "RNAi"

stats <- bind_rows(crispr_frac,rnai_frac)

total_plots <- list()
for (t in c("CRISPR","RNAi")){
  test_stats <- subset(stats,tech == t)
  
  mypal = c("nd"="white","hvd"="#F37021","ssd"="#175B33","pd"="#C62026","wsd"="Gray")

  test_stats$variable <- factor(test_stats$variable,levels=rev(c("nd","ssd","wsd","hvd","pd")))

  #get dep frac order
  tmp_pr <- pr_datasets[[tolower(t)]]
  total_deps <- rowSums(tmp_pr > .5,na.rm=T)
  total_dep_frac <- total_deps/total_screened 
  total_dep_frac <- sort(total_dep_frac,decreasing=T)
  test_stats$CL <- factor(test_stats$CL,levels=names(total_dep_frac))
  
  u_total_pd <- mean(subset(test_stats,variable == "pd")$value)
  
  total_plots[[t]] <- ggplot(test_stats,aes(x=CL,y=value)) +
    geom_bar(aes(fill = variable), width = 1, position = position_stack(reverse = TRUE),stat="identity") +
    geom_hline(yintercept=u_total_pd,color=mypal["pd"]) +
    scale_fill_manual(values=mypal) +
    theme_classic(base_size=11) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("")
  
}

total_row <- cowplot::plot_grid(total_plots[["CRISPR"]] + xlab("CRISPR") + ylab("Dep. fraction per CL") + scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits=c(0,.3)) + expand_limits(x = -20, y = 0),
                       total_plots[["RNAi"]] + xlab("RNAi") + ylab("Dep. fraction per CL") + scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits=c(0,.3)) + expand_limits(x = -20, y = 0),
                       align = 'h',
                       nrow=1,
                       rel_widths=c(1,1),
                       axis="tb")

ggsave(plot=total_row,file.path("figures","gene_deps_profiles_totalfrac_perCL_stackedbar.pdf"),height=2,width=3)

u_dep_crispr_pd <- mean(subset(stats, (tech == "CRISPR") & (variable == "pd"))$value)
sd_dep_crispr_pd <- sd(subset(stats, (tech == "CRISPR") & (variable == "pd"))$value)

u_dep_rnai_pd <- mean(subset(stats, (tech == "RNAi") & (variable == "pd"))$value)
sd_dep_rnai_pd <- sd(subset(stats, (tech == "RNAi") & (variable == "pd"))$value)


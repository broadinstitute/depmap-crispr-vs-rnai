library(stringr) 
library(ggrepel)

source("src/packages_paths.R")

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
t2$CDS_ID <- entrez_to_cds(t2$entrez_id,hgnc)

t2$pandep_group <- "init"
t2$pandep_group[t2$CRISPR_PD & (!t2$RNAi_PD)] <- "CRISPR-specific"
t2$pandep_group[t2$CRISPR_PD & t2$RNAi_PD] <- "shared"
t2 %<>% dplyr::select(.,CDS_ID,entrez_id,pandep_group,symbol)

msigdb_sets <- readRDS(file.path(data_raw,"gene-set-msigdb-onehot.rds"))
msigdb_sets <- msigdb_sets[rownames(msigdb_sets) %in% t2$entrez_id,]

shared_pd <- subset(t2,pandep_group == "shared")$entrez_id
notshared_pd <- subset(t2,pandep_group != "shared")$entrez_id
distinct_pd <- subset(t2,pandep_group == "CRISPR-specific")$entrez_id
notdistinct_pd <- subset(t2,pandep_group != "CRISPR-specific")$entrez_id

df_multilib <- data.frame(gene_set=colnames(msigdb_sets),
                          shared_oddsratio=NA,
                          shared_pval=NA,
                          distinct_oddsratio=NA,
                          distinct_pval=NA,
                          stringsAsFactors = F)

#For each set, in group is shared or distinct and out group is the remainder of the high confidence dependency set
for (i in 1:ncol(msigdb_sets)){
  
  gset = colnames(msigdb_sets)[i]
  
  in_set = rownames(msigdb_sets)[msigdb_sets[,gset]]
  out_set = rownames(msigdb_sets)[!msigdb_sets[,gset]]
  
  in_shared = length(intersect(in_set,shared_pd))
  in_notshared = length(intersect(in_set,notshared_pd))
  
  out_shared = length(intersect(out_set,shared_pd))
  out_notshared = length(intersect(out_set,notshared_pd))
  
  #IN/OUT are the columns
  #Shared/!shared are rows
  m = matrix(c(in_shared,out_shared,in_notshared,out_notshared),byrow=T,nrow=2,ncol=2) + 1
  colnames(m) <- c("IN","OUT")
  rownames(m) <- c("Shared","Not Shared")
  f_test <- fisher.test(m,alternative="greater")
  df_multilib[i,"shared_pval"] = f_test$p.value
  df_multilib[i,"shared_oddsratio"] = f_test$estimate
  
  in_distinct = length(intersect(in_set,distinct_pd))
  in_notdistinct = length(intersect(in_set,notdistinct_pd))
  
  out_distinct = length(intersect(out_set,distinct_pd))
  out_notdistinct = length(intersect(out_set,notdistinct_pd))
  
  m = matrix(c(in_distinct,out_distinct,in_notdistinct,out_notdistinct),byrow=T,nrow=2,ncol=2) + 1
  colnames(m) <- c("IN","OUT")
  rownames(m) <- c("Distinct","Not Distinct")
  f_test = fisher.test(m,alternative="greater")
  df_multilib[i,"distinct_pval"] = f_test$p.value
  df_multilib[i,"distinct_oddsratio"] = f_test$estimate
}

df_multilib$OR_diff <- df_multilib$shared_oddsratio - df_multilib$distinct_oddsratio
df_multilib$shared_qval <- p.adjust(df_multilib$shared_pval,method="BH")
df_multilib$distinct_qval <- p.adjust(df_multilib$distinct_pval,method="BH")

write_csv(df_multilib,file.path(data_processed,"pandependency_agreement_enrichment_of_shared_or_distinct.csv"))

#Individual Set Collection Plots
c_names <- c("BIOCARTA","GO","HALLMARK","KEGG","REACTOME")
c_pal <- c("shared"="#A44E78","distinct"="darkblue","none"="grey70")

for (i in 1:length(c_names)){
  tmp_c <- c_names[i] 
  
  plot_df <- subset(df_multilib,grepl(paste0("^",tmp_c),df_multilib$gene_set))
  plot_df$hit <- "none"
  plot_df$hit[plot_df$shared_qval < .05] <- "shared"
  plot_df$hit[plot_df$distinct_qval < .05] <- "distinct"
  
  plot_df$clean_names <- gsub(paste0(tmp_c,"_"),"",plot_df$gene_set)
  plot_df$clean_names <- gsub("_"," ",plot_df$clean_names)
  plot_df$clean_names <- tolower(plot_df$clean_names)
  plot_df$clean_names <- str_to_title(plot_df$clean_names)
  
  p <- ggplot(plot_df,aes(x=shared_oddsratio,y=distinct_oddsratio,color=hit)) +
    geom_point(data=subset(plot_df,(hit=="none")),size=3) +
    geom_point(data=subset(plot_df,hit=="shared"),size=3) +
    geom_point(data=subset(plot_df,hit=="distinct"),size=3) +
    scale_color_manual(values=c_pal) +
    xlab("Shared Pan-Dep. Odds Ratio") +
    ylab("Distinct CRISPR Pan-Dep. Odds Ratio") +
    ggtitle(tmp_c) +
    theme_bw(base_size=9) +
    theme(legend.position = "none")
  
  if (tmp_c %in% c("KEGG","BIOCARTA")){
    p <- p + geom_text_repel(data=subset(plot_df,(hit != "none") & ((distinct_oddsratio > 20) | (shared_oddsratio > 15))),aes(label=clean_names)) 
  } else if (tmp_c %in% c("REACTOME","GO")) {
    p <- p + geom_text_repel(data=subset(plot_df,(hit != "none") & ((distinct_oddsratio > 30) | (shared_oddsratio > 60))),aes(label=clean_names))
  }
  ggsave(plot=p,file.path("figures",paste0("pandependency_agreement_",tmp_c,"_odds_ratio_scatter.pdf")),height=2.5,width=3.1)
}



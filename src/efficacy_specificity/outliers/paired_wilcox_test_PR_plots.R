library(tidyverse)
library(magrittr)
library(data.table)
library(ggrepel)

t3 <- fread("tables/Supplemental-Table-1.csv")
t3 %<>% subset(.,CRISPR_agreement & RNAi_agreement)
t3 %<>% subset(.,!(CRISPR_shared_nondep & RNAi_shared_nondep))
conf_genes <- paste0(t3$symbol," (",t3$entrez_id,")")

crispr.gs <- fread("data/raw/gene-effect-scaled-crispr-matched.csv") %>% 
  column_to_rownames(.,var="Row.name")
rnai.gs <- fread("data/raw/gene-effect-scaled-rnai-matched.csv") %>%
  column_to_rownames(.,var="Row.name")

crispr.gs <- crispr.gs[,colnames(crispr.gs) %in% conf_genes]
rnai.gs <- rnai.gs[,colnames(rnai.gs) %in% conf_genes]

CRISPR.PR <- fread("data/processed/dependency-probability-crispr-matched.csv") %>% 
  column_to_rownames(.,var="Row.name")
RNAi.PR <- fread("data/processed/dependency-probability-rnai-matched.csv") %>%
  column_to_rownames(.,var="Row.name")


CRISPR.PR <- CRISPR.PR[,colnames(CRISPR.PR) %in% conf_genes]
RNAi.PR <- RNAi.PR[,colnames(RNAi.PR) %in% conf_genes]

conf_genes <- intersect(conf_genes,colnames(CRISPR.PR))


stats <- list()
for (g in conf_genes){
  
  x <- RNAi.PR[,g]
  y <- CRISPR.PR[,g]
  
  rest <- t.test(x,y,paired=T,var.equal=F)
  resw <- wilcox.test(x,y,paired=T)
  
  stats[[g]] <- data.frame(gene=g,
                           statistic=as.numeric(resw$statistic),
                           p.value=as.numeric(resw$p.value),
                           estimate=as.numeric(rest$estimate),
                           stringsAsFactors = F)
}
stats <- bind_rows(stats)

stats$Y <- -log10(stats$p.value)

mypal <- c("#A44E78","#4186C4","#E49042")

point_label_font = 8
quadrant_label_font = 11

stats$labels <- gsub(" .*","",stats$gene)
stats$labels[(stats$estimate > -.9) & (stats$estimate < .1)] <- NA
stats$labels[(stats$Y < 10)] <- NA

# write_csv(stats,file="src/efficacy_specificity/outliers/mean_diff_stats.csv")

#### Volcano figure

exponents <- c(20,40,60)
b <- 10^-(exponents)
l <- paste0("10^{-",exponents,"}")
b <- -log10(b)

ggplot(stats,aes(estimate,Y)) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill=mypal[2],alpha=0.2) +
  annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill=mypal[1],alpha=0.2) +
  geom_vline(xintercept = 0,linetype="dotted",color="gray") +
  geom_point(alpha=.7,color="darkblue",size=1.5) +
  geom_text_repel(data = subset(stats,!is.na(labels)),
                  size=(point_label_font / ggplot2:::.pt),
                  aes(estimate,Y,label = labels),
                  segment.colour="grey50",
                  segment.size = 0.25,
                  force = 20,
                  min.segment.length = 0) +
  theme_bw(base_size=11) +
  xlim(c(min(-1,min(stats$estimate)),max(1,max(stats$estimate)))) +
  xlab("Mean Difference in Dependency Probability\n(RNAi - CRISPR)") +
  ylab("P Value") +
  scale_y_continuous(breaks=b,labels=parse(text=l)) +
  theme(legend.position="none") +
  annotate("text", x = -1, y = 5, label = "Stronger CRISPR\nDependency",size=(quadrant_label_font / ggplot2:::.pt),hjust=0) + 
  annotate("text", x = 1, y = 5, label = "Stronger RNAi\nDependency",size=(quadrant_label_font / ggplot2:::.pt),hjust=1) 

ggsave("figures/efficacy_specificity_PR_volcano.pdf",width=6.8,height=5.5)
 
##### mean CRISPR -- mean RNAi scatter

crispr_gs <- fread("data/raw/gene-effect-scaled-crispr-matched.csv") %>% column_to_rownames(.,var="Row.name")
rnai_gs <- fread("data/raw/gene-effect-scaled-rnai-matched.csv") %>% column_to_rownames(.,var="Row.name")

crispr_pd <- fread("data/processed/pandependency-score-crispr-matched.csv") %>% subset(.,Common_Essential)
rnai_pd <- fread("data/processed/pandependency-score-rnai-matched.csv") %>% subset(.,Common_Essential)

cross_df <- data.frame(gene=colnames(crispr_gs),
                       crispr_mean=colMeans(crispr_gs,na.rm=T),
                       rnai_mean=colMeans(rnai_gs,na.rm=T))

sig_stats <- stats

sig_stats %<>% left_join(.,cross_df,by="gene")

sig_stats$pd_status <- "Neither"
sig_stats$pd_status[sig_stats$gene %in% crispr_pd$Row.name] <- "CRISPR"
sig_stats$pd_status[sig_stats$gene %in% rnai_pd$Row.name] <- "RNAi"
sig_stats$pd_status[sig_stats$gene %in% intersect(rnai_pd$Row.name,crispr_pd$Row.name)] <- "Both"

ggplot(sig_stats,aes(x=crispr_mean,y=rnai_mean,color=pd_status)) +
  geom_point() +
  geom_abline(slope=1,linetype="dashed") +
  xlab("CRISPR mean gene effect") +
  ylab("RNAi mean gene effect") +
  geom_text_repel(data=subset(sig_stats,!is.na(sig_stats$labels)),
                  aes(x=crispr_mean,y=rnai_mean,label=labels),
                  min.segment.length = 0) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom")
ggsave("figures/efficacy_specificity_CRISPR_RNAi_uu_scatter.pdf",
       height=6,width=6)


##### mean KY -- mean RNAi scatter

crispr_gs <- fread("data/raw/gene-effect-scaled-crispr-ky.csv") %>% column_to_rownames(.,var="V1")
rnai_gs <- fread("data/raw/gene-effect-scaled-rnai-matched.csv") %>% column_to_rownames(.,var="Row.name")

cls <- intersect(rownames(crispr_gs),rownames(rnai_gs))
genes <- intersect(colnames(crispr_gs),colnames(rnai_gs))

crispr_gs <- crispr_gs[,genes]
rnai_gs <- rnai_gs[,genes]

crispr_pd <- fread("data/processed/pandependency-score-crispr-ky.csv") %>% subset(.,Common_Essential)
rnai_pd <- fread("data/processed/pandependency-score-rnai-matched.csv") %>% subset(.,Common_Essential)

cross_df <- data.frame(gene=colnames(crispr_gs),
                       crispr_mean=colMeans(crispr_gs,na.rm=T),
                       rnai_mean=colMeans(rnai_gs,na.rm=T))

sig_stats <- stats

sig_stats %<>% left_join(.,cross_df,by="gene")

sig_stats$pd_status <- "Neither"
sig_stats$pd_status[sig_stats$gene %in% crispr_pd$Row.name] <- "CRISPR"
sig_stats$pd_status[sig_stats$gene %in% rnai_pd$Row.name] <- "RNAi"
sig_stats$pd_status[sig_stats$gene %in% intersect(rnai_pd$Row.name,crispr_pd$Row.name)] <- "Both"

sig_stats$labels[sig_stats$estimate < 0] <- NA

ggplot(sig_stats,aes(x=crispr_mean,y=rnai_mean,color=pd_status)) +
  geom_point() +
  geom_abline(slope=1,linetype="dashed") +
  xlab("CRISPR mean gene effect") +
  ylab("RNAi mean gene effect") +
  geom_label_repel(data=subset(sig_stats,!is.na(sig_stats$labels)),
                   aes(x=crispr_mean,y=rnai_mean,label=labels),
                   min.segment.length = 0) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom")
ggsave("figures/efficacy_specificity_KY_RNAi_uu_scatter.pdf",
       height=6,width=6)



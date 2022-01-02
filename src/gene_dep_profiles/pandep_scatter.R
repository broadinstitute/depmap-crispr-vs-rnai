
library(ggrepel)

source("src/packages_paths.R")

hgnc <- fread(file.path(data_raw,"hgnc-complete-set.csv"))
hgnc$entrez_id %<>% as.character(.)
hgnc %<>% subset(.,locus_group == "protein-coding gene")

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

t2 %<>% dplyr::select(.,entrez_id,CRISPR_pandependency_score,CRISPR_PD,RNAi_pandependency_score,RNAi_PD)

crispr_thresh_1 <- max(subset(t2,CRISPR_PD)$CRISPR_pandependency_score)
crispr_thresh_2 <- min(subset(t2,!CRISPR_PD)$CRISPR_pandependency_score)
crispr_thresh <- mean(c(crispr_thresh_1,crispr_thresh_2))

rnai_thresh_1 <- max(subset(t2,RNAi_PD)$RNAi_pandependency_score)
rnai_thresh_2 <- min(subset(t2,!RNAi_PD)$RNAi_pandependency_score)
rnai_thresh <- mean(c(rnai_thresh_1,rnai_thresh_2))

mypal <- c("#A44E78","#4186C4","#E49042")

t2$pandep_group <- "none"
t2$pandep_group[t2$CRISPR_PD & !t2$RNAi_PD] <- "CRISPR-specific"
t2$pandep_group[!t2$CRISPR_PD & t2$RNAi_PD] <- "RNAi-specific"
t2$pandep_group[t2$CRISPR_PD & t2$RNAi_PD] <- "shared"

t2 %<>% add_column(.,CDS_ID=entrez_to_cds(t2$entrez_id,hgnc),.before=1)

crispr_only_count <- nrow(subset(t2,pandep_group == "CRISPR-specific"))
rnai_only_count <- nrow(subset(t2,pandep_group == "RNAi-specific"))
both_count <- nrow(subset(t2,pandep_group == "shared"))

quadrant_label_font <- 11

t2$color_group <- "hit"
t2$color_group[!t2$CRISPR_PD & !t2$RNAi_PD] <- "nohit"

t2$symbol <- entrez_to_symbol(t2$entrez_id,hgnc)

color_pal <- c("hit"="darkblue","nohit"="grey50")

ggplot(t2,aes(y=CRISPR_pandependency_score,x=RNAi_pandependency_score,color=color_group)) +
  annotate("rect", xmin=0, xmax=rnai_thresh_1, ymin=0, ymax=crispr_thresh, size=1,color="grey50",fill=NA) +
  annotate("rect", xmin=rnai_thresh_1, xmax=Inf, ymin=0, ymax=crispr_thresh, fill=mypal[2],alpha=0.2) +
  annotate("rect", xmin=0, xmax=rnai_thresh_1, ymin=crispr_thresh, ymax=Inf, fill=mypal[1],alpha=0.2) +
  geom_point(alpha=.7,size=1.5) +
  geom_text_repel(data=subset(t2,!CRISPR_PD & RNAi_PD),aes(y=CRISPR_pandependency_score,x=RNAi_pandependency_score,label=symbol),
                  color="black",
                  xlim=c(.18,0)) +
  scale_color_manual(values=color_pal) +
  scale_x_reverse() +
  scale_y_reverse() +
  xlab("RNAi Pan-dependency Score") +
  ylab("CRISPR Pan-dependency Score") +
  theme_bw(base_size=11) +
  theme(legend.position = "none") +
  annotate("text", x = .7, y = -.05, label = paste0("CRISPR specific (",crispr_only_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=1) +
  annotate("text", x = .15, y = .8, label = paste0("RNAi\nspecific\n(",rnai_only_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=0) +
  annotate("text", x = .1, y = -.05, label = paste0("Shared (",both_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=.5)
ggsave(file.path("figures","gene_deps_profiles_pandependency_scatter.pdf"),width=3.5,height=2.5)


source("src/packages_paths.R")

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

point_pal <- c("grey50"="grey50","darkblue"="darkblue")
t2$color_group <- "grey50" 
t2$color_group[t2$CRISPR_SSD | t2$RNAi_SSD] <- "darkblue"

normLRT <- dplyr::select(t2,symbol,
                         CRISPR_LRT,RNAi_LRT,color_group)


max_LRT <- max(c(normLRT$CRISPR_LRT,normLRT$RNAi_LRT),na.rm=T)
min_LRT <- min(c(normLRT$CRISPR_LRT,normLRT$RNAi_LRT),na.rm=T)
normLRT$X <- normLRT$CRISPR_LRT
normLRT$X[normLRT$X < 0] <- 0
normLRT$Y <- normLRT$RNAi_LRT
normLRT$Y[normLRT$Y < 0] <- 0
normLRT$X <- log2(normLRT$X + 100)
normLRT$Y <- log2(normLRT$Y + 100)
min_LRT <- log2(min_LRT+100)
max_LRT <- log2(max_LRT+100)
mypal <- c("#A44E78","#4186C4","#E49042")

se_labels <- subset(normLRT,normLRT$CRISPR_LRT > 170 | normLRT$RNAi_LRT > 160)

crispr_only <- subset(normLRT,normLRT$CRISPR_LRT > 100 & normLRT$RNAi_LRT < 100)
crispr_only_count <- nrow(crispr_only)
rnai_only <- subset(normLRT,normLRT$CRISPR_LRT < 100 & normLRT$RNAi_LRT > 100)
rnai_only_count <- nrow(rnai_only)
both <- subset(normLRT,normLRT$CRISPR_LRT > 100 & normLRT$RNAi_LRT > 100)
both_count <- nrow(both)

point_label_font = 6
quadrant_label_font = 11

breaks <- c(0,100,200,400,600)
labels <- as.character(breaks)
breaks <- log2(breaks + 100)

threshold_val <- 100
threshold <- log2(threshold_val + 100)

ggplot(normLRT,aes(X,Y)) +
  annotate("rect", xmin=threshold, xmax=Inf, ymin=threshold, ymax=Inf, size=1,color="grey50",fill=NA) +
  annotate("rect", xmin=threshold, xmax=Inf, ymin=breaks[1], ymax=threshold, fill=mypal[2],alpha=0.2) +
  annotate("rect", xmin=breaks[1], xmax=threshold, ymin=threshold, ymax=Inf, fill=mypal[1],alpha=0.2) +
  geom_abline(slope=1,intercept = 0,linetype="dotted",color="gray") +
  geom_point(aes(color=color_group),alpha=.7,size=1.5) +
  geom_text_repel(data = se_labels, 
                  size=(point_label_font / ggplot2:::.pt),
                  aes(X,Y,label = symbol),
                  segment.colour="grey50",
                  segment.size = 0.25,
                  force = 20,
                  min.segment.length = 0) +
  theme_bw(base_size=11) +
  scale_x_continuous(breaks=breaks,labels=labels,limits=c(min_LRT,max_LRT)) +
  scale_y_continuous(breaks=breaks,labels=labels,limits=c(min_LRT,max_LRT)) +
  xlab("CRISPR LRT Score") +
  ylab("RNAi LRT Score") +
  scale_colour_manual(values=point_pal) +
  theme(legend.position="none") +
  annotate("text", x = max_LRT, y = log2(50 + 100), label = paste0("CRISPR\n(",crispr_only_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=1) +
  annotate("text", x = log2(25 + 100), y = log2(500+100), label = paste0("RNAi\n(",rnai_only_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=0) +
  annotate("text", x = log2(500+100), y = log2(500+100), label = paste0("Shared\n(",both_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=.5)
ggsave(file.path("figures","gene_deps_profiles_LRT_scatter.pdf"),width=3.5,height=2.5)

       
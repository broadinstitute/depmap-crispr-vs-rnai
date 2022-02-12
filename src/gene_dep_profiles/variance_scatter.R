
source("src/packages_paths.R")

library(ggrepel)

t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

crispr_var_thresh <- min(subset(t2,CRISPR_HVD)$CRISPR_dependency_variance)
rnai_var_thresh <- min(subset(t2,RNAi_HVD)$RNAi_dependency_variance)

var_df <- dplyr::select(t2,entrez_id,symbol,CRISPR_HVD,RNAi_HVD,CRISPR_dependency_variance,RNAi_dependency_variance)

var_df$color_group <- "grey50" 
var_df$color_group[var_df$CRISPR_HVD | var_df$RNAi_HVD] <- "darkblue"

max_var <- max(c(var_df$CRISPR_dependency_variance,var_df$RNAi_dependency_variance))

point_labels <- subset(var_df,(CRISPR_dependency_variance > .12) | (RNAi_dependency_variance > .12) | ((CRISPR_dependency_variance > .09) & (RNAi_dependency_variance > .09)))

point_label_font = 6
quadrant_label_font = 11

crispr_only <- subset(var_df,CRISPR_HVD & !RNAi_HVD)
crispr_only_count <- nrow(crispr_only)
rnai_only <- subset(var_df,!CRISPR_HVD & RNAi_HVD)
rnai_only_count <- nrow(rnai_only)
both <- subset(var_df,CRISPR_HVD & RNAi_HVD)
both_count <- nrow(both)

point_pal <- c("grey50"="grey50","darkblue"="darkblue")
mypal <- c("#A44E78","#4186C4","#E49042")
ggplot(var_df,aes(x=CRISPR_dependency_variance,y=RNAi_dependency_variance)) +
  annotate("rect", xmin=crispr_var_thresh, xmax=Inf, ymin=rnai_var_thresh, ymax=Inf, size=1,color="black",fill=NA) +
  annotate("rect", xmin=crispr_var_thresh, xmax=Inf, ymin=0, ymax=rnai_var_thresh, fill=mypal[2],alpha=0.2) +
  annotate("rect", xmin=0, xmax=crispr_var_thresh, ymin=rnai_var_thresh, ymax=Inf, fill=mypal[1],alpha=0.2) +
  geom_abline(slope=1,intercept = 0,linetype="dotted",color="gray") +
  geom_point(aes(color=color_group),alpha=.7,size=1.5) +
  geom_text_repel(data = point_labels,
                  size=(point_label_font / ggplot2:::.pt),
                  aes(CRISPR_dependency_variance,RNAi_dependency_variance,label = symbol),
                  segment.colour="grey50",
                  segment.size = 0.25,
                  force = 20,
                  min.segment.length = 0) +
  theme_bw(base_size=11) +
  xlim(c(0,max_var)) +
  ylim(c(0,max_var)) +
  xlab("CRISPR Variance") +
  ylab("RNAi Variance") +
  scale_colour_manual(values=point_pal) +
  theme(legend.position="none") +
  annotate("text", x = .15, y = .025, label = paste0("CRISPR\n(",crispr_only_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=.5) +
  annotate("text", x = .02, y = .15, label = paste0("RNAi\n(",rnai_only_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=0) +
  annotate("text", x = .15, y = .15, label = paste0("Shared\n(",both_count,")"),size=(quadrant_label_font / ggplot2:::.pt),hjust=.5)
ggsave(file.path("figures","gene_deps_profiles_variance_scatter.pdf"),width=3.5,height=2.5)

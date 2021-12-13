
source("src/packages_paths.R")

library(ggalluvial)

high_conf <- fread(file.path("tables","Supplemental-Table-1.csv"))
high_conf$entrez_id %<>% as.character(.)
high_conf %<>% subset(.,high_confidence)

dep_class <- fread(file.path("tables","Supplemental-Table-2.csv"))
dep_class$entrez_id %<>% as.character(.)

#Filter for high-confidence dependencies
dep_class %<>% subset(.,entrez_id %in% high_conf$entrez_id)

sum((dep_class$CRISPR_class == "Non-dependency") & (dep_class$RNAi_class == "Non-dependency"))
dep_class %<>% subset(.,!((CRISPR_class == "Non-dependency") & (RNAi_class == "Non-dependency")))

#####################  Remove genes that are non-dependency using both CRISPR and RNAi 
common_alluvial <- dep_class
common_alluvial %<>% dplyr::select(.,entrez_id,CRISPR_class,RNAi_class)
ce_union_n <- nrow(common_alluvial)
common_alluvial <- melt(common_alluvial,id="entrez_id")
common_alluvial$variable <- gsub("RNAi_class","RNAi",common_alluvial$variable)
common_alluvial$variable <- gsub("CRISPR_class","CRISPR",common_alluvial$variable)
common_alluvial$variable <- factor(common_alluvial$variable,levels=c("RNAi","CRISPR"))
common_alluvial$value <- factor(common_alluvial$value,levels=c("Non-dependency","Weakly Selective Dependency","Strongly Selective Dependency","High-variance Dependency","Pan-dependency"))

groups_pal = c("Non-dependency"="#00204DFF","Strongly Selective Dependency"="#145A32","Pan-dependency"="#C71000B2","Weakly Selective Dependency"="grey","High-variance Dependency"="#FF6F00B2")


ce_plot <- ggplot(common_alluvial,
                  aes(x = variable, stratum = value, alluvium = entrez_id,
                      fill = value, label = value)) +
  scale_fill_manual(values=groups_pal) +
  geom_flow(stat = "alluvium", lode.guidance = "leftright",
            alpha=1,width=1/12) +
  geom_stratum(width = 1/12) +
  theme_classic(base_size=11) +
  theme(legend.position = "none") +
  xlab("") +
  scale_x_discrete(expand=c(.1,.1)) 

ggsave(plot=ce_plot,file.path("figures","gene_deps_profiles_flow_diagram.pdf"),width=2.5,height=2.5)


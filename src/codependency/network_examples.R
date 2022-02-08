library(igraph)

related <- fread(file.path(data_raw,"gene-set-related-features.csv"))
gs_df <- fread(file.path(data_raw,"gene-set-kegg.csv"))
related <- bind_rows(related,gs_df)
related %<>% subset(.,target != partner) #remove self 

res<- fread(file.path(data_processed,"codependency_network_enrichment.csv"))

res %<>% dplyr::select(.,target,tech,ks_ABS_pvalue)
res %<>% spread(.,key=tech,value=ks_ABS_pvalue) %>% column_to_rownames(.,var="target")
res$RNAi_imp <- -log10(res$RNAi) - -log10(res$CRISPR)
res$RNAi_SNF_imp <- -log10(res$SNF) - -log10(res$RNAi)

res$CRISPR_imp <- -log10(res$CRISPR) - -log10(res$RNAi)
res$CRISPR_SNF_imp <- -log10(res$SNF) - -log10(res$CRISPR)

d1 <- "RPS21 (6227)"
d2 <- "MED16 (10025)"

cor_data <- list(CRISPR="codependency_CRISPR_pearson_baseline.rds",
                 RNAi="codependency_RNAi_pearson_baseline.rds",
                 SNF="codependency-CRISPR-RNAi-SNF.rds")
cor_data <- lapply(cor_data,function(x){readRDS(file.path(data_processed,x))})

zscore_similarity <- function(x){
  diag(x) = NA
  s = sd(x,na.rm=T)
  x = x / s
  x = x / max(x,na.rm=T)
  diag(x) = 2
  return(x)
}

cor_data <- lapply(cor_data,zscore_similarity)

network_dfs_1 <- list()
network_dfs_2 <- list()
for (dset in names(cor_data)){
  tmp_cors <- cor_data[[dset]]
  
  top10 <- sort(abs(tmp_cors[,d1]),decreasing = T)[1:11]
  top10 <- data.frame(gene=names(top10),edges=top10,stringsAsFactors = F) %>% remove_rownames(.)
  top10$related <- top10$gene %in% subset(related,target == d1)$partner
  top10$label <- gsub(" .*","",top10$gene)
  top10$edge_norm[2:11] <-  top10$edges[2:11] / sum(top10$edges[2:11])
  top10$edge_norm[1] <- 1
  top10$annotation <- "non-related"
  top10$annotation[top10$related] <- "related"
  top10$annotation[1] <- "target"
  network_dfs_1[[dset]] <- top10
  
  top10 <- sort(abs(tmp_cors[,d2]),decreasing = T)[1:11]
  top10 <- data.frame(gene=names(top10),edges=top10,stringsAsFactors = F) %>% remove_rownames(.)
  top10$related <- top10$gene %in% subset(related,target == d2)$partner
  top10$label <- gsub(" .*","",top10$gene)
  top10$edge_norm[2:11] <-  top10$edges[2:11] / sum(top10$edges[2:11])
  top10$annotation <- "non-related"
  top10$annotation[top10$related] <- "related"
  top10$annotation[1] <- "target"
  top10$edge_norm[1] <- 1
  
  network_dfs_2[[dset]] <- top10
  
}


node_dfs <- list("RPS21 (6227)"=network_dfs_1,
                 "MED16 (10025)"=network_dfs_2)
link_dfs <- list("RPS21 (6227)"=list(),
                 "MED16 (10025)"=list())

for (target in c(d1,d2)){
  
  #Get top 10 connections for each gene in first degree nodes
  for (dset in names(cor_data)){
    tmp_cors <- cor_data[[dset]]
    tmp_nodes <- node_dfs[[target]][[dset]]
    colnames(tmp_nodes)[1] <- "id"
    node_dfs[[target]][[dset]] <- tmp_nodes
    
    all_nodes <- tmp_nodes$id
    
    #compile links for each node
    edge_list <- list()
    for (codep_id in all_nodes){
      second_top10 <- sort(abs(tmp_cors[,codep_id]),decreasing = T)[1:11]
      second_top10 <- second_top10[names(second_top10) %in% all_nodes]
      tmp_links <- data.frame(from=codep_id,to=names(second_top10),type="abs_pearson",weight=as.numeric(second_top10),stringsAsFactors = F)
      edge_list[[codep_id]] <- tmp_links
    }
    edge_df <- bind_rows(edge_list)
    edge_df %<>% subset(., from != to)
    link_dfs[[target]][[dset]] <- edge_df
  }
  
}

for (target in c(d1,d2)){
  
  #Get top 10 connections for each gene in first degree nodes
  for (dset in names(cor_data)){
    
    nodes <- node_dfs[[target]][[dset]]
    links <- link_dfs[[target]][[dset]]
    
    net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    
    colrs <- c("tomato","gold","gray50")
    V(net)$color <- colrs[as.numeric(factor(V(net)$annotation,levels=c("target","related","non-related")))]
    
    deg <- degree(net, mode="all")
    V(net)$size <- deg*3
    
    E(net)$arrow.size <- .3
    E(net)$edge.color <- "gray80"
    E(net)$width <- E(net)$weight * 5
    
    graph_attr(net, "layout") <- layout_as_star
    
    filename <- paste0("codependency_",dset,"_",gsub(" .*","",target),"_network.pdf")
    pdf(file=file.path("figures",filename),width=5,height=5)
    plot(net,vertex.label=V(net)$label, vertex.label.color="black") 
    dev.off()
    
  }
}


# [1] "layout_as_star"       "layout_components"    "layout_in_circle"    
# [4] "layout_nicely"        "layout_on_grid"       "layout_on_sphere"    
# [7] "layout_randomly"      "layout_with_dh"       "layout_with_drl"     
# [10] "layout_with_fr"       "layout_with_gem"      "layout_with_graphopt"
# [13] "layout_with_kk"       "layout_with_lgl"      "layout_with_mds" 




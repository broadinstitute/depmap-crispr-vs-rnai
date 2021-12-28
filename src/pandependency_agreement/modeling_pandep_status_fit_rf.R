
library(viridis)

# require(randomForest)
# require(pROC)
# require(rpart)
# require(rpart.plot)
# require(randomForestExplainer)

source("src/packages_paths.R")

# Fit a model to classify the genes into shared pandep or CRISPR-specific pandep
# y: df with 1 numeric binary column to indicate class (referred to by position), rownames must match X
# X: df with same rownames as X, each column represents a feature
# tree_explainer_plots: T/F value to indicate whether to plot an example tree and feature depth
fit_rf_model <- function(y,X,tree_explainer_plots){
  
  y_ordered <- rownames(y)[order(y[,1])]
  split_vec <- 1:length(y_ordered) %% 5
  
  imp_list <- list()
  pred_list <- list()
  for (k in 0:5){
    
    if (k < 5){
      test_genes <- y_ordered[split_vec == k]
      train_genes <- y_ordered[split_vec != k]
      
      X_train <- X[train_genes,]
      X_test <- X[test_genes,]
      
      y_train <- y[train_genes,1]
      names(y_train) <- train_genes
      y_train <- factor(y_train,levels=c("0","1"))
      y_test <- y[test_genes,1]
      names(y_test) <- test_genes
      y_test <- factor(y_test,levels=c("0","1"))
      
      rf <- randomForest::randomForest(x=X_train,y=y_train,xtest=X_test,ytest=y_test,importance=T,keep.forest=T,localImp = TRUE)
      pred_list[[as.character(k)]] <- rf$test$votes
      imp_list[[as.character(k)]] <- rf$importance
      
    } else if (tree_explainer_plots & (k == 5)) {
      
      train_genes <- y_ordered[split_vec != k]
      X_train <- X[train_genes,]
      y_train <- y[train_genes,1]
      names(y_train) <- train_genes
      y_train <- factor(y_train,levels=c("0","1"))
      rf <- randomForest::randomForest(x=X_train,y=y_train,importance=T,keep.forest=T,localImp = TRUE)
      
      forest_frame <- randomForestExplainer::min_depth_distribution(rf)
      randomForestExplainer::plot_min_depth_distribution(forest_frame) +
        xlab("") +
        ggtitle("") +
        theme_bw(base_size=11) +
        theme(axis.text=element_text(size=9)) +
        scale_fill_viridis(discrete=TRUE,option = "D")
      ggsave(file.path("figures","pandependency_agreement_model_allfeat_depth_distribution.pdf"),width=6,height=4)
      
      
      model_df <- cbind(y,X)
      single_tree <- rpart::rpart(pandep_group~., data=model_df, method="class")
      pdf(file.path("figures","pandependency_agreement_model_allfeat_example_tree.pdf"),width=4.75,height=3)
      rpart.plot::rpart.plot(single_tree, box.palette="RdBu", shadow.col="gray", nn=TRUE,yesno=2,
                             type=2,
                             clip.right.labs = F,
                             branch=1)
      dev.off()
      
    }
    
  }
  
  y_pred <- do.call(rbind,pred_list)
  
  feat_imp <- Reduce('+', imp_list)
  feat_imp %<>% as.data.frame(.)
  feat_imp$norm <- feat_imp$MeanDecreaseGini / sum(feat_imp$MeanDecreaseGini)
  
  #Calculate ROC AUC
  y_pred <- y_pred[rownames(y),]
  rf.roc <- pROC::roc(y[,1],y_pred[,"1"])
  
  return(list("pred"=y_pred,"obs"=y,"roc"=rf.roc,"auc"= auc(rf.roc),"feat_imp"=feat_imp))
  
}

###################### Get y vector of pandep class ###################### 

#### Get CRISPR and RNAi pandep status for high-confidence dependencies
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

t2 %<>% subset(.,pandep_group %in% c("CRISPR-specific","shared"))
scatter_data <- dplyr::select(t2,entrez_id,pandep_group)

###################### Construct X matrices ###################### 
stats_df <- readRDS(file.path(data_processed,"modeling_pandep_status_predictive_features.rds"))
stats_df <- stats_df[!is.na(stats_df$CRISPR_gene_id),]
stats_df %<>% add_column(.,entrez_id=extract_entrez(stats_df$CRISPR_gene_id),.before=1)
stats_df <- stats_df[,-grep("gene_id",colnames(stats_df))]
dep_columns <- sort(unique(c(grep("CRISPR",colnames(stats_df)),grep("RNAi",colnames(stats_df)))))
stats_dep <- stats_df[,c(1,dep_columns)]
stats_nondep <- stats_df[,-dep_columns]
stats_nondep <- stats_nondep[,-grep("Expression_mult_iso",colnames(stats_nondep))]
stats_nondep$Protein_mult_iso <- stats_nondep$Protein_mult_iso + 0
stats_nondep <- stats_nondep[,-grep("symbol",colnames(stats_nondep))]

full_df <- left_join(scatter_data,stats_nondep,by="entrez_id")
full_df %<>% subset(.,!is.na(Expression_Protein_cor))

# y target: pandep class, diff in CRISPR-RNAi mean, RNAi mean
dups <- subset(full_df,full_df$entrez_id %in% full_df$entrez_id[duplicated(full_df$entrez_id)])
not_dup <- subset(full_df,!(full_df$entrez_id %in% full_df$entrez_id[duplicated(full_df$entrez_id)]))

#when using protein features
dups <- dups[order(dups$Protein_missing),]
dups <- dups[!duplicated(dups$entrez_id),]

full_df <- rbind(dups,not_dup)
full_df %<>% remove_rownames(.) %>% column_to_rownames(.,var="entrez_id")

y <- full_df[,"pandep_group",drop=F]
y$pandep_group <- y$pandep_group == "CRISPR-specific"
y$pandep_group <- y$pandep_group + 0

print(paste0("CRISPR-specific pandep targets: ",sum(y$pandep_group)))
print(paste0("shared pandep targets: ",nrow(y) - sum(y$pandep_group)))

X <- full_df[,!grepl("pandep_group",colnames(full_df))]

rnaseq_feats <- c("Expression_u","Expression_q1","Expression_q2","Expression_sd","Expression_missing")
prot_feats <- c("Protein_u","Protein_q1","Protein_q2","Protein_sd","Protein_missing","Protein_mult_iso")

feat_sets <- list(
  "mRNA"=X[,rnaseq_feats],
  "Protein"=X[,prot_feats],
  "mRNA+prot+cor"=X[,c(rnaseq_feats,prot_feats,"Expression_Protein_cor")])

###################### Fit models ###################### 

rf_models <- list()
for (set_name in names(feat_sets)){
  if (set_name == "mRNA+prot+cor"){
    rf_models[[set_name]] <- fit_rf_model(y=y,X=feat_sets[[set_name]],tree_explainer_plots=T)
  } else {
    rf_models[[set_name]] <- fit_rf_model(y=y,X=feat_sets[[set_name]],tree_explainer_plots=F)
  }
}

saveRDS(rf_models,file.path(data_processed,"rf_pandep_class_models_wFeats.rds"))


library(Rcpp)
library(RcppEigen)

source("src/packages_paths.R")

sourceCpp("src/codependency/matrix_mult_functions.cpp")
source('src/codependency/SNF_mod.R')

eigen_mult <- function(A,B){
  return(eigenMapMatMult(eigenMapMatMult(A, B), t(A)))
}

hgnc <- fread(file.path(data_raw,"hgnc-complete-set-2019.csv"))
hgnc$entrez_id %<>% as.character(.)

combinations <- list()

#### Get high-confidence dependencies with at least 3 dependent cell lines using CRISPR or RNAi ####
t1 <- fread(file.path("tables","Supplemental-Table-1.csv"))
t1$entrez_id %<>% as.character(.)
t1 %<>% subset(.,high_confidence)

t2 <- fread(file.path("tables","Supplemental-Table-2.csv"))
t2$entrez_id %<>% as.character(.)
t2 %<>% subset(.,entrez_id %in% t1$entrez_id)

t2 %<>% subset(.,(`CRISPR_numDeps(50-100)` >= 3) | (`RNAi_numDeps(50-100)` >= 3))

#Matched genetic screens
genetic_data <- list("crispr"="gene-effect-scaled-crispr-matched.csv",
                     "rnai"="gene-effect-scaled-rnai-matched.csv"
)
genetic_data <- lapply(genetic_data,function(x){fread(file.path(data_raw,x)) %>% column_to_rownames(.,var="Row.name") %>% as.matrix(.)})

avana <- genetic_data[["crispr"]]
d2_comb <- genetic_data[["rnai"]]

keep_ids <- colnames(avana)[extract_entrez(colnames(avana)) %in% t2$entrez_id]

avana <- avana[,keep_ids]
d2 <- d2_comb[,keep_ids]

test_cors <- cor(avana,use="pairwise.complete")
stopifnot(!any(is.na(test_cors)))
saveRDS(test_cors,file.path(data_processed,"codependency_CRISPR_pearson_baseline.rds"))
test_cors <- cor(d2,use="pairwise.complete")
stopifnot(!any(is.na(test_cors)))
saveRDS(test_cors,file.path(data_processed,"codependency_RNAi_pearson_baseline.rds"))

combinations[["CRISPR-RNAi"]] <- list("CRISPR"=avana,
                                      "RNAi"=d2)


####################
##### Run SNF ######
####################

for (combo in names(combinations)){
  print(paste0("Combination: ",combo))
  datasets <- combinations[[combo]]
  
  datasets <- lapply(datasets,function(x){t(x)})
  print("computing W")
  W_all <- pearsonWeightMatrix(datasets)
  stopifnot(!any(is.na(W_all[[1]])))
  stopifnot(!any(is.na(W_all[[2]])))
  for (d in names(W_all)){
    saveRDS(W_all[[d]],file=file.path(data_processed,paste0(combo,"_",d,"-weight-matrix.rds")))
  }
  
  
  print("computing P")
  P_all <- lapply(W_all, function(X) normalize2(X))
  K = ceiling(log(ncol(W_all[[1]])))
  print(paste0("k: ",K))
  print("computing S")
  S_all <- lapply(W_all, function(X) knearestNeighborABS(X,K))
  print("running SNF2")
  
  #SNF2 <- function (P_all, S_all, t = 20, norm_func){
  t = 20 
  norm_func = normalize2
  
  LW = length(P_all)
  nextP <- vector("list", LW)
  
  for (i in 1:t) { # for each iteration
    print(paste0("Iteration: ",i))
    for (j in 1:LW) { #for each dataset
      
      sumPJ = matrix(0, dim(P_all[[j]])[1], dim(P_all[[j]])[2])
      for (k in 1:LW) {
        if (k != j) {
          sumPJ = sumPJ + P_all[[k]]
        }
      }
      
      #Original matrix multiplication -- too slow
      #nextP[[j]] = S_all[[j]] %*% (sumPJ/(LW - 1)) %*% t(S_all[[j]])
      
      nextP[[j]] <- eigen_mult(A=S_all[[j]],B=(sumPJ/(LW - 1)))
      
    } #end dataset
    
    #apply normalization
    P_all <- lapply(nextP, function(X) norm_func(X))
    
    
  } #end iteration
  
  #Average the P matrices to get W
  W = matrix(0, nrow(P_all[[1]]), ncol(P_all[[1]]))
  for (i in 1:LW) {
    W = W + P_all[[i]]
  }
  W = W/LW
  
  #normalize
  W = norm_func(W)
  W = (W + t(W) + diag(nrow(W)))/2
  
  colnames(W) <- colnames(W_all[[1]])
  rownames(W) <- rownames(W_all[[1]])
  
  print("writing output")
  saveRDS(W,file=file.path(data_processed,paste0("codependency-",combo,"-SNF.rds")))
  
  
}


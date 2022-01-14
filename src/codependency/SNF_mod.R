library(matrixStats)
library(spatstat)
library(igraph)
forceGraph <- function(S,groups,labels){
  #S <- normalize2(S)
  #S <- (S + t(S)) / 2
  diag(S) <- 0
  g <- graph.adjacency(S,mode="upper",weighted=TRUE)
  V(g)$label.cex <- .5
  coords <- layout_with_kk(g, weights=E(g)$weight)
  plot(g, layout=coords, vertex.color=groups, vertex.size=10, vertex.label=labels)
}

first.degree.neighbors <- function(Wg,genes){
  
  genes <- genes[genes %in% rownames(Wg)]
  
  Wg_logical <- Wg > 0
  
  neighbors <- list()
  for (g in genes){
    neighbors[[g]] <- colnames(Wg_logical)[Wg_logical[g,]]
  }
  
  return(neighbors)
}

#Takes list of raw datasets with samples as rows
#Converts variables to zscore and computes pearson correlation matrix
pearsonWeightMatrix <- function(datasets){
  datasets <- lapply(datasets,function(X) standardNormalization2(X))
  Wall <- lapply(datasets,function(X) {
    if (any(is.na(X))){
      X <- cor(t(X),method="pearson",use="pairwise.complete.obs")
    } else {
      X <- cor(t(X),method="pearson")
    }
    
    #X[X < 0] <- 0
    X  <- (X + 1)/2
    return(X)
  })
}

standardNormalization2 <- function (x) 
{
  x = as.matrix(x)
  m = colMeans(x,na.rm=T)
  s = colSds(x,na.rm=T)
  s[s == 0] = 1
  xNorm = t((t(x) - m)/s)
  return(xNorm)
}

knearestNeighbor <- function(xx,KK) {
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,])
  }
  
  return(normalize(A))
}

knearestNeighborABS <- function(xx,KK) {
  zero <- function(x) {
    s = sort(abs(x), index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,])
  }
  
  return(normalize(A))
}

dominateset2 <- function(xx,KK=20,norm_func) {
  #if ((length(xx)-KK) <= 0){
  #  return(xx)
  #}
  zero <- function(x) {
      s = sort(x, index.return=TRUE)
      x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  #normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,])
  }
  
  return(norm_func(A))
}

CS_Prediction <- function(W,Y0,method,alpha=0.9){
  
  P= W/rowSums(W)
  
  if(method == 0){
    Y <- (1-alpha) * solve(diag(dim(P)[1]) - alpha * P) %*% Y0
  } else {
    has_label <- which(rowSums(Y0) != 0)
    Y <- Y0
    for (i in 1:1000){
      Y <- P %*% Y
      Y[has_label,] <- Y0[has_label,]
    }
  }
  return(Y)
}

normalize <- function(X) X/rowSums(X)

normalize2 <- function(X) {
  diag(X) <- 0
  X <- X/(2*rowSums(X))
  diag(X) <- .5
  return(X)
}

SNF2 <- function (P_all, S_all, t = 20, norm_func){
  
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
      nextP[[j]] = S_all[[j]] %*% (sumPJ/(LW - 1)) %*% t(S_all[[j]])
      
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
  
  rownames(W) <- rownames(S_all[[1]])
  colnames(W) <- colnames(S_all[[1]])
  return(W)
}

groupPredict2 <- function(W,groups,method=1){
  
  #Transform label information into binary matrix
  Y0 <- matrix(0,nrow(W), max(groups))
  for (i in 1:length(groups)) {
    Y0[i,groups[i]] <- 1
  }
  
  #Apply method
  Y <- CS_Prediction(W,Y0,method)
  
  #Discretize result predictions back into labels
  newgroups <- rep(0,nrow(Y))
  for (i in 1:nrow(Y)){ 
    newgroups[i] <- which(Y[i,] == max(Y[i,]))
  }
  
  return (newgroups)
}
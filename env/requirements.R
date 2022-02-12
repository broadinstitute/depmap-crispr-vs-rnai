
install.packages('tidyverse') #1.3.1
install.packages('magrittr') #2.0.1
install.packages('data.table') #1.14.2

install.packages('Rcpp') #1.0.7
install.packages('RcppEigen') #0.3.3.9.1

install.packages('matrixStats') #0.61.0
install.packages('spatstat') #2.3-0
install.packages('MASS') #7.3-54
install.packages('sn') #2.0.1
install.packages('e1071') #1.7-9

install.packages('igraph') #1.2.11

install.packages('plyr') #1.8.6

install.packages('purrr') #0.3.4

install.packages('reshape') #0.8.8
install.packages('reshape2') #1.4.4

install.packages('hexbin') #1.28.2

install.packages('cowplot') #1.1.1
install.packages('ggalluvial') #0.12.3
install.packages('ggrepel') #0.9.1
install.packages('gridExtra') #2.3
install.packages('ggmosaic')
install.packages('eulerr') #6.1.1
install.packages('viridis') #0.6.2
install.packages('ggsci') #2.9
install.packages('ggridges') #0.5.3
install.packages('ggbeeswarm') #0.6.0
install.packages('ggExtra') #0.9

install.packages('randomForest') #4.6-14
install.packages('pROC') #1.18.0
install.packages('rpart.plot') #3.1.0
install.packages('randomForestExplainer') #0.10.1


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea",update=F)

#For WGCNA
BiocManager::install("impute",update=F) #1.68.0
BiocManager::install("preprocessCore",update=F) #1.56.0
BiocManager::install("GO.db",update=F) 
install.packages('WGCNA') #1.70-3

BiocManager::install("GSVA",update=F) #1.42.0

#Already included
# library(tools) #4.1.2
# library(stringr) #1.4.0
# require(grid) #4.1.2
# require(scales) #1.1.1
# library(jsonlite) #1.7.2
# require(rpart) #4.1-15

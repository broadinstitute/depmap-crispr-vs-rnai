
# require(purrr)
# library(readr)

source("src/packages_paths.R")

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
out_file <- args[2]
data <- fread(data_file) %>% column_to_rownames(.,var="Row.name")

get_bins <- function(gene){
  #hist(gene, breaks=c(0.0,0.25,0.5,0.75,1.0),plot = FALSE)$counts
  hist(gene, breaks=c(0.0,0.5,0.9,1),plot = FALSE)$counts
}

res <- purrr::map_df(data,get_bins) %>% 
  t %>% 
  as.data.frame %>%
  #set_colnames(.,paste0("numDeps(",c("0-25","25-50","50-75","75-100"),")")) %>% 
  set_colnames(.,paste0("numDeps(",c("0-50","50-90","90-100"),")")) %>% 
  rownames_to_column(.,var="Row.name")

res$`numDeps(50-100)` <- rowSums(select(res,one_of(c("numDeps(50-90)","numDeps(90-100)"))))
res$depCL_frac <- res$`numDeps(50-100)` / (res$`numDeps(0-50)` + res$`numDeps(50-100)`)

# res %<>% select(.,`Row.name`,depCL_frac,strong_depCL_count=`numDeps(90-100)`)

write_csv(res,path=out_file)

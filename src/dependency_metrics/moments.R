
# require(purrr)
# require(e1071)
# library(readr)

source("src/packages_paths.R")

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
out_file <- args[2]
data <- fread(data_file) %>% column_to_rownames(.,var="Row.name")

fun <- function(f) purrr::pmap_dbl(list(x = data, na.rm = TRUE), f)
param <- list(list(mean), 
              list(var),
              list(e1071::skewness))

res <- purrr::invoke_map(.f = fun, .x = param)

df <- data.frame(Mean=res[[1]],Variance=res[[2]],Skewness=res[[3]])
df %<>% rownames_to_column(.,var="Row.name")

write_csv(df,path=out_file)
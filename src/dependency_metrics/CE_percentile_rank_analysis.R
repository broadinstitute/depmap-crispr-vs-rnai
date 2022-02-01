
source("src/packages_paths.R")

require(plyr)

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
out_file <- args[2]
out_file_list <- args[3]
data <- fread(data_file) %>% column_to_rownames(.,var="Row.name")
data <- t(data)

rank_data <- plyr::aaply(data, .margins=2, rank, ties.method="min")
rank_data %<>% t(.)
rank_data[is.na(data)] <- NA
rank_data %<>% t(.)
rank_data <- rank_data / colSums(!is.na(data))
percentile <- plyr::aaply(rank_data, .margins=2, quantile, probs=.9, na.rm=T)

dens <- density(percentile)
df <- data.frame(x=dens$x, y=dens$y)
df_range <- subset(df,df$x > .1 & df$x < .9)
threshold <- df_range$x[which.min(df_range$y)]

ce_list <- data.frame(`Row.name`=names(percentile),CE_percentile=percentile)
ce_list$Common_Essential <- ce_list$CE_percentile <= threshold
write_csv(ce_list,path=out_file)

ce_achilles_format <- subset(ce_list, Common_Essential) %>% select(.,gene=`Row.name`)
write_csv(ce_achilles_format,path=out_file_list)

# require(matrixStats)

source("src/packages_paths.R")

args <- commandArgs(trailingOnly = TRUE)
pr_file <- args[1]
exp_file <- args[2]
out_file <- args[3]

pr <- fread(pr_file) %>% column_to_rownames(.,var="Row.name")

var_df <- data.frame(CDS_ID=colnames(pr),dep_var=matrixStats::colVars(as.matrix(pr),na.rm=T))
var_df %<>% add_column(.,entrez_id=extract_entrez(var_df$CDS_ID),.before=1)

exp <- fread(exp_file) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
cls <- rownames(pr)
genes <- colnames(pr)
exp <- exp[cls,colnames(exp) %in% genes]
d <- density(unlist(exp))
d_df <- data.frame(x=d$x,y=d$y)
d_df %<>% subset(.,(d_df$x > 0) & (d_df$x < 5))
exp_cutoff <- d_df$x[which.min(d_df$y)]
nonexp_frac <- colSums(exp < exp_cutoff) / nrow(exp)
nonexp <- data.frame(genes=colnames(exp),frac=nonexp_frac)
nonexp_genes <- nonexp$genes[nonexp$frac == 1]

nonexp_var <- subset(var_df,CDS_ID %in% nonexp_genes)
nonexp_var$var_percentile <- percentile_rank(nonexp_var$dep_var)

var_thresh <- min(subset(nonexp_var,var_percentile >= .99)$dep_var)

var_df$high_variance <- var_df$dep_var > var_thresh

var_df$nonexpressed <- var_df$CDS_ID %in% nonexp_genes
var_df %<>% rename(.,`Row.name`="CDS_ID")

write_csv(var_df,file=out_file)


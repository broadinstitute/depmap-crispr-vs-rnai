
source("src/packages_paths.R")

library(sn)
# require(MASS)

g <- "WRN (7486)"

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-avana.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
rnai_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
crispr_gs <- crispr_gs[rownames(crispr_gs) %in% rownames(rnai_gs),]

crispr_vec <- crispr_gs[,g]
rnai_vec <- rnai_gs[,g]
rnai_vec <- rnai_vec[!is.na(rnai_vec)]

plot_LRT <- function(vec,gene){
  skewt_fit <- data.frame(data = vec) %>% sn::selm(data ~ 1, family = "ST", data = .)
  y <- skewt_fit@residuals.dp + skewt_fit@fitted.values.dp
  dp0 <- skewt_fit@param$dp
  n <- skewt_fit@size["n.obs"]
  h <- hist(y, plot = FALSE)
  extr <- extendrange(x = h$breaks)
  x.pts <- seq(max(extr), min(extr), length = 501)
  d.fn <- get(paste("d", tolower(skewt_fit@family), sep = ""), 
              inherits = TRUE)
  pdf <- d.fn(x.pts, dp = dp0)
  norm_fit <- MASS::fitdistr(vec, 'normal')
  pdf_norm <- dnorm(x.pts, mean=norm_fit$estimate["mean"], sd=norm_fit$estimate["sd"])
  
  hist_df <- data.frame(values=y)
  skew_df <- data.frame(x=x.pts,y=pdf,type="skewed-t",stringsAsFactors = F)
  norm_df <- data.frame(x=x.pts, y=pdf_norm,type="gaussian",stringsAsFactors = F)
  lines_df <- rbind(skew_df,norm_df)
  
  lines_pal <- c("skewed-t"="#3A53A4","gaussian"="#6ABD45")
  
  p <- ggplot(data=hist_df,aes(values)) +
    geom_histogram(aes(y=..density..),bins = 25,fill="gray95",col="gray60") +
    geom_rug() +
    geom_line(data=lines_df,aes(x=x,y=y,color=type),size=1) +
    scale_color_manual(values=lines_pal) +
    theme_classic(base_size=11) +
    xlab(paste0(gene," Gene Score")) +
    ylab("Probability Density") +
    theme(legend.position = c(0.25, 0.6)) +
    guides(color=guide_legend(title="Distribution"))
  
  return(p)
}

rnai_plot <- plot_LRT(rnai_vec,g)
ggsave(plot=rnai_plot,file.path("figures","gene_deps_profiles_WRN_LRT_RNAi_illustration.pdf"),height=3,width=3.5)

crispr_plot <- plot_LRT(crispr_vec,g)
ggsave(plot=crispr_plot,file.path("figures","gene_deps_profiles_WRN_LRT_CRISPR_illustration.pdf"),height=3,width=3.5)




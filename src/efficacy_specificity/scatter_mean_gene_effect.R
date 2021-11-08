library(tidyverse)
library(magrittr)
library(data.table)

#requires(wesanderson)
#requires(hexbin)

data_raw <- file.path("data","raw")
data_processed <- file.path("data","processed")

crispr_gs <- fread(file.path(data_raw,"gene-effect-scaled-crispr-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)
rnai_gs <- fread(file.path(data_raw,"gene-effect-scaled-rnai-matched.csv")) %>% column_to_rownames(.,var="V1") %>% as.matrix(.)

mean_per_gene <- data.frame(CRISPR=colMeans(crispr_gs,na.rm=T),RNAi=colMeans(rnai_gs,na.rm=T))

pal <- wesanderson::wes_palette("Royal1", 200, type = "continuous")
pal <- c("#899DA4","#C93312","#DC863B","#FAEFD1")

p1 <- ggplot(mean_per_gene) +
  geom_hex(aes(x=CRISPR,y=RNAi,colour = ..count..),bins=100) +
  scale_color_gradientn(colors=pal,name = "count", trans = "log",
                        breaks = 10^(0:6)) +
  scale_fill_gradientn(colors=pal,name = "count", trans = "log",
                       breaks = 10^(0:6)) +
  xlab("CRISPR Gene Effect") +
  ylab("RNAi Gene Effect") +
  geom_abline(slope=1,intercept=0,linetype=2,color="black",size=.5) +
  theme_classic(base_size=11) +
  theme(legend.position="bottom")
ggsave(plot=p1,file.path("figures","efficacy_specificity_scatter_mean_gene_effect_legend.pdf"),height=4.5,width=4.5)


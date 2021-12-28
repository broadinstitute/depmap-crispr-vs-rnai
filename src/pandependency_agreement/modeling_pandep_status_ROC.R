
source("src/packages_paths.R")

pandep_models <- readRDS(file.path(data_processed,"rf_pandep_class_models_wFeats.rds"))

plot_df <- list()
for (m in names(pandep_models)){
  tmp_roc <- pandep_models[[m]]$roc
  plot_df[[m]] <- data.frame("model"=m,"Specificity"=tmp_roc$specificities,"Sensitivity"=tmp_roc$sensitivities,stringsAsFactors = F)
}
plot_df <- bind_rows(plot_df)
plot_df$model <- factor(plot_df$model,levels=c("mRNA","Protein","mRNA+prot+cor"))

mypal <- c("#FD6467","#5B1A18","#D67236")
names(mypal) <- c("mRNA","mRNA+prot+cor","Protein")
ggplot(plot_df, aes(x=Specificity, y=Sensitivity,color=model)) +
  geom_line(size=1) +
  geom_abline(intercept = 1,slope=1,col="gray",linetype="dashed") +
  scale_color_manual(values=mypal) +
  scale_x_reverse() +
  theme_classic(base_size=11) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.75,.25))
ggsave(file.path("figures","pandependency_agreement_model_ROC.pdf"),height=2.5,width=2.2)

       
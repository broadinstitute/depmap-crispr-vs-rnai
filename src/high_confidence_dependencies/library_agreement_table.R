
source("src/packages_paths.R")

results <- fread(file.path(data_processed,"library_agreement.csv"))
stable <- dplyr::select(results,symbol,entrez_id)

#### CRISPR

stable$CRISPR_agreement <- results$CRISPR_agreement

stable %<>% add_column(.,CRISPR_shared_target=results$CRISPR_avana_target & results$CRISPR_ky_target)
stable$CRISPR_shared_target[is.na(stable$CRISPR_shared_target)] <- F

stable %<>% add_column(.,CRISPR_shared_pandep=results$CRISPR_avana_pandep & results$CRISPR_ky_pandep)

stable$CRISPR_cor_rank <- results$CRISPR_cor_rank

stable %<>% add_column(.,CRISPR_shared_nondep=results$CRISPR_avana_nodep & results$CRISPR_ky_nodep)

stable$CRISPR_no_shared_depCLs <- results$CRISPR_no_shared_depCLs

#### RNAi

stable$RNAi_agreement <- results$RNAi_agreement

stable %<>% add_column(.,RNAi_shared_target=results$RNAi_achilles_target & results$RNAi_drive_target)
stable$RNAi_shared_target[is.na(stable$RNAi_shared_target)] <- F

stable %<>% add_column(.,RNAi_shared_pandep=results$RNAi_achilles_pandep & results$RNAi_drive_pandep)

stable$RNAi_cor_rank <- results$RNAi_cor_rank

stable %<>% add_column(.,RNAi_shared_nondep=results$RNAi_achilles_nodep & results$RNAi_drive_nodep)

stable$RNAi_no_shared_depCLs <- results$RNAi_no_shared_depCLs

#### High-confidence dependencies
stable %<>% add_column(.,high_confidence=stable$CRISPR_agreement & stable$RNAi_agreement,.before=3)

write_csv(stable,file.path("tables","Supplemental-Table-1.csv"))

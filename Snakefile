
rule get_figshare_data:
	output:
		"data/raw/hgnc-complete-set.csv",
		"data/raw/lfc-unscaled-rnai-achilles.csv",
		"data/raw/lfc-unscaled-rnai-drive.csv",
		"data/raw/lfc-unscaled-crispr-avana.csv",
		"data/raw/lfc-unscaled-crispr-ky.csv",
		"data/raw/reagent-to-gene-map-sgrna.csv",
		"data/raw/reagent-to-gene-map-shrna.csv",
		"data/raw/control-essential-genes-core.csv",
		"data/raw/control-essential-genes-unbiased.csv",
		"data/raw/control-nonessential-genes.csv",
		"data/raw/gene-effect-unscaled-crispr-avana.csv",
		"data/raw/gene-effect-unscaled-crispr-ky.csv",
		"data/raw/gene-effect-unscaled-rnai-drive.csv",
		"data/raw/gene-effect-unscaled-rnai-achilles.csv",
		"data/raw/depmap-omics-expression-rnaseq-tpm.csv"
	shell:
		"Rscript src/figshare_downloader.R"

rule reagent_level_SSMD:
	input:
		"data/raw/lfc-unscaled-rnai-achilles.csv",
		"data/raw/lfc-unscaled-rnai-drive.csv",
		"data/raw/lfc-unscaled-crispr-avana.csv",
		"data/raw/lfc-unscaled-crispr-ky.csv",
		"data/raw/reagent-to-gene-map-sgrna.csv",
		"data/raw/reagent-to-gene-map-shrna.csv",
		"data/raw/control-essential-genes-core.csv",
		"data/raw/control-nonessential-genes.csv"
	output:
		"figures/selecting_datasets_reagents_SSMD_per_dataset.pdf",
		"figures/selecting_datasets_reagents_SSMD_CL_example.pdf"

	shell:
		"Rscript src/selecting_datasets/unprocessed_reagents/SSMD_across_CLs.R"

rule reagent_correlation_per_gene:
	input:
		"data/raw/hgnc-complete-set.csv",
		"data/raw/lfc-unscaled-rnai-achilles.csv",
		"data/raw/lfc-unscaled-rnai-drive.csv",
		"data/raw/lfc-unscaled-crispr-avana.csv",
		"data/raw/lfc-unscaled-crispr-ky.csv",
		"data/raw/reagent-to-gene-map-sgrna.csv",
		"data/raw/reagent-to-gene-map-shrna.csv"
	output:
		"figures/selecting_datasets_reagents_pairwise_cors_boxplot.pdf"
	shell:
		"Rscript src/selecting_datasets/unprocessed_reagents/pairwise_reagent_cor.R"

rule reagent_gene_mean:
	input:
		"data/raw/hgnc-complete-set.csv",
		"data/raw/lfc-unscaled-rnai-achilles.csv",
		"data/raw/lfc-unscaled-rnai-drive.csv",
		"data/raw/lfc-unscaled-crispr-avana.csv",
		"data/raw/lfc-unscaled-crispr-ky.csv",
		"data/raw/reagent-to-gene-map-sgrna.csv",
		"data/raw/reagent-to-gene-map-shrna.csv"
	output:
		"data/processed/mean_reagent_lfc.rds"
	shell:
		"Rscript src/selecting_datasets/processed_reagents/mean_reagent_lfc.R"

rule gene_effect_AUC:
	input:
		"data/raw/hgnc-complete-set.csv",
		"data/processed/mean_reagent_lfc.rds",
		"data/raw/gene-effect-unscaled-crispr-avana.csv",
		"data/raw/gene-effect-unscaled-crispr-ky.csv",
		"data/raw/gene-effect-unscaled-rnai-drive.csv",
		"data/raw/gene-effect-unscaled-rnai-achilles.csv",
		"data/raw/control-essential-genes-unbiased.csv",
		"data/raw/depmap-omics-expression-rnaseq-tpm.csv"
	output:
		"figures/selecting_datasets_processed_vs_unproc_AUC.pdf"
	shell:
		"Rscript src/selecting_datasets/processed_reagents/AUC_per_cell_line.R"

rule dataset_correlation_per_gene:
	input:
		"data/raw/hgnc-complete-set.csv",
		"data/processed/mean_reagent_lfc.rds",
		"data/raw/gene-effect-unscaled-crispr-avana.csv",
		"data/raw/gene-effect-unscaled-crispr-ky.csv",
		"data/raw/gene-effect-unscaled-rnai-drive.csv",
		"data/raw/gene-effect-unscaled-rnai-achilles.csv",
		"data/raw/control-nonessential-genes.csv",
		"data/raw/control-essential-genes-core.csv"
	output:
		"figures/selecting_datasets_proc_vs_unproc_dataset_cor_boxplot.pdf",
		"data/processed/processed_unprocessed_scaled_gene_effects.rds"
	shell:
		"Rscript src/selecting_datasets/processed_reagents/gene_cor_byVar.R"

rule dataset_correlation_gene_recall:
	input:
		"data/raw/hgnc-complete-set.csv",
		"data/processed/processed_unprocessed_scaled_gene_effects.rds"
	output:
		"figures/selecting_datasets_proc_vs_unprc_topcor_recall.pdf"
	shell:
		"Rscript src/selecting_datasets/processed_reagents/numgene_topcor_recall.R"

rule pandependency_per_library:
	input:
		"data/raw/gene-effect-scaled-crispr-avana.csv",
		"data/raw/gene-effect-scaled-crispr-ky.csv",
		"data/raw/gene-effect-scaled-rnai-achilles.csv",
		"data/raw/gene-effect-scaled-rnai-drive.csv"
	output:
		"data/processed/multilib_ce_percentile_results.csv",
		"figures/pandependency_90th_percentile_ranks_CRISPR.pdf",
		"figures/pandependency_90th_percentile_ranks_RNAi.pdf"
	shell:
		"Rscript src/high_confidence_dependencies/ce_percentile_ranks_and_hists.R"

rule pandependency_per_library_benchmarks:
	input:
		"data/raw/hgnc-complete-set.csv",
		"data/raw/control-essential-genes-core.csv",
		"data/raw/control-nonessential-genes.csv",
		"data/raw/control-essential-genes-CEGv2.csv",
		"data/raw/depmap-omics-expression-rnaseq-tpm.csv",
		"data/processed/multilib_ce_percentile_results.csv"
	output:
		"figures/pandependency_90th_percentile_ranks_RNAi_benchmark.pdf",
		"figures/pandependency_90th_percentile_ranks_CRISPR_benchmark.pdf"
	shell:
		"Rscript src/high_confidence_dependencies/pos_neg_cntrl_ce_percentile.R"

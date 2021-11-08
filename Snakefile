
include: "selecting_datasets.snake"
include: "td_metrics.snake"
include: "high_conf_deps.snake"
include: "efficacy_specificty.snake"

rule get_figshare_data:
	output:
		"data/raw/hgnc-complete-set.csv",
		"data/raw/collection.zip"
	shell:
		"Rscript src/figshare_downloader.R"

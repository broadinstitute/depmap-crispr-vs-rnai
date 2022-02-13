docker run -v $(pwd):/tmp/pipeline \
-it depmap-crispr-vs-rnai:latest

snakemake -s ensemble_prediction.snake --configfile snake_config.json --cores 1


Set jobs and variance count to proper numbers. They have been set low for testing
fusion removed for internal. Need to be recreated


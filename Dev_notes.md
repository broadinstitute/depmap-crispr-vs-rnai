Update the dockerfile present for local run
Create docker file for snakemake for remote execution

### efficacy_specificity.snake
The snakemake DAG is not proper. The input of the first rule does not use other rules so a lot of the other runs are never called

### gene_dep_profiles.snake
Similar as efficacy_specificity.snake

### pandependency_agreement.snake
These inputs are missing:
tables/Supplemental-Table-2.csv
tables/Supplemental-Table-1.csv

### td_metrics.snake
Several input files missing
Missing input files for rule gene_deps_table:
pandependency-rnai-matched.csv
pandependency-rnai-drive.csv
LRT-rnai-drive.csv
pandependency-crispr-avana.csv
dependency-probability-rnai-drive.csv
LRT-rnai-achilles.csv
LRT-rnai-matched.csv
dependency-probability-rnai-matched.csv
dependency-probability-crispr-avana.csv
pandependency-crispr-matched.csv
pandependency-rnai-achilles.csv
pandependency-crispr-ky.csv
LRT-crispr-matched.csv
dependency-probability-rnai-achilles.csv
LRT-crispr-ky.csv
dependency-probability-crispr-ky.csv
LRT-crispr-avana.csv
dependency-probability-crispr-matched.csv

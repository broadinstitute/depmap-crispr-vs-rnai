
<!-----

Conversion notes:

* Docs to Markdown version 1.0β33
* Mon Jan 31 2022 02:16:32 GMT-0800 (PST)
* Source doc: Untitled document
----->


# Comparing DepMap CRISPR and RNAi datasets


## Motivation

Reproduce the analyses, figures, and tables presented in Krill-Burger et al., “Partial gene suppression improves identification of cancer vulnerabilities when CRISPR-Cas9 knockout is pan-lethal”. 


## Repository structure



* Organization of R and Python scripts into Snakefiles
    * Snakemake is a workflow language to orchestrate a series of tasks (rules) where each task has defined inputs, outputs, and a shell command to execute (indicates the script being run). 
    * Each Snakefile (designated by the ‘.snake’ suffix) corresponds to a section of related results
        * selecting_datasets.snake: Screen quality comparisons across reagent-level and processed CERES or DEMETER2 gene effects for all large genetic perturbation datasets (Supplemental Fig. 1, 2)
        * high_conf_deps.snake: Defining a high-confidence dependency set (Supplemental Fig. 3, 4, Supplemental Table 1)
        * efficacy_specificity.snake: Benchmarking efficacy and specificity (Fig. 1a-c, Supplemental Fig. 5)
        * td_metrics.snake: Compute metrics used bin genes into the dependency classes, such as pan-dependency, high-variance, strongly selective, or non-dependency (Supplemental Table 2)
        * gene_dep_profiles.snake: Analysis of gene dependency classes between perturbation types (Fig. 1d-f, Supplemental Fig. 6-8,12-13)
        * pandependency_agreement.snake: Analysis of pan-dependencies that are shared or distinct between perturbation types (Fig. 1g-h, Supplemental Fig. 9-11)
        * ensemble_prediction.snake: Ensemble prediction pipeline to fit random forest models to all gene effect profiles using multi-omics and cell line annotation predictive features
        * predictive_markers.snake: Summary of ensemble prediction results and CYCLOPS analysis (Fig. 2, Supplemental Table 3)
        * drug_response.snake: Drug-gene target associations (Fig. 3, Supplemental Fig. 16)



## Installing dependencies
    
Required R packages are in `env/general/requirements.R`. However, we suggest using Docker for reproducibility. The image can be built from the Dockerfile:
```
$cd ./env/general
$docker build -t depmap-crispr-vs-rnai .
```



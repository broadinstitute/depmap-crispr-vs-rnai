
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
    * [Snakemake](https://snakemake.readthedocs.io/) is a workflow language to orchestrate a series of tasks (rules) where each task has defined inputs, outputs, and a shell command to execute (indicates the script being run).
    * Each Snakefile (designated by the ‘.snake’ suffix) corresponds to a section of related results and can be found in the root directory of this code repository
        * **Supplemental Fig. 1, 2**: Screen quality comparisons across reagent-level and processed CERES or DEMETER2 gene effects for all large genetic perturbation datasets 
            * selecting_datasets.snake
        * **Supplemental Fig. 3, 4, Supplemental Table 1**: Defining a high-confidence dependency set
            * high_conf_deps.snake: 
        * **Fig. 1a-c, Supplemental Fig. 5**: Benchmarking efficacy and specificity 
            * efficacy_specificity.snake
        * **Supplemental Table 2**: Compute metrics used bin genes into the dependency classes, such as pan-dependency, high-variance, strongly selective, or non-dependency 
            * td_metrics.snake
        * **Fig. 1d-f, Supplemental Fig. 6-8,12-13**: Analysis of gene dependency classes between perturbation types 
            * gene_dep_profiles.snake
        * **Fig. 1g-h, Supplemental Fig. 9-11**: Analysis of pan-dependencies that are shared or distinct between perturbation types 
            * pandependency_agreement.snake
        * Ensemble prediction pipeline to fit random forest models to all gene effect profiles using multi-omics and cell line annotation predictive features
            * ensemble_prediction.snake
        * **Fig. 2, Supplemental Table 3**: Summary of ensemble prediction results and CYCLOPS analysis 
            * predictive_markers.snake
        * **Fig. 3, Supplemental Fig. 16**: Drug-gene target associations 
            * drug_response.snake
        * **Fig. 4**: Functional relationships identified by co-dependency network
            * codependency.snake



## Installing dependencies
    
Required R packages are in `env/requirements.R`. We suggest using [Docker](https://docs.docker.com/) for reproducibility. A Docker image with the required R packages and Snakemake can be built from the Dockerfile:

```
$docker build -f env/Dockerfile -t jkrillbu/depmap-crispr-vs-rnai:1 .
```

Alternatively, a pre-built Docker image can be downloaded from DockerHub with the following command
```
$docker pull jkrillbu/depmap-crispr-vs-rnai:1
```

One way to use the Docker image locally is to create an interacive container. Executing the following commands will create a link between `/tmp/pipeline` within the container and the root of the code repository:
```
$git clone https://github.com/broadinstitute/depmap-crispr-vs-rnai.git
$cd depmap-crispr-vs-rnai
$docker run -v $(pwd):/tmp/pipeline -it jkrillbu/depmap-crispr-vs-rnai:1
```

Now individual scripts or snakefiles can be executed. For example, the `src/figshare_downloader.R` R script is run with the following command
```
root@id:/tmp/pipeline# Rscript src/figshare_downloader.R
```

## Loading data

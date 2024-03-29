
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

## Running the complete pipeline

The following commands will download the raw data inputs from Figshare and run the complete analysis pipeline. However, certain workflows, such as the `ensemble_prediction.snake`, will be slow to execute without the use of highly parallel computation. 
```
$Rscript src/figshare_downloader.R
$snakemake --configfile snake_config.json --cores 1
```
See the following section for ways to experiment with specific workflows and/or download computationally expensive results from Figshare.

## Running individual workflows

Individual workflows group the scripts that are used for related analyses into snakefiles ('.snake') as defined above. Most workflows depend on data files generated by other workflows so some output data will probably need to be downloaded in order to run a workflow independently. We suggest downloading the entire collection of processed data from [Figshare](https://figshare.com/s/94f7ce1b3820c2c7d64d) and moving the files into the `data/processed` directory. This serves to decouple the workflows from each other, but the downloaded files will be automatically overwritten when executing a workflow that produces them. 

To run and individual workflow, such as `selecting_datasets.snake`:

```
$snakemake -s selecting_datasets.snake --configfile snake_config.json --cores 1
```

## License

BSD 3-Clause License

Copyright (c) 2023, John Michael Krill-Burger

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

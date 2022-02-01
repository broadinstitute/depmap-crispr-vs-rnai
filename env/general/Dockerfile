FROM bioconductor/bioconductor_docker:RELEASE_3_14

LABEL maintainer="mburger@broadinstitute.org"

COPY ./requirements.R /tmp/requirements.R 
RUN Rscript /tmp/requirements.R 

RUN apt-get install python3
RUN pip3 install snakemake

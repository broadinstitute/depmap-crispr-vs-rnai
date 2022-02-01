FROM rocker/r-ver:3.4.1

LABEL maintainer="mburger@broadinstitute.org"

RUN mkdir -p /usr/tda
WORKDIR /usr/tda

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y emacs
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y curl

ENTRYPOINT ["/usr/bin/env"]
CMD ["bash"]

RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3.6.5 \
    && conda update conda \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

RUN conda config --prepend channels conda-forge \ 
    && conda config --prepend channels bioconda \
    && conda config --append channels r

RUN conda install -y -c r r=3.4.1 \
    && conda install -y r-essentials \
    && conda install -y r-e1071 \
    && conda install -y r-matrixstats \
    && conda install -y r-sn \
    && conda install -y r-ggrepel \
    && conda install -c bioconda bioconductor-limma

RUN conda install -y r-diptest

RUN conda install -y scikit-learn=0.19.1 \
    && conda install -y numpy=1.13.1 \
    && conda install -y rpy2=2.9.4 \
    && conda install -y pandas=0.23.4

RUN conda install -y snakemake

RUN echo v2 && mkdir /tmp/sparklespray && \
    cd /tmp/sparklespray && \
    curl 'https://codeload.github.com/broadinstitute/sparklespray/zip/master' --output master.zip && \
    unzip master.zip && \
    cd sparklespray-master && \
    conda create -n sparkles -y python=3.5 && \
    /usr/local/envs/sparkles/bin/python -m pip install -r requirements.txt && \
    /usr/local/envs/sparkles/bin/python setup.py install

RUN apt-get update && apt-get install -y lsb-release && \
    apt-get install -my wget gnupg

RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y --allow-unauthenticated google-cloud-sdk

COPY packages /usr/tda/packages



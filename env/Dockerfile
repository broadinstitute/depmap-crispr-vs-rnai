FROM ubuntu:focal

LABEL maintainer="mburger@broadinstitute.org, aborah@broadinstitute.org"

RUN mkdir -p /usr/tda
WORKDIR /usr/tda

ENTRYPOINT ["/usr/bin/env"]
CMD ["bash"]

RUN apt-get update

RUN apt-get install -y libcurl4-openssl-dev \
      libssl-dev \
      bzip2 \
      curl \
      lsb-release \
      wget \
      gnupg \
      software-properties-common \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

#RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
#    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
#    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
#    apt-get update && apt-get install -y --allow-unauthenticated google-cloud-sdk \
#    && apt-get -qq -y autoremove \
#    && apt-get autoclean \
#    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

#RUN conda config --prepend channels conda-forge \
#    && conda config --prepend channels bioconda \
#    && conda config --append channels r

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
   add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
   apt-get update && \
   apt-get install -y r-base r-base-dev

#COPY env/install-R-deps.R /root/install-R-deps.R
#RUN mkdir -p /install && \
#    Rscript /root/install-R-deps.R


RUN R -e "install.packages('ggrepel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('feather',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('essentials',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringi',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('cli',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('e1071',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('matrixStats',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('diptest',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('sn',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('wesanderson',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggridges',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('limma', ask = F)"

RUN apt install -y python3-pip && \
    pip3 install --upgrade pip

RUN  apt-get install python3

RUN pip3 install scikit-learn
RUN pip3 install numpy==1.18.5
RUN pip3 install rpy2==2.9.4
RUN pip3 install pandas
RUN pip3 install snakemake
RUN apt-get install -y libncurses5-dev libncursesw5-dev
RUN pip3 install boto3

RUN apt-get update && apt-get install -y unzip gcc tzdata \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ENV TZ=America/New_York

#RUN conda create -n  sparkles python=3.6 -y

#RUN conda upgrade -n base -c defaults --override-channels conda

#RUN cd /tmp/ \
#&& wget https://github.com/broadinstitute/sparklespray/releases/download/v3.13.3/sparklespray-3.13.3.tar.gz \
#&& tar -xvzf sparklespray-3.13.3.tar.gz \
#&& cd sparklespray-3.13.3 \
#&& /usr/local/envs/analysis_env/bin/python setup.py install

RUN apt-get update && apt-get install -y libxml2-dev libxml2

RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

COPY packages/Python env/packages/Python

WORKDIR /tmp/pipeline

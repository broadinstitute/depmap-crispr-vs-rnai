FROM bioconductor/bioconductor_docker:RELEASE_3_14

LABEL maintainer="mburger@broadinstitute.org"

RUN mkdir -p /usr/tda
WORKDIR /usr/tda

ENTRYPOINT ["/usr/bin/env"]
CMD ["bash"]

RUN apt-get update

COPY env/requirements.R /tmp/requirements.R
RUN Rscript /tmp/requirements.R

RUN apt-get install python3

RUN pip3 install scikit-learn
RUN pip3 install numpy==1.18.5
RUN pip3 install rpy2==2.9.4
RUN pip3 install pandas
RUN pip3 install snakemake
RUN apt-get install -y libncurses5-dev libncursesw5-dev

RUN apt-get update && apt-get install -y unzip gcc tzdata \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ENV TZ=America/New_York

COPY src/ensemble_prediction_pipeline/env/packages/Python env/packages/Python

WORKDIR /tmp/pipeline

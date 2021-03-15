FROM bioconductor/bioconductor_docker:latest

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    libglpk-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e "devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos=BiocManager::repositories())"

# Getting base image Tidyverse
FROM rocker/tidyverse:4.0.3

MAINTAINER Glenn <glchang@bcgsc.ca>

# Instal R CRAN package
RUN apt-get update && \
        apt-get install -y build-essential libglpk40 && \
        install2.r --error --skipinstalled \
        optparse && \
        rm -rf /tmp/downloaded_packages

# Install bioconductor packages
RUN install2.r --error --skipinstalled BiocManager && \
        R -e 'BiocManager::install(ask = F)' && \
        R -e 'BiocManager::install(c("polyester","Biostrings", "tximport", "biomaRt", ask = F))'

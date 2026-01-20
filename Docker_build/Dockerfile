## This is a two stage build to reduce image size ##

## -- Stage 1 -- ##
FROM rocker/r-ver:4.5.2 AS builder

RUN apt-get update && apt-get install -y \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libproj-dev \
    proj-bin \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install R dependencies
# -- CRAN packages --
RUN R -e "install.packages(c('stringr', 'dplyr', 'ggplot2', 'tidyr', 'readr', 'optparse', 'this.path', 'remotes', 'outliers', 'ggpubr', 'ggsci', 'cli'))"

# -- GitHub packages --
RUN R -e "remotes::install_github('hrbrmstr/ggalt')"



## -- Stage 2 -- ##
# Base image for R
FROM rocker/r-ver:4.5.2 AS runtime

# Container directory
RUN mkdir /home/XenoR
WORKDIR /home/XenoR

# Script I/O directories
RUN mkdir /home/XenoR/Input_files && \
    mkdir /home/XenoR/Output_files

# Install system dependecies
RUN apt-get update && apt-get install -y \
    libcurl4 \
    libssl3 \
    libxml2 \
    libproj25 \
    proj-bin \
    && rm -rf /var/lib/apt/lists/*

#  Copy the installed libraries and binaries from the builder image
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

# Move the XenoScript tools
COPY SynthXenoGen_v1.0.1.r .
COPY XenoVol_v1.0.1.r .

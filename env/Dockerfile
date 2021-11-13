 FROM ubuntu:20.04
 COPY environment.yml .
 ARG DEBIAN_FRONTEND=noninteractive
 
 # Updating repository and getting some standard packages
 RUN apt-get -y update && apt-get install -y \
        git \
        software-properties-common \
        tar \
        wget
        
# Installing R 4.1.1
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt install -y r-base r-base-core r-recommended r-base-dev

# libtbb2 is a requirement for sina
RUN add-apt-repository universe && \
    apt-get install libtbb2

# Installing the required packages for tidyverse
RUN apt install -y libssl-dev libxml2-dev libcurl4-openssl-dev
# Installing the R packages: tidyverse and ape
RUN R -e 'install.packages(c("tidyverse", "ape"))'

# Installing conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -bfp /usr/local/miniforge && \
    rm Miniforge3-Linux-x86_64.sh 

# Using conda to install snakemake, vsearch, sina, and fasttree
RUN /usr/local/miniforge/bin/conda env create -f environment.yml

# Appending PATH so all installed programs are easily available
ENV PATH=$PATH:/usr/local/miniforge/envs/snakemake/bin:/usr/local/miniforge/bin

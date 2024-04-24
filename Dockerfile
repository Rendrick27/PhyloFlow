# Use snakemake as base image
FROM snakemake/snakemake:latest

# Install requirements via conda
RUN mamba install mafft=7.520 -c bioconda
RUN mamba install conda-forge::biopython
RUN mamba install -c bioconda -c conda-forge modeltest-ng=0.1.7
RUN mamba install -c bioconda -c conda-forge raxml-ng=1.2.0

# Install python packages to create the final phylo trees
RUN pip install --force-reinstall -v "toyplot==1.0.3"
RUN pip install --force-reinstall -v "toytree==2.0.1"

WORKDIR /App
COPY . /App


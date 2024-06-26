FROM snakemake/snakemake:v8.11.4

# Install requirements via conda
RUN mamba install mafft=7.520 -c bioconda
RUN mamba install conda-forge::biopython
RUN mamba install -c bioconda -c conda-forge modeltest-ng=0.1.7
RUN mamba install -c bioconda -c conda-forge raxml-ng=1.2.0

# Install python packages
RUN pip install --force-reinstall -v "toyplot==1.0.3"
RUN pip install --force-reinstall -v "toytree==2.0.1"

WORKDIR /phylo_flow
COPY . /phylo_flow
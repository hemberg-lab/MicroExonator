# Introduction

MicroExonator is a fully-integrated computational pipeline that allows for systematic de novo discovery and quantification of microexons using raw RNA-seq data for any organism with a gene annotation. Compared to other available methods MicroExonator is more sensitive for discovering smaller microexons and it provides higher specificity for all lengths. Moreover, MicroExonator provides integrated downstream comparative analysis between cell types or tissues using [Whippet](https://github.com/timbitz/Whippet.jl) ([Sterne-Weiler et al. 2018](https://doi.org/10.1016/j.molcel.2018.08.018)).


# Installation

Start by cloning MicroExonator

    git clone https://github.com/hemberg-lab/MicroExonator

Install [Miniconda 3](https://docs.conda.io/en/latest/miniconda.html)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh 
    ./Miniconda3-latest-Linux-x86_64.sh


Finally, create an enviroment to run [snakemake](https://snakemake.readthedocs.io/en/stable/)

    conda create -n snakemake_env -c bioconda -c conda-forge snakemake
    
    
# Documentation 

Extended documentation can be found at https://microexonator.readthedocs.io.

The default confidence filter is now the quantitative robustness route. Bulk
RNA-seq runs must declare biological sample groups in a tab-separated
`bulk_samples.tsv` file and set `bulk_samples` in `config.yaml`. Pure
single-cell runs may omit this file: cluster metadata defines the groups when
provided, otherwise all cells are treated as one `all_cells` group. See the
Discovery and Quantification documentation for the full configuration and the
explicit `legacy_mixture` compatibility option.


# Contact

For questions, ideas, feature requests and potential bug reports please contact gp7@sanger.ac.uk.

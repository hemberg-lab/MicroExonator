# MicroExonator

MicroExonator is a Snakemake workflow for the de novo discovery and quantification of microexons (exons of 30 nt or less) from raw RNA-seq data, for any organism with a genome and a gene annotation. Compared with other methods, it is more sensitive for the smallest microexons and more specific at all lengths. It reports PSI values per sample and can test differential inclusion between groups of samples, either with its own statistical test or with [Whippet](https://github.com/timbitz/Whippet.jl) ([Sterne-Weiler et al. 2018](https://doi.org/10.1016/j.molcel.2018.08.018)). Single-cell RNA-seq data are supported through pseudo-bulk comparisons between cell types.

## MicroExonator 2.0.0 (beta)

This branch holds MicroExonator 2.0.0, now in beta. It runs under Python 3 and fixes several quantification problems that affected a sizeable share of microexons: read length handling, microexons without an annotated skipping junction, consecutive microexons, microexons that share a splice site with a longer exon, and junctions with a copy in the genome. It also adds differential inclusion without Whippet and new diagnostic reports. See [What's new](https://microexonator.readthedocs.io/en/latest/whats_new.html) for details.

Coming from an earlier version: the default confidence filter is now the robustness filter, so bulk RNA-seq runs need a `bulk_samples` file that assigns each sample to a condition. The original Gaussian mixture filter is still available as `filter_method : legacy_mixture`, and needs no `bulk_samples` file.

The *Methods in Molecular Biology* protocol corresponds to the `MiMB` branch, which is kept unchanged for reproducibility.

## Installation

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (or any conda distribution), then create an environment with Snakemake 7.32.4, the version MicroExonator is tested with:

    conda install -n base -c conda-forge mamba
    mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=7.32.4

Clone MicroExonator and check out the 2.0.0 beta:

    git clone https://github.com/hemberg-lab/MicroExonator
    cd MicroExonator
    git checkout refurbishment/legacy-fixes

Every other dependency is installed by Snakemake through conda the first time each step runs.

## Documentation

The full documentation, including input files, configuration, cluster setup and output formats, is at https://microexonator.readthedocs.io.

## Citation

Parada GE, Munita R, Georgakopoulos-Soares I, Fernandes HJR, Kedlian VR, Metzakopian E, Andres ME, Miska EA, Hemberg M (2021). MicroExonator enables systematic discovery and quantification of microexons across mouse embryonic development. *Genome Biology* 22:43. https://doi.org/10.1186/s13059-020-02246-2

## Contact

For questions, ideas, feature requests and bug reports, please open an [issue](https://github.com/hemberg-lab/MicroExonator/issues) or write to guillermo_parada@kcl.ac.uk.

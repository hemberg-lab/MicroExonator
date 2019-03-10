# MicroExonator
reproducible discovery and quantification of microexons using RNA-seq data


# Install

clone MicroExonator

    git clone https://github.com/geparada/MicroExonator

Then, install Miniconda 2.

    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    chmod +x Miniconda2-latest-Linux-x86_64.sh
    ./Miniconda2-latest-Linux-x86_64.sh -b -p /cvmfs/softdrive.nl/$USER/Miniconda2
    
Finnaly create an enviroment to run snakemake
    
    conda create -n snakemake snakemake python=3.6 pandas cookiecutter 

# Configure

Before runnig MicroExonator you need to have certain files at `MicroExonator\`. You need to create `config.yaml` and `cluster.json` inside . Finnaly, to declare the input, you need to create `desing.tsv` (for fastq.gz files that are locally in your machine) and/or `NCBI_accession_list.txt`(for SRA accession names). 

The input need to be declared at

# Run

If you are working remotelly or even in your own machine, we higly recomment create an screen before runing MicroExonator

    screen -S test_run

Then activate snakemake enviroment

    source activate snakemake

Then run

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster [cluster system params] --use-conda -k  -j 40` .
    
Where in the case of lsf, cluster system params will be "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"
    
For large dataset runnig the pipeline in two strands:
    
    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster [cluster system params] --use-conda -k  -j 40 discovery
    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster [cluster system params] --use-conda -k  -j 40 quant
  


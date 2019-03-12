# Introduction

MicroExonator is a fully-integrated computational pipeline that allows for systematic de novo discovery and quantification of microexons using raw RNA-seq data for any organism. Compared to other available methods MicroExonator is more sensitive to the discovery of smaller microexons. Moreover, MicroExonator provides integrated downstream comparative analysis between cell types or tissues using Whippet (Sterne-Weiler et al. 2018). As a proof of principle MicroExonator identified X novel microexons in Y RNA-seq samples from mouse and systematically characterised microexons in terms of tissue and cell type specificity.


# Install

clone MicroExonator.

    git clone https://github.com/geparada/MicroExonator

Then, install [Miniconda 2](https://docs.conda.io/en/latest/miniconda.html)

    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    chmod +x Miniconda2-latest-Linux-x86_64.sh
    ./Miniconda2-latest-Linux-x86_64.sh -b -p /cvmfs/softdrive.nl/$USER/Miniconda2

Finally create an enviroment to run snakemake

    conda create -n snakemake snakemake python=3.6 pandas cookiecutter

# Configure

Before running MicroExonator you need to have certain files at `MicroExonator/`. You need to create `config.yaml` and `cluster.json` inside . Finally, to declare the input, you need to create `desing.tsv` (for fastq.gz files that are locally in your machine) and/or `NCBI_accession_list.txt`(for SRA accession names). Take a look to Examples folder.

The `config.yaml` file is an standar yaml file that should have the path of the input files and certain paramethers:

    Genome_fasta : /path/to/danRer11.fa
    Gene_anontation_bed12 : /path/to/danRer11.ensembl.bed12
    GT_AG_U2_5 : /path/to/danRer11_GT_AG_U2_5.good.matrix
    GT_AG_U2_3 : /path/to/Zebrafish/Data/danRer11_GT_AG_U2_3.good.matrix
    vertebrates_phylop : /path/to/danRer11.chrom.sizes.bw.sep  
    working_directory : /path/to/Zebrafish/
    ME_DB : /path/to/VastDb.bed12
    ME_len : 30

Whereas `Genome_fasta` is a multifaste file containg the genome cromosomes. `Gene_anontation_bed12` is a BED file containing the transcript annotation, which can be found at [UCSC genome browser](http://genome.ucsc.edu/cgi-bin/hgTables). `GT_AG_U2_5` and `GT_AG_U2_5` are splice site PWMs that can be from [SpliceRack](http://katahdin.cshl.edu/SpliceRack/poster_data.html) (as this server is currently down, we provisionally provide the PWM for human and mouse), but if you do not have these PWM for the specie for which you want to conduct the analysis, you can leave it as `NA` and MicroExonator will generate these PWMs internaly.   

# Run

If you are working remotely or even in your own machine, we highly recommend create an screen before running MicroExonator

    screen -S session_name

Then activate snakemake enviroment

    source activate snakemake

Then run

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs}

Notice that you should use `--cluster` only if you are working in a computer cluster that works with queue systems such as lsf, qsub, SLURM, etc. We provide an example of `cluster.json` to work with lsf and in that case the cluster system params should be `"bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"`. The number of parallel jobs can have any int number, this depend of your machine work load capacity. 

If you want to process a large dataset, we recommend to run MicroExonator in two stages:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} discovery
    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} quant
    
By doing this, you will optimise disk space, which is often a restrictive resource for running large data sets. 

# Troubleshooting

Before running it is recommended to see if SnakeMake can corretly generate all the steps given your input. For this, you can do a dry-run using `-np` parameters:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} -np

If the dry-run cannot be iniciated, make sure you are running MicroExonator from inside the folder you cloned from this repository. Make sure you have the right configuration inside `config.yaml`. 

The current version of MicroExonator does not support chromosome names that has `_` or `|`, for example some chromosome names can be `chr1_KZ111v2_alt` and in this case you migth have some errors that will prevent you to complete the run. For now we recommed to replace these caracters by any string you can recognise, for example:

    sed 's/_/SEP/g' genome.fa > genome.fa.sed
    
And then use this modified genome (`genome.fa.sed`) as an input. Future versions of MicroExonator will overcome this issue withouth requering this step.


# Contact

For questions, ideas, feature requests and potential bug reports please contact gp7@sanger.ac.uk.

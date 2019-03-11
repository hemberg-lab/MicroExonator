environment# MicroExonator
reproducible discovery and quantification of microexons using RNA-seq data


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

# Run

If you are working remotely or even in your own machine, we highly recommend create an screen before running MicroExonator

    screen -S session_name

Then activate snakemake enviroment

    source activate snakemake

Then run

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs}

Notice that you should use `--cluster` only if you are working in a computer cluster that works with queue systems such as lsf, qsub, SLURM, etc. We provide an example of `cluster.json` to work with lsf and in that case the cluster system params should be `"bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"`. The number of parallel jobs can have any int number, this depend of your machine work load capacity. 

For large dataset running the pipeline in two strands:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} discovery
    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} quant

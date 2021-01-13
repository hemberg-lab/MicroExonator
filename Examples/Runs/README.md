Here we provide example runs that we have implemented for two different projects:

# Zebrafish

Small project that were ran using SRA accession codes as input. These accession codes are inputed inside `NCBI_accession_list.txt` file.

# COSMIC

Large cancer cell-lines project, where we used a local copy of the input fastq.gz files as an input. The paths and the name of the samples needs to be provided inside a `desing.tvs` file.


# Running under lsf

`snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" --use-conda -k  -j 1000000`

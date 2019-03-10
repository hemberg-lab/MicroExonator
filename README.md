# MicroExonator
reproducible discovery and quantification of microexons using RNA-seq data

# Instructions

* git clone https://github.com/geparada/MicroExonator
* create `config.yaml` and `cluster.json` inside `MicroExonator\`
* run `snakemake  --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" -k -s MicroExonator.skm --use-conda  -j 40`


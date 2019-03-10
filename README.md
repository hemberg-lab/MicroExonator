# MicroExonator
reproducible discovery and quantification of microexons using RNA-seq data

# Instructions

* git clone https://github.com/geparada/MicroExonator
* create `config.yaml` and `cluster.json` inside `MicroExonator\`
* run `snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster [cluster system params] --use-conda -k  -j 40` . Where in the case of lsf, cluster system params will be `"bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"`


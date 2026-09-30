# Example runs

These are the configurations of runs from the original publication and related projects. Paths point to the cluster where they were run and must be replaced with your own, and some files use settings from earlier versions (for example, `whippet_delta` given inside `config.yaml`; it is now the path of a separate YAML file, as in `Parada_et_al/whippet_delta.yaml`). Use them as a guide to the input files; the documentation at https://microexonator.readthedocs.io describes the current settings.

- `Parada_et_al`: the mouse embryonic development runs of the publication, with SRA accessions, URL inputs and the `whippet_delta.yaml` comparisons.
- `Autism`, `Zebrafish`: runs from SRA accession codes, listed in `NCBI_accession_list.txt`.
- `ENCODE`: FASTQ files downloaded from URLs listed in `sample_url.tsv`.
- `COSMIC`: a large set of cancer cell lines from local FASTQ files, listed with their sample names in `local_samples.tsv`.
- `C_elegans`: a configuration for *C. elegans*.

Cluster configuration examples for LSF and PBS are in `Examples/Cluster_config/`; the command to submit jobs is described in the documentation.

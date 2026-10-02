# Example runs

`Template/` is the starting point for a new MicroExonator 2.0.0 run: a `config.yaml` with the default robustness filter and the differential inclusion test that needs no Whippet, plus example `local_samples.tsv`, `bulk_samples.tsv` and `whippet_delta.yaml` files. Replace the `/path/to/` paths with your own.

The other folders are the configurations of runs from the original publication and related projects, updated to the current settings. Their paths are placeholders, they use `filter_method : legacy_mixture` as the original runs did, the ones with Whippet set `delta_method : whippet`, and comparisons are in a separate `whippet_delta.yaml`. Bulk runs also need a `bulk_samples` file that assigns each sample to a condition; these examples do not include one yet, except that `Autism` runs without one because it uses `downstream_only`.

- `Parada_et_al`: the mouse embryonic development runs of the publication, with SRA accessions, URL inputs and the `whippet_delta.yaml` comparisons.
- `Autism`, `Zebrafish`: runs from SRA accession codes, listed in `NCBI_accession_list.txt`.
- `ENCODE`: FASTQ files downloaded from URLs listed in `sample_url.tsv`.
- `COSMIC`: a large set of cancer cell lines from local FASTQ files, listed with their sample names in `local_samples.tsv`.
- `C_elegans`: a configuration for *C. elegans*.

Cluster configuration examples for LSF and PBS are in `Examples/Cluster_config/`. The documentation at https://microexonator.readthedocs.io describes every setting.

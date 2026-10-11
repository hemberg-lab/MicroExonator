# Cached gene-expression QC

Optional standalone Snakemake 7 target, using retained comparison count matrices.
It never acquires FASTQs or produces missing counts, alignments or processing QC.
Existing `quant_umbrella` execution and expression results are unchanged.

Use a separate `expression-qc.yaml`:

```yaml
umbrella_expression_qc_sources: [featurecounts, salmon]
umbrella_expression_qc_root: comparisons
umbrella_expression_qc_cooks: true
# Optional restriction; omit to process all cached comparisons:
# umbrella_expression_qc_comparisons:
#   - reference-id/project-id/comparison-id
```

Source settings:

| Setting | Reports |
| --- | --- |
| `[]` (default) or `none` | Disabled; no matrix files required |
| `[featurecounts]` | featureCounts only |
| `[salmon]` | Salmon only |
| `[featurecounts, salmon]` | Both, separate reports |

Comma-separated CLI syntax is also accepted:

```bash
snakemake -s GeneExpressionQC.smk gene_expression_qc --use-conda --cores 1 \
  --config umbrella_expression_qc_sources=featurecounts,salmon -n
```

Or use YAML, first dry-run and then execute:

```bash
snakemake -s GeneExpressionQC.smk gene_expression_qc --configfile expression-qc.yaml --use-conda --cores 1 -n
snakemake -s GeneExpressionQC.smk gene_expression_qc --configfile expression-qc.yaml --use-conda --cores 1
```

Normal conda frontend/prefix options may be added. `salmon` maps internally to
`tximport`, matching existing joined count/TPM/effective-length matrices and
output names. Root must point to the actual comparison directory (for example
`umbrella/comparisons` for a downloaded pilot). Every selected comparison must
have the requested source caches. Missing required matrices fail at DAG creation;
missing optional processing QC is reported, never regenerated.

## One report per comparison and source

`gene_expression_qc_<route>.pdf` contains seven pages:

1. Original PC1-PC2 plot, with actual comparison-group colours and ggrepel labels.
2. Additional PCs and gene-selection sensitivity: up to 500, 2,000, all variable
   filtered genes. Small panels deduplicate equivalent selections.
3. Clustered RMS expression-distance heatmap, all filtered VST genes.
4. Pearson correlation heatmap, up to 500 most variable VST genes.
5. Gene-assigned counts, detected genes, size factors, normalized distributions.
6. Retained mapping/assignment/strandedness metrics, if available.
7. Cook's gene-level influence diagnostics, or an explicit unavailable/disabled page.

Two compact compressed tables retain PCA coordinates and per-sample summary
metrics/availability. Technical runs are collapsed into biological replicates
using the existing comparison definition. Processing-QC ratios are simple means
across available processing runs of a replicate, not pooled read-weighted ratios;
partial coverage is explicitly marked. Inferred strandedness categories are
retained in the summary table. Gene-body coverage is not included in this version.

Preprocessing is the same as the original PCA: at least 10 counts in two samples,
DESeq2 size factors (`poscounts`) and blind VST with mean dispersion fit. Gene
selection never uses differential significance. Salmon retains length-aware
normalization using the existing tximport convention. At least three samples and
two variable genes are required for the descriptive report.
Salmon gene/sample offsets do not provide scalar size factors; those values are
marked not applicable, rather than substituting a misleading scalar estimate.

## Interpretation and resource use

There are **no automatic sample-outlier shapes or exclusions** in this report.
The old PCA-only workflow remains available for compatibility, including its
explicitly exploratory flags; this report does not use those flags to classify
samples. Separated samples may reflect library protocol, tissue composition,
batch, disease heterogeneity or technical failure. Expression-derived plots are
not independent confirmations of one another.

Cook's diagnostics fit `~ group` without count replacement; they do not change
stored DE results. They require at least two replicates per group and a supported
comparison preflight. Fit failures are recorded in the summary and PDF. The
99th-percentile F reference cutoff applies to gene-level Cook's distances, not
sample-level significance. Set `umbrella_expression_qc_cooks: false` to skip this
extra model fit while retaining every descriptive plot.

The rule requests one CPU and 4 GB RAM. ggplot2, ggrepel, gridExtra, DESeq2 and
jsonlite are conda-managed in the plotting environment. No separate installation
is required. Heatmap/sample-distance memory grows quadratically with sample
number; unusually large comparisons may need a resource override.

See [DESeq2 QC and influence guidance](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
and [ggrepel](https://ggrepel.slowkow.com/articles/examples.html).

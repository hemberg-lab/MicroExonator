# Cached gene-expression PCA

Run after expression comparison joins have completed. This optional standalone
workflow contains no count-matrix producers, alignment rules or read acquisition
rules. Missing cached inputs fail before execution; existing analyses are not
modified. It does not change `quant_umbrella` or its conda environments.

```bash
snakemake -s ExpressionPCA.smk expression_pca --use-conda --cores 1 -n
snakemake -s ExpressionPCA.smk expression_pca --use-conda --cores 1
```

`quant_umbrella` also makes these plots itself, in the same call that quantifies
the runs: every comparison gets `pca_tximport.*` when Salmon is on and
`pca_featurecounts.*` when HISAT2 is on, with the same script and conda env. Set
`umbrella_pca_routes` in the main config to choose routes, or `[]` for none:

```yaml
umbrella_pca_routes: [tximport]   # tximport needs salmon, featurecounts needs hisat2
```

The default uses featureCounts for every cached comparison under `comparisons/`.
Use a separate optional YAML to select comparisons/routes:

```yaml
umbrella_pca_root: comparisons
umbrella_pca_routes: [featurecounts, tximport]
umbrella_pca_comparisons:
  - hg38/project-id/comparison-id
umbrella_pca_min_group: 5
```

Pass it with `--configfile pca.yaml`. Selection keys are
`reference_id/project_id/comparison_id`, not just comparison names. Omit
`umbrella_pca_comparisons` to select all cached preflights. Requested routes must
have all their cached joined matrices; absent expression modules cannot provide
PCA. At least three distinct biological replicates and two variable genes are
required. The plot does not automatically remove batch effects.

Each comparison/route emits only `pca_<route>.pdf` and
`pca_<route>.scores.tsv.gz`. Points are coloured by actual preflight comparison
labels, labelled with biological replicate IDs, and positioned using
`ggrepel::geom_text_repel(seed=42, max.overlaps=Inf)`.
The blue/orange palette distinguishes the first/second comparison side; colours
are not global condition identities across unrelated comparisons. Dense figures
can still need larger page dimensions despite label repulsion.

Counts are filtered without differential results (at least 10 counts in at least
two replicates), size-factor normalized, and variance stabilized with DESeq2
(`blind=TRUE`, mean dispersion fit). Salmon uses the existing length-aware
tximport convention, including technical-run effective-length averages. PCA uses
up to 500 highest-variance genes, without unit-variance scaling. Each point is a
biological replicate, with technical runs already collapsed by the comparison
join. Raw matrices remain unchanged.

## Exploratory within-group flags

Using all retained transformed genes, calculate each sample's median RMS
expression distance to the other members of its own group. Flag `Review` when
that score exceeds the group's median score plus three scaled MADs
(`stats::mad`, default normal-consistency scaling). These are transparent
heuristic flags, **not calibrated p-values, proof of technical failure, or
automatic exclusions**. Multiple outliers can mask one another.

With fewer than five biological replicates in a group, or zero/near-zero score
dispersion, assessment is disabled. The minimum of five is a pragmatic guard,
not a validated statistical guarantee; it may be increased but not lowered.
Circles mean unflagged, triangles mean review, open circles mean not assessed.
The compressed table records scores, thresholds, reasons and PC variance.

The PCA conda environment contains DESeq2, ggplot2, ggrepel and jsonlite. It is
separate to avoid invalidating existing expression jobs when installing plotting
dependencies. Nothing needs to be installed outside conda.

References: [ggrepel](https://ggrepel.slowkow.com/articles/examples.html),
[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

# MegaSearch umbrella mode

Umbrella mode runs several quantification and differential tools side by side
on the same samples, keeps only compact group-by-batch tables, and supports
adding batches of samples over time. It is opt-in: without `umbrella_manifest`
in the config, MicroExonator behaves exactly as before.

## Inputs

**Manifest** (`umbrella_manifest`, TSV). Required columns:

| Column | Values |
|---|---|
| `sample_id`, `run_id`, `biological_replicate_id` | IDs (letters, digits, `_`, `.`, `-`) |
| `project_id`, `batch_id`, `group` | IDs |
| `source_type` | `sra`, `fastq`, `cram`, `bam` |
| `source_1`, `source_2` | accession or path; `source_2` only for local PE FASTQ |
| `layout` | `SE`, `PE` |
| `strandedness` | `unstranded`, `firststrand`, `secondstrand`, `auto` |
| `reference_id` | the content-derived ID of the reference bundle |
| `include` | `true`, `false` (excluded runs need a `reason` column) |

Runs of one `biological_replicate_id` are technical runs: count tables are
summed per replicate and they never count as independent replicates.

**Reference bundle** (`umbrella_reference`, mapping). `genome_fasta`,
`annotation_gtf` (full annotation: featureCounts, rMATS, splice hints),
`whippet_gtf` (the fixed Whippet annotation), `me_db`, `transcriptome_fasta`, `decoys`, `salmon_gtf`,
`splice_sites` (HISAT2 format, full annotation plus any extra junctions).
Optional prebuilt indexes: `hisat2_index_prefix` (with `hisat2_index_type:
large` for `.ht2l`), `whippet_index`, `salmon_index`. The workflow hashes the
whole bundle into `umbrella/reference/<id>/manifest.json` and refuses a
manifest whose `reference_id` does not match.

**ME_DB microexons in the GTFs.** The configured GTFs are never changed.
The workflow writes derived copies under `umbrella/reference/<id>/`
(`src/insert_microexons_gtf.py`): each ME_DB microexon that is not already an
exon is added to one copy of one host transcript whose intron contains it
(ranked MANE_Select, Ensembl_canonical, basic, then more exons). One host is
enough for Whippet, which joins any annotated donor to any annotated acceptor
of a gene. A report lists each microexon as present, inserted or without host.
Whippet, the junction labels, the HISAT2 splice hints, featureCounts and rMATS
use the derived copies; Salmon and SUPPA2 keep `salmon_gtf`, whose transcripts
must match `transcriptome_fasta`. A prebuilt `whippet_index` must have been
built from `whippet.microexons.gtf.gz`. Set `umbrella_insert_microexons:
false` to use the configured GTFs as they are.

**Comparisons** (`umbrella_comparisons`, YAML): a list under `comparisons`
with `comparison_id`, `project_id`, `group_a`, `group_b` and optional
`exclude`. Effects are reported as A - B.

**Optional analyses** (`umbrella_optional`, all `true` by default): `rmats`,
`coverage`, `leafcutter` (also needs `leafcutter_dir`, a LeafCutter checkout,
and optionally `leafcutter_rscript`), `suppa2`.

## Targets

| Target | Reaches |
|---|---|
| `quant_microexonator` | the legacy `quant` outputs and MicroExonator group shards |
| `quant_whippet` | the above, plus native Whippet quantification and Whippet shards |
| `quant_umbrella` | everything: alignment shards, Salmon shards, joins, comparisons, synthesis |

`quant` and `get_whippet_psi` keep their legacy meaning.

## What is kept

Everything per run is temporary, including the staged reads, except the
rMATS prep files that `--task post` needs later. Kept per
`{reference_id}/{project_id}/{group}/{batch_id}`:

| Directory | Files |
|---|---|
| `junctions/` | junction counts with index labels; per-site index capture; capture rate per run |
| `genes/` | featureCounts; Salmon gene counts, TPM, effective length; Salmon transcript counts and TPM |
| `qc/` | one row of raw QC per run (mapping, strandedness, spliced reads, featureCounts assignment) |
| `splicing/` | MicroExonator corrected PSI and Whippet PSI, every column the delta tools read |
| `rmats/` | one prep file per run and an inventory |
| `coverage/` | summed CPM coverage (bigWig) and the number of runs, so group means stay exact across batches |

Every kept file is written once, read-only, next to a checksum guard. Writing
different content to an existing shard fails instead of overwriting it. Per
comparison, `comparisons/{reference_id}/{project_id}/{comparison_id}/` holds
the preflight, joined matrices, each tool's output and `synthesis.tsv`.

## Comparisons

The preflight (`preflight.json`) decides whether inference is supported:
at least two biological replicates per group, and group not confounded with
batch, layout or strandedness. When it is not, every tool writes one comment
line with the reason and runs no test; joined matrices are still made. rMATS
also needs one layout across the comparison.

Tools: DESeq2 `~ group` on featureCounts and on Salmon
(`DESeqDataSetFromTximport`), MicroExonator delta, Whippet delta, rMATS-turbo
(`inte` then `post`), LeafCutter and SUPPA2 (pilots). `synthesis.tsv` lists,
per tool, its status, calls, direction, calling rule and shared upstream data,
plus per-run index capture and QC outliers computed across each group. It
does not vote.

## Operation

- **Batches.** Add runs with a new `batch_id` to the manifest and rerun
  `quant_umbrella`. Existing shards are checked, not rewritten; joins and
  comparisons rerun from shards.
- **Checking temporaries.** `snakemake --delete-temp-output --dry-run`
  lists what would be removed.
- **Local smoke test.** `tests/integration/mm10_simulated/run_umbrella_smoke.py`
  stages the mm10 fixture (two groups, two batches, SE and PE) against
  prebuilt full-mm10 indexes; it never builds an index. `--dry-run` needs no
  tools. The full run uses conda; rMATS and SUPPA2 packages are Linux-first.
- **Before a cluster pilot.** Record the real manifest and sample count, build
  the fixed Whippet GTF with every ME_DB microexon in it, compute the
  reference ID, and pin tool versions in `umbrella_reference: versions`.

## Known limits

- Whippet psi path columns (`Inc_Paths`, `Exc_Paths`, `Edges`) are not kept.
- Technical runs are summed in count tables, but the delta tools still see
  each run's PSI table separately (flagged as `technical_run_replicates`).
- Junction capture uses splice-site sets per chromosome and strand, not per
  gene.
- The planned Salmon versus Whippet TPM concordance is not implemented.

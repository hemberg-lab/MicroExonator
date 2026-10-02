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

**Libraries** (optional column `library_id`). Runs sharing a `library_id`,
such as the lanes or resequencing runs of one library, are merged into one
unit before any tool runs:
- **Staging** appends them in run order, both mates of a run together. The
  unit is named by the `library_id` everywhere: work folder, shard columns,
  MicroExonator sample, comparisons.
- **Leave it empty** for a run that is its own library; manifests without the
  column behave as before.
- **Members must agree** on project, batch, group, layout, strandedness,
  reference, source type and `biological_replicate_id`. Excluded runs are left
  out of their library. A `library_id` may not equal another run's `run_id`.
  Several BAM/CRAM runs cannot be merged. Runs whose first reads are
  identical are refused as a duplicate submission.
- **Members may differ** in `sample_id` (lanes can have their own GEO
  samples) and other metadata: those values are kept joined with `;`, and
  `samples:`/`runs:` in a comparison match a library through any of its runs.
  Preflights list each library's runs under `library_members`.
- **Decide the grouping upstream.** The workflow merges exactly what the
  column says. Group runs of one SRA experiment, or curated lane pairs; a
  shared BioSample or donor alone is not evidence of lanes. Set `library_id`
  before a run is processed: a processed run cannot be merged into a library
  afterwards without downloading it again, and the shard guard refuses it.

**Reference bundle** (`umbrella_reference`, mapping). Supply `genome_fasta`,
`annotation_gtf` (full annotation: featureCounts, rMATS, splice hints),
`whippet_gtf` (the fixed Whippet annotation) and `me_db`. The workflow derives
the transcript FASTA with gffread, the genome-contig decoy list and HISAT2
splice hints from those inputs. `salmon_gtf` defaults to `annotation_gtf` and
may be supplied separately. `transcriptome_fasta`, `decoys` and `splice_sites`
may still be supplied as overrides for an existing reference bundle.
Optional prebuilt indexes: `hisat2_index_prefix` (with `hisat2_index_type:
large` for `.ht2l`), `whippet_index`, `salmon_index`. The `reference_id` is a
hash of the configured files (prebuilt indexes included), `versions`,
`hisat2_index_type` and `umbrella_insert_microexons`; files the workflow
builds are recorded in `umbrella/reference/<id>/manifest.json` but are not
part of the ID. Write `auto` in the manifest's `reference_id` column and the
workflow computes it at start-up and prints it (`umbrella reference_id
(auto): ...`). File digests are cached in
`umbrella/reference/checksum_cache.json` by path, size and modification time,
so only the first start hashes the genome (about 3 s per GB); later parses,
including every cluster job, take milliseconds. An explicit ID also works
(`python3 src/shard_guard.py reference-id --configfile config.yaml` prints it);
`umbrella_reference_identity` checks it before any index is built.

**ME_DB microexons in the GTFs.** The configured GTFs are never changed.
The workflow writes derived copies under `umbrella/reference/<id>/`
(`src/insert_microexons_gtf.py`): each ME_DB microexon that is not already an
exon is added to one copy of one host transcript whose intron contains it
(ranked MANE_Select, Ensembl_canonical, basic, then more exons). One host is
enough for Whippet, which joins any annotated donor to any annotated acceptor
of a gene. A report lists each microexon as present, inserted or without host.
Whippet, the junction labels, the HISAT2 splice hints, featureCounts and rMATS
use the derived copies; Salmon and SUPPA2 keep `salmon_gtf`, whose transcripts
must match the derived transcript FASTA. A prebuilt `whippet_index` must have been
built from `whippet.microexons.gtf.gz`. Set `umbrella_insert_microexons:
false` to use the configured GTFs as they are.

**Comparisons** (`umbrella_comparisons`, YAML): a list under `comparisons`
with `comparison_id`, `project_id`, two sides and optional `exclude` (run or
biological replicate IDs). Effects are reported as A - B. A side is either a
manifest group (`group_a: case`) or a selection of runs with a label:

```yaml
comparisons:
  - comparison_id: asd_vs_ctrl
    project_id: PRJNA1120182
    group_a: cerebellum_asd
    group_b: cerebellum_unaffected
  - comparison_id: asd_vs_ctrl_female
    project_id: PRJNA1120182
    a: {label: asd_female, where: {phenotype: ASD, sex: F}}
    b: {label: ctrl_female, groups: [cerebellum_unaffected], where: {sex: F}}
```

A selection keeps the runs matching every key it gives: `groups`, `samples`
(sample IDs), `runs` (run IDs) and `where` (any manifest column, including
extra metadata columns; a value or a list of values).

**Tool selection** (`umbrella_optional`). Every tool can be switched off:

| Key | Default | Controls |
|---|---|---|
| `microexonator` | on | MicroExonator per run, its shards and `microexonator_delta` |
| `whippet` | on | Whippet per run, its shards and `whippet_delta` |
| `salmon` | on | Salmon per run, its shards and DESeq2 on tximport |
| `hisat2` | on | the analyses of the shared HISAT2 alignment: junctions and capture rates, featureCounts and its DESeq2, the compact QC table |
| `rmats`, `coverage`, `leafcutter` | on (need `hisat2`) | rMATS, summed coverage, LeafCutter |
| `suppa2` | on (needs `salmon`) | SUPPA2 |
| `multiqc` | on | FastQC plus the kept summaries of the selected tools, one MultiQC report per project |
| `microexonator_whippet_delta` | off (needs `microexonator` and `whippet`) | the legacy ME-in-Whippet delta, run next to the Whippet-free one for validation |

`quant_umbrella`, the comparison joins, the per-comparison tools and the
synthesis follow the selection. An extra left unset is switched off with the
tool it needs; one set to `true` without it, an unknown key, or switching off
all four tools stops the run before any job. Modules that need an alignment
(MAJIQ, DaPars2) still align with `hisat2: false`. The reference bundle
(`umbrella/reference/<id>/manifest.json`) lists every index, so a run that
uses any of `whippet`, `salmon` or `hisat2` builds all three indexes once;
MicroExonator alone builds none. For example, MicroExonator plus its delta
without anything else:

```yaml
umbrella_optional: {whippet: false, salmon: false, hisat2: false, multiqc: false}
```

LeafCutter uses the conda-managed Python
`leafcutter-cluster` and `leafcutter-ds` commands; it needs no separate source
checkout, R installation or scheduler submission. It is pinned in
`envs/umbrella-leafcutter.yaml`. The first real installation and comparison
remain a cluster validation step.

**Software environments.** Run with `--use-conda`. The new umbrella tools are
declared in rule-local environments: read staging (`sra-tools`, `samtools`),
alignment and coverage (HISAT2, samtools, featureCounts, rMATS, bedtools,
bedGraphToBigWig), reference and Salmon quantification (HISAT2, Salmon,
gffread), DESeq2/SUPPA2, LeafCutter, and NumPy for MicroExonator delta.
The Snakemake launcher itself remains in the base `snakemake7` environment.
Whippet is the exception: the umbrella index, quantification and delta rules
still use the configured `julia` and `whippet_bin_folder` installation. Upstream
Whippet v1.6 documents Julia 1.6.7 plus `Pkg.instantiate()` rather than a
Bioconda package, so this checkout is not yet wholly managed by a Snakemake
Conda environment. Do not remove the working cluster installation for the
first pilot. Its replacement needs a pinned, tested Julia/Whippet packaging
path. The rule-local environments have not yet been solved on the cluster.

## Targets

| Target | Reaches |
|---|---|
| `quant_microexonator` | legacy per-run quantification outputs and MicroExonator group shards, without cohort detection filtering |
| `quant_whippet` | the above, plus native Whippet quantification and Whippet shards |
| `quant_umbrella` | everything selected in `umbrella_optional` (all by default): alignment shards, Salmon shards, joins, comparisons, synthesis, and the robustly detected microexon list |
| `differential_inclusion` | MicroExonator per run, its group shards and its Whippet-free delta (`src/me_delta.py`) for every comparison in `umbrella_comparisons`; no HISAT2, Salmon or Whippet index. With `delta_method: whippet`, Whippet and its delta instead |

`quant` and `get_whippet_psi` keep their legacy meaning. In umbrella mode the
comparisons come from `umbrella_comparisons`; the legacy `whippet_delta` file
is refused.

The umbrella delta uses a fixed microexon list generated from
`Round2/TOTAL.ME_centric.txt`, which comes from the configured annotation and
ME_DB in the no-discovery run. It does **not** wait for the cohort-wide
`Report/out.robustly_detected.txt` filter. The latter is requested by
`quant_umbrella` as a parallel report, but is not a delta prerequisite.
MicroExonator delta still
applies its per-comparison read-support thresholds.

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
| `qc/{...}/runs/{run_id}/` | (`multiqc`) FastQC `fastqc_data.txt` for the first `umbrella_fastqc_reads` reads per file (default 2,000,000), HISAT2, featureCounts, Salmon and Whippet summaries: about 50 KB per run |

With `multiqc`, `multiqc/{reference_id}/{project_id}/` holds one MultiQC
report per project, built from the kept summaries, plus `runs_without_qc.txt`.
Only runs processed after this was added get summaries: a run whose QC shard
already exists has no reads left, and never triggers FastQC.

**Protected (write-protected once written)**, besides the shards above: the
per-run Whippet outputs (`umbrella/work/.../whippet/quant.{psi,jnc,gene.tpm,isoform.tpm,map}.gz`;
junction counts and gene/isoform TPM exist nowhere else), the per-run
`reads.valid` markers, the reference identity and the built HISAT2, Whippet and
Salmon indexes, and the kept QC summaries. Do not delete them: the Whippet
outputs are targets, so a missing one is remade from downloaded reads.
Comparison outputs are not protected and are recomputed freely.

**Restaging guard.** Staging refuses a run whose outputs are all kept already
(`reads.valid`, Whippet PSI, and its QC, MicroExonator and Whippet shards),
because nothing should need its reads again: something upstream changed, and
`snakemake -n -r quant_umbrella` shows what. Set `umbrella_allow_restage: true`
to reprocess on purpose. Unfinished runs are staged normally.

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

## Regrouping

`group` is the storage label a run's shards are filed under, and `batch_id`
the unit of addition; both are fixed once a run's shards exist. A run's
identity, which decides whether it is staged, aligned or quantified again, is
only its processing fields (`run_id`, sources, `layout`, `strandedness`,
`reference_id`, `project_id`, `batch_id`). So these edits reprocess nothing:

- adding metadata columns or changing their values;
- changing `sample_id` or `biological_replicate_id`;
- adding comparisons, including selections that mix groups or pick runs by
  metadata. Only the preflight, joins, comparison tools and synthesis run.

These are refused at start-up with a message saying what to do instead,
because the shard could only be rebuilt by downloading the reads again:

- moving a run to another `group` or `batch_id` (add a metadata column and
  select on it);
- adding a run to a `group`/`batch_id` whose shards exist (use a new
  `batch_id`);
- excluding (`include: false`) a run whose shards exist (use the comparison's
  `exclude`).

Coverage tracks are sums per group and batch, so they exist only for the
original groups.

**Updating a workflow started before this change.** The first start after
updating drops Snakemake's recorded params of the existing staging and shard
outputs, once (`umbrella/.run_identity_v2` marks it done), because those
params held the old whole-row hashes. Start and end times are kept. Every
preflight is rewritten, so joins, comparison tools and syntheses rerun once;
nothing per run or per shard does. Check with a dry run first.

## Operation

- **Batches.** Add runs with a new `batch_id` to the manifest and rerun
  `quant_umbrella`. A run is fingerprinted by its processing fields and a
  shard by those of the runs in its reference/project/group/batch; appending
  unrelated rows or editing labels and metadata invalidates neither. Joins and comparisons rerun only when
  their selected runs change. Keep the reference ID and existing manifest
  rows unchanged. The per-comparison outputs for a comparison expanded to
  include the new batch are recomputed, but existing per-batch shards remain
  reusable. `Report/out.robustly_detected.txt` is regenerated after its
  group-detection inputs change, so it reflects all included batches and can
  be used to assess the fixed-universe delta results afterward. An already
  materialized shard from a checkout using the older
  whole-manifest checksum guard needs review before reuse; do not delete or
  overwrite it blindly.
- **SRA downloads.** Staging an SRA run takes one `get_data` resource, as
  the legacy download rules do: `--resources get_data=10` allows ten at once.
- **Checking temporaries.** `snakemake --delete-temp-output --dry-run`
  lists what would be removed.
- **Local smoke test.** `tests/integration/mm10_simulated/run_umbrella_smoke.py`
  stages the mm10 fixture (two groups, two batches, SE and PE) against
  prebuilt full-mm10 indexes; it never builds an index. `--dry-run` needs no
  tools. Most new tools use Conda; Whippet remains the documented exception.
  rMATS and SUPPA2 packages are Linux-first.
- **Before a cluster pilot.** Record the real manifest and sample count, pin
  tool versions in `umbrella_reference: versions`, compute the reference ID,
  and check the microexon report (`whippet.microexons.report.tsv`).

## Optional modules: MAJIQ, DaPars2 and QAPA

Three opt-in tools that can run on their own, in any subset, without
MicroExonator, Whippet or the other comparison tools:

| Module | Per run (kept) | Per comparison |
|---|---|---|
| `majiq` (MAJIQ v3) | SJ file from the HISAT2 BAM | splicegraph from the selected runs, PsiCoverage per biological replicate, HET |
| `dapars2` (DaPars2) | raw coverage inside merged 3' UTR windows, mapped read count | DaPars2 by chromosome, PDUI per replicate, Welch t-test with BH (exploratory) |
| `qapa` (QAPA) | Salmon quantification against QAPA's 3' UTR library | PAU with `qapa quant`, DEXSeq site-usage test |

Config:

```yaml
umbrella_modules: [majiq, dapars2, qapa]   # default: none
umbrella_modules_mode: ingest              # default: cached_only
```

- **Targets.** `prepare_<tool>` makes the per-run caches, `quant_<tool>` adds
  every configured comparison, `quant_umbrella_modules` runs the selected
  modules plus an inventory, and `quant_umbrella` includes the selected ones
  and their inventory.
  A named target implies its tool even when it is not selected.
- **Modes.** In `cached_only`, every run a comparison uses must already have
  its caches: the workflow stops at start-up otherwise, and never stages,
  aligns or quantifies for a module. Use `ingest` once to make missing caches;
  runs whose reads are gone also need `umbrella_allow_restage: true`, for
  that invocation only. `cached_only` together with restaging is refused.
- **Identities.** Each module has a `module_reference_id` (tool, pinned
  software, module settings and the base reference ID), each cache a
  `cache_id` (the run's processing fields plus that ID), and each comparison
  output an `analysis_id` (its runs, replicates and settings). The base
  reference ID never changes, so a new module version only invalidates that
  module. Labels, comparisons and the selection list are in none of them, so
  regrouping never re-ingests.
- **Layout.** `umbrella/modules/reference/<reference_id>/<tool>/<module_reference_id>/`,
  `umbrella/modules/cache/<reference_id>/<project>/<batch>/<run>/<tool>/<cache_id>/`
  (protected), and `comparisons/<reference_id>/<project>/<comparison>/modules/<tool>/<analysis_id>/`.
- **Results.** Every comparison writes `results.tsv.gz` with shared columns
  (effect as A minus B, its definition, native statistic and its type, p and q
  only where the method gives them, and a status), the tool's native table,
  `status.json` (`ok`, `unsupported` for designs the preflight rejects,
  `native_only`, `empty` when the tool reported no features, or `no_tests`
  when nothing could be tested) and a Snakemake benchmark. Effects are not interchangeable:
  MAJIQ reports junction inclusion within an LSV, DaPars2 distal poly(A)
  usage of a 3' UTR, QAPA the usage of one poly(A) site within its gene.
  MAJIQ's effect is median PSI(A) minus median PSI(B) of one LSV connection
  over the per-replicate posterior means; p is the raw p-value of the HET test
  chosen with `majiq_het_statistic` (`ttest`, the default, `mannwhitneyu`,
  `tnom` or `infoscore`), q is Benjamini-Hochberg over the tested connections,
  and the TNOM score is kept as the native statistic.
- **Software.** DaPars2 runs from a pinned upstream commit, unpacked once;
  QAPA installs from its v1.4.1 tag in `envs/umbrella-qapa.yaml`; DEXSeq is in
  `envs/umbrella-apa-stats.yaml`. MAJIQ v3 (`rna_majiq`) is licensed: after
  registering with BioCiphers, either install it yourself (a Python 3.12
  environment with HTSlib, then `pip install ./moccasin ./majiq` from a clone
  of `majiq_academic`) and point `majiq_bin_folder` at that environment's
  `bin`, or set `majiq_source` to the clone (or a `.tar.gz` of it) to have the
  workflow build the environment. MAJIQ also needs its licence file, free for
  academic use (`https://majiq.biociphers.org/app_download/majiq_license_academic_official.lic`):
  set `majiq_license` to its path, or leave a `majiq_license*` file in `$HOME`.
  MAJIQ targets stop before any job when no licence is found. MAJIQ and DaPars2 reuse the
  existing HISAT2 index; QAPA builds its own Salmon 3' UTR index once.
- **Annotations.** DaPars2 and QAPA locate 3' UTRs from coding ends, so they
  read `apa_annotation_gtf` (default: `annotation_gtf`), which must have CDS
  lines. An exon-only GTF, such as a Whippet-oriented build, is refused, and a
  reference with fewer than `apa_min_utrs` 3' UTRs (default 5000) fails
  instead of giving an empty or tiny library. Changing it changes only the
  DaPars2 and QAPA module references (and their per-run caches), never the
  umbrella `reference_id`. `apa_transcripts: basic` (default) keeps transcripts
  tagged `basic`, QAPA's documented GENCODE input; `all` keeps every
  transcript; `megasearch` (recommended with GENCODE) first builds the
  MegaSearch APA set, once, and both tools read every transcript of it. A
  transcript is in the set when it is protein-coding with its coding end in
  the last exon, has no partial-model tag (`cds_start_NF`, `cds_end_NF`,
  `mRNA_start_NF`, `mRNA_end_NF`), is not a readthrough, and its 3' end is
  supported: a Tier-1 model (MANE, Ensembl canonical, GENCODE Primary, or
  basic with TSL 1-2), or within `apa_set_slop` nt (default 50) of a
  PolyASite 2.0 cluster supported by at least 2 protocols or of a GENCODE
  `polyA_site`. It needs `qapa_gencode_polya` and `qapa_polyasite`, the same
  files QAPA uses for its sites. On GENCODE v50 it keeps 202,029 transcripts
  of 19,111 genes (72,332 DaPars2 3' UTRs; basic gives 84,358). With any
  setting, transcripts without CDS lines or with an incomplete 3' end
  (`cds_end_NF`, `mRNA_end_NF`) are left out. DaPars2 uses protein-coding transcripts whose coding end lies in
  the last exon (3' UTRs with introns are skipped, as in QAPA). QAPA is
  annotation-only (`qapa build -N`) unless poly(A) sites are given: either
  `qapa_gencode_polya` (GENCODE's `polyAs` GTF, converted to QAPA's BED) with
  `qapa_polyasite` (a PolyASite atlas BED), QAPA's standard `-g`/`-p` route,
  or `qapa_polya_sites` for a custom BED (`-o`). Set `qapa_decoys: true` for a
  genome-decoy index. MAJIQ uses the microexon-inserted Whippet GTF, so
  microexons are in its splicegraph.
- **Options.** `dapars2_coverage_threshold` (default 10);
  `majiq_update_args`, `majiq_psicov_args`, `majiq_heterogen_args` replace
  the default arguments of those MAJIQ commands (each command's `--help` is
  written to the comparison log); `majiq_strandness` overrides the manifest's
  strandedness for MAJIQ; `majiq_het_statistic` picks the HET test for p and q.
  Each run's SJ file is named by run (`--prefix`), and the technical runs of
  a biological replicate are summed into one PsiCoverage with rna_majiq's
  `PsiCoverage.sum`, run with the `python` of the MAJIQ installation. Each SJ
  cache records the MAJIQ version that wrote it; a comparison refuses SJ files
  from another version than the installed one, and reports the version in its
  results and `status.json`.

## Known limits

- Whippet psi path columns (`Inc_Paths`, `Exc_Paths`, `Edges`) are not kept.
- Technical runs are summed in count tables, but the delta tools still see
  each run's PSI table separately (flagged as `technical_run_replicates`).
- Junction capture uses splice-site sets per chromosome and strand, not per
  gene.
- The planned Salmon versus Whippet TPM concordance is not implemented.
- The MAJIQ v3 commands were checked against the `rna_majiq` source
  (majiq_academic main, 2026-10-01) and tested with argument-checking fakes,
  but not yet run on real data; check the first comparison log.
- Module-only targets on a brand-new reference also build the Whippet and
  Salmon indexes, because HISAT2 waits for the full reference manifest; its
  inputs cannot change without re-aligning existing runs.

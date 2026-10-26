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

## Optional modules: MAJIQ, DaPars2, QAPA, vast-tools, conservative rMATS, PSI-Sigma and circRNA

Six opt-in tools that can run on their own, in any subset, without
MicroExonator, Whippet or the other comparison tools:

| Module | Per run (kept) | Per comparison |
|---|---|---|
| `majiq` (MAJIQ v3) | SJ file from the HISAT2 BAM | splicegraph from the selected runs, PsiCoverage per biological replicate, HET |
| `dapars2` (DaPars2) | raw coverage inside merged 3' UTR windows, mapped read count | DaPars2 by chromosome, PDUI per replicate, Welch t-test with BH (exploratory) |
| `qapa` (QAPA) | Salmon quantification against QAPA's 3' UTR library | PAU with `qapa quant`, DEXSeq site-usage test |
| `vasttools` (vast-tools) | compressed `to_combine` cache from staged reads and VASTDB indexes | technical-run merge, combine, tidy, Bayesian diff or Mann-Whitney, microexon support mapping |
| `vasttools_turbo` (vast-tools, chunked) | the same cache as `vasttools`, made by independent chunk jobs and pooled with `vast-tools merge` | identical to `vasttools` (same compare step) |
| `psisigma` (PSI-Sigma v2.3) | junction counts and coverage at retained-intron positions, from the HISAT2 BAM | replicate sums, PSI-Sigma database, PSI and filter steps, Welch t-test with BH |
| `rmats_conservative` (rMATS-turbo 4.3.0) | `prep` summary from the HISAT2 BAM, without `--novelSS` | `inte` + `post` on the kept summaries, JC events normalized, native tables archived |

Config:

```yaml
umbrella_modules: [majiq, dapars2, qapa, vasttools, vasttools_turbo, rmats_conservative, psisigma]   # any subset; default: none
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

### vast-tools

Select `vasttools` explicitly; it is not added to existing runs by default.
`prepare_vasttools` produces its retained caches and `quant_vasttools` adds
the configured comparisons. The standard `quant_umbrella` target includes it
only when selected in `umbrella_modules`. First ingestion needs staged reads;
if those reads were already deleted, explicitly allow restaging using the
same guard as the other modules. Later cached-only comparisons need no FASTQs.

```yaml
umbrella_modules: [vasttools]
umbrella_modules_mode: ingest
vasttools_species: hg38
```

The runtime is managed by `envs/umbrella-vasttools.yaml` (Bowtie 1, Perl, R
and dependencies). Snakemake downloads the pinned MIT-licensed tool source
and the matching prebuilt VASTDB library once per run directory; no manual
software checkout or genome re-index is needed. `vasttools_species` supports
hg38 (default), hg19, mm10, mm9 and dm6; explicitly match it to the input
assembly. Downloads run locally through Snakemake's local rules. hg38's
compressed archive is about 5 GB; allow additional space for the unpacked
permanent library and the temporary split reads. A mirror can be supplied
with `vasttools_vastdb_url`, but must serve the same library release/content.

For each comparison:

- Technical runs are merged within biological replicates, never across them.
- `tidy` requires coverage in each group, with `--noVLOW --p_IR --noB3`.
- When either group has at most `vasttools_diff_max_n` replicates (default 5),
  calls require a tidy pass and Bayesian MV >= `vasttools_min_dpsi` (0.1), at
  `vasttools_diff_prob` (0.95). Minimum corrected reads are
  `vasttools_min_reads` (10); minimum covered samples per group are
  `vasttools_min_samples` (2), bounded by group size.
- When both groups exceed that limit, calls use two-sided Mann-Whitney on
  tidy PSI values (at least three valid values per group), BH q below
  `vasttools_mwu_q` (0.05), and absolute mean dPSI >=
  `vasttools_mwu_min_dpsi` (0.15). The unused Bayesian calculation is skipped;
  `diff.tsv.gz` contains a clear not-run marker. Effects are always A minus B
  as fractions, not percentages. Native inclusion/tidy tables retain percentages.

`microexons.tsv.gz` checks exact coordinates and strand against the installed
FULL and MIC catalogue templates, resolving old-to-new event IDs. Statuses:

| Status | Interpretation |
|---|---|
| `called` / `tested` / `filtered` | A reported matching exon event, with the stated comparison outcome |
| `not_quantified` | The event is reported but has no PSI in the compared samples |
| `not_reported` | Present in the installed catalogue templates but not emitted as a matching exon event; may be a legacy MIC entry outside the active FULL catalogue |
| `mapping_ambiguous` | Coordinate candidate exists but strand cannot be resolved safely |
| `not_in_vastdb` | No exact exon-coordinate/strand match in the inspected installed templates |

Only measured matching events provide comparison evidence. None of the last
four statuses is evidence contradicting MicroExonator. Catalogue membership
does not imply that the installed tool can quantify that event. Canonical
microexon evidence remains the standalone MicroExonator Delta result; this
module supplies orthogonal support where VASTDB measures a matching event.

The guarded two-pass launcher accepts an explicit `vasttools` entry in a
complement execution's `tools`; its default tool selection is unchanged.
Comparison implementation changes alter the analysis ID, not per-run cache
IDs. Do not weaken packet-identity guards to apply a code update to an
already-frozen experimental packet; use its documented explicit revision.

Optional native smoke test (reuse an installed source and matching VASTDB):

```bash
VASTTOOLS_TEST_DIR=/path/to/vast-tools \
VASTTOOLS_TEST_DB=/path/to/VASTDB \
VASTTOOLS_TEST_SPECIES=dm6 \
python -m unittest tests.test_run_vasttools.RealVastToolsTests -v
```

Use the conda environment's Python and put its `bin` on PATH. The test samples
at most 60 junction events and 24,000 paired fragments per processing unit;
it tests seven units, including a split technical replicate, using existing
indexes. It tests tool wiring and catalogue statuses, not sensitivity on
biological RNA-seq. Native `vast-tools` uses Bash process substitution;
restrictive filesystem sandboxes can block Bowtie's `/dev/fd` writes.

### vast-tools-turbo

`vasttools_turbo` makes the same per-run cache as `vasttools`, but splits a run
into chunks that align as independent cluster jobs. It is a separate optional
module, like `rmats_conservative`: the native `vasttools` rules, caches and
results are untouched, and both can be selected in one execution to compare
them.

```yaml
umbrella_modules: [vasttools, vasttools_turbo]   # native and turbo side by side
umbrella_modules_mode: ingest
vasttools_species: hg38
vasttools_turbo_chunks: 8      # chunks per run (changes the cache identity), default 8
vasttools_turbo_threads: 4     # Bowtie threads per chunk job (no identity change), default 4
```

`prepare_vasttools_turbo` makes the caches and `quant_vasttools_turbo` adds the
comparisons. Both modules share the pinned vast-tools source, VASTDB and the
staged reads, so nothing is downloaded or staged twice.

Per run, three kinds of jobs (`src/run_vasttools_turbo.py`):

1. `umbrella_vasttools_turbo_split` decides strandedness exactly as
   `vast-tools align` does (first 500K reads against VASTDB's mRNA library,
   0.35 / 0.65 thresholds, each mate separately), then deals the reads
   round-robin into N chunks. Mates always stay in the same chunk.
2. `umbrella_vasttools_turbo_chunk` (N jobs) runs `vast-tools align` on one
   chunk with that decision forced (`--ns`, or `--rc1/--nrc1` and
   `--rc2/--nrc2`), so no chunk can decide differently from its siblings.
3. `umbrella_vasttools_turbo_merge` pools the chunks with `vast-tools merge`
   (the tool's own way to deepen a sample) and writes `to_combine.tar` and
   `cache.json` in the native layout. It refuses to write a cache if a chunk
   is missing, a chunk's strandedness or read count disagrees with the split,
   a required file is missing, or the merged `.eej2` is empty.

The chunk FASTQs and chunk outputs are temporary. `compare` is the native
step (`src/run_vasttools.py compare`), so results have the same columns; the
`tool` column says `vasttools_turbo`.

Differences from native, by design:

- Raw counts are the same. Corrected counts and PSI can differ in the last
  rounding digits, because vast-tools rounds per file before `merge` sums
  them. Measured on the 250,000-read Hsa test FASTQ (hg19 VASTDB, 2 and 3
  chunks, single-end and its halves as paired-end; opt-in
  `RealTurboEquivalenceTests`): single-end gives the same events, identical raw
  EEJ counts, and every number within 0.06 (vast-tools' rounding), with
  identical labels. Paired-end differs a little more: `Trim.pl` numbers its
  50-mers with one counter over the whole run, which sets the order of a pair's
  fragments after sorting and so which of them vast-tools' adjacent-line logic
  counts once; chunks renumber, so up to 3.6% of EEJ rows move by at most 2
  reads (totals within 0.1%; intron-retention totals within 0.5%, one intron
  corrected count by 5%). A native run repeated twice is identical, so these
  differences come from chunking, not noise. Compare the two modules on your
  own paired-end data before relying on turbo for it.
- vast-tools `merge` writes the literal `S` in the exskX complexity column,
  losing the `C` annotation of 2,900 of 18,612 events; turbo restores it from
  the chunks (it does not depend on reads, and chunks that disagree are
  refused). Native technical-run merges in `compare` still have this upstream
  limitation; it is left as is.
- Only `vasttools_ir_version: 2` is supported (`merge` does not rebuild the
  version 1 intron retention summary); version 1 is refused at start-up.
- With `vasttools_expr: true` the cRPKM files are merged, but vast-tools'
  `.3bias` file is not (`compare` uses neither).
- Each chunk loads the Bowtie indexes again (about 26 s of CPU for the
  COMBI library), so total CPU grows slightly while wall time per run drops
  to roughly one chunk. Each run needs at least 1000 reads per chunk.
- The cache identity includes the chunk count; changing
  `vasttools_turbo_chunks` makes new caches.

To compare the two after a run, read the benchmark files:
`umbrella/modules/benchmarks/<reference>/<project>/<batch>/<run>.vasttools.<id>.tsv`
for native, and the `...vasttools_turbo_split`, `..._chunk<k>` and
`..._merge` files for turbo (wall time of a run is about split + the slowest
chunk + merge). The two comparison folders under `modules/vasttools/` and
`modules/vasttools_turbo/` can be diffed on `results.tsv.gz`.

Cluster resources are in `Examples/Cluster_config/pbs/cluster.PBS.json`
(`umbrella_vasttools_turbo_split`, `_chunk`, `_merge`, `_compare`). Chunk jobs
use 4 cores and 32 GB: memory is set by the Bowtie index, not by the thread
count.

### Conservative rMATS

`rmats_conservative` is a second rMATS profile next to the original rMATS
comparison tool, which is unchanged and still runs with `--novelSS`. Select it
explicitly (`umbrella_modules: [rmats_conservative]`); `prepare_rmats_conservative`
makes the per-run `run.rmats` caches and `quant_rmats_conservative` adds every
comparison.

- **No novel splice sites.** `prep` and `post` run without `--novelSS`, so
  events are built from annotated splice sites only. This is not a guarantee of
  correctness: rMATS can still report novel combinations of annotated sites.
- **Same library policy as the original route:** `fr-unstranded`,
  `--variable-read-length`; when prep read lengths differ, `post` uses the
  maximum (recorded with a warning in `status.json`).
- **One tool family.** Both profiles are rMATS; the module inventory lists
  `rmats_conservative` under the `rmats` family. Count them as one tool, not
  two independent votes.
- **Effect** is rMATS' `IncLevelDifference`, PSI(A) minus PSI(B), with its
  likelihood-ratio p-value and per-event-class FDR. `native.tar.gz` keeps the
  native JC/JCEC tables and individual counts.
- **Replication.** A comparison needs two biological replicates per side, and
  one processing unit per replicate: merge technical runs with `library_id`
  first. Otherwise the module reports `unsupported` instead of treating runs as
  replicates.
- **Reuse.** The per-run prep reads the shared temporary HISAT2 BAM in
  `ingest` mode; later comparisons in `cached_only` mode reuse `run.rmats`
  without reads or BAMs. A cache whose reference, module profile, cache ID or
  layout differs from the comparison is refused.

### PSI-Sigma

PSI-Sigma v2.3 (Lin & Krainer, *Bioinformatics* 2019) quantifies exon
skipping, alternative 5'/3' sites, mutually exclusive exons and intron
retention. Select it explicitly (`umbrella_modules: [psisigma]`);
`prepare_psisigma` makes the per-run caches and `quant_psisigma` adds every
comparison. Its licence allows non-commercial use by individuals at academic or
non-profit institutions; the workflow downloads the pinned source
(`umbrella_psisigma_software`), checks every script against its SHA-256 and
keeps PSI-Sigma's LICENSE next to it.

- **Installation downloads from CPAN.** PSI-Sigma's statistics use PDL,
  PDL::Stats, PDL::GSL and Statistics::Multtest, which conda does not carry.
  `src/install_psisigma.py` builds pinned, checksummed CPAN releases once into
  `umbrella/modules/software/psisigma/<commit>/perl5`, against the Perl and GSL
  of `envs/umbrella-psisigma.yaml` (cpanm also fetches their dependencies). It
  runs on the login node, like the other downloads.
- **HISAT2 instead of STAR.** PSI-Sigma was written for STAR. The module feeds
  it the shared HISAT2 alignment with three adaptations:
  - unique alignments are HISAT2's `NH:i:1`, given to PSI-Sigma as MAPQ 255
    (its unique threshold) in the stream; no BAM is rewritten;
  - spliced records are grouped by read name (`samtools collate`) before
    PSI-Sigma's own junction counter, so mates crossing one junction count once,
    as PSI-Sigma's `samtools sort -n` does;
  - non-canonical junctions (not GT-AG, GC-AG or AT-AC) are dropped, as
    PSI-Sigma's STAR protocol does with `--outFilterIntronMotifs
    RemoveNoncanonical` (`psisigma_canonical_only: true`, the default). This
    filters junction counts, not whole alignments, so it approximates STAR.
  STAR's two-pass mode has no equivalent here. Soft clipping differs between
  aligners and PSI-Sigma ignores clipped reads for intron coverage, so each
  run's `qc.json` records its soft-clipped fraction: check it is similar
  between groups.
- **Intron retention without keeping BAMs.** Only PSI-Sigma's default
  `--irmode 0` is allowed (`psisigma_irmode` must be 0). In that mode the
  retained-intron events and their five coverage positions come from the
  annotation alone, so each run keeps its coverage at those positions, counted
  with PSI-Sigma's read rules. A comparison sums a replicate's runs (exactly
  what a merged BAM gives), takes PSI-Sigma's median per intron, and checks
  that its own database has the same retained introns. Comparisons never read
  a BAM, so `cached_only` regrouping keeps intron retention.
- **Annotation.** The umbrella's `annotation_gtf` with the ME_DB microexons
  inserted, restricted to primary chromosomes (PSI-Sigma splits IDs on `_`),
  with gene and transcript lines added where missing and transcripts sorted
  before their exons. PSI-Sigma quantifies annotated events only: microexons
  absent from that annotation are not seen, and very short exons align poorly,
  so their inclusion tends to read low.
- **Statistics.** PSI-Sigma's `--variance 0` (the default) is labelled
  Student's t but runs Welch's t-test; the FDR is Benjamini-Hochberg over all
  events (`--adjp 1`). PSI-Sigma needs two biological replicates per side.
  When it cannot test an event (equal group means, or identical replicates) it
  prints dPSI 0; the module reports that effect as blank and the event as
  `untested`.
- **Effect** is PSI(A) minus PSI(B) as a fraction: B is passed as PSI-Sigma's
  `groupa` (N) and A as `groupb` (T), so its dPSI (T - N) is already A - B.
  `native.tar.gz` keeps PSI-Sigma's database, unfiltered and sorted tables.
- **Options** (defaults are PSI-Sigma's driver defaults): `psisigma_nread`
  (10), `psisigma_skipratio` (0.05), `psisigma_irrange` (5),
  `psisigma_variance` (0), `psisigma_irclean` (0), `psisigma_canonical_only`
  (true).

### circRNA modules: CIRIquant, CIRCexplorer2 and DCC

Three opt-in modules quantify and detect back-spliced junctions (BSJ):

| Module | Per run (kept) | Per comparison |
|---|---|---|
| `ciriquant` (CIRIquant 1.1.3) | `counts.tsv.gz`: BSJ and FSJ read pairs for every catalogue circRNA, from the staged FASTQs and the shared HISAT2 BAM | replicate BSJ/FSJ matrices, CIRIquant's junction ratio 2*BSJ/(2*BSJ+FSJ), A minus B of the replicate means |
| `circexplorer2` (2.3.8) | `calls.tsv.gz` from the shared STAR chimeric cache | detection support per replicate |
| `dcc` (0.5.0) | none (reads the shared STAR caches) | DCC detection with its replicate filter |

All comparison outputs are descriptive or detection-only (`status: no_tests`): no
differential test is claimed for circRNAs in this release. Circular RNAs are depleted
from poly(A)-selected libraries; select comparisons with documented rRNA-depleted
libraries for inference, and record the library evidence in the selection.

**Per-comparison scope.** These modules default to `umbrella_module_scope: selected`:
they run only on the comparisons listed for them in `umbrella_module_selection`, and
make per-run caches only for those comparisons' runs, in the same pass as the modules
that run on every comparison. The selection is provenance, not identity: changing it
never changes a cache, reference or analysis ID.

```yaml
umbrella_modules: [vasttools, rmats_conservative, psisigma, ciriquant, circexplorer2, dcc]
umbrella_modules_mode: ingest
umbrella_module_selection: circrna_selection.tsv
# umbrella_module_scope: {circexplorer2: all}   # optional override per module
circrna_catalogue_sources:
  - {path: Homo-sapiens.txt, format: circpedia, tissue: Central nervous system, name: CIRCpedia_v3_CNS}
  - {path: human_bed_v3.0.txt, format: circatlas, role: annotate, name: circAtlas_v3}
  - {path: hsa_hg19_circRNA.txt, format: circbase, assembly: hg19, role: annotate, name: circBase}
circrna_liftover_chain: hg19ToHg38.over.chain.gz
circrna_max_span: 100000
```

`circrna_selection.tsv` (tab-separated; `#` lines are comments) needs the columns
`comparison_id, tool, role, library_evidence, decided_by, approved_by, reason`. Rows
without `approved_by`, unknown comparisons and modules whose scope is `all` are refused.
`role` (e.g. `test`, `polya_reference`, `library_unknown`) and the evidence are copied
into each comparison's `status.json` and the module inventory.

- **Catalogue sources.** `format` reads each database as it is published: `bed6`
  (0-based start), `circatlas` (circAtlas 3.0 `human_bed_v3.0.txt`, whose Start is
  1-based), `circbase` (`hsa_hg19_circRNA.txt`, 0-based) and `circpedia` (CIRCpedia v3
  `Homo-sapiens.txt`; `tissue` keeps rows whose long-read evidence names that tissue).
  `role: quantify` sources make the catalogue; `role: annotate` sources only add their
  names to the same BSJs. Whole databases are far too large for CIRIquant (circAtlas 3.0
  human: 32 Gb of pseudo-reference), so quantify a tissue-relevant subset.
- **Catalogue.** Sources are on `circrna_assembly` (hg38 by default) or lifted to it: each BSJ boundary is lifted separately through the UCSC
  chain, both must land on the same contig and orientation, and failed rows are kept
  beside the lifted file with their reason. The build merges names of the same BSJ,
  drops rows outside the genome, ambiguous strands and mostly-N spans, reports splice
  motifs, and refuses catalogues whose CIRIquant pseudo-reference (each circRNA's
  genomic span written twice) exceeds `circrna_max_pseudo_bases` (2e9). `circrna_max_span`
  (0 = no cap) is optional. The frozen catalogue and its content-derived ID live in
  `umbrella/modules/reference/<reference_id>/circrna_catalogue/<catalogue_id>/`; its ID
  enters CIRIquant's identity only.
- **CIRIquant** needs paired-end runs; a selection asking it for a single-end or mixed
  comparison stops at start-up. Its library type follows the umbrella strandedness
  (`auto` reads the inferred one). `src/run_ciriquant.py` calls CIRIquant's own index,
  realignment and BSJ/FSJ counting functions; it replaces only two helpers that still
  carry Python 2 arithmetic in 1.1.3 (`load_fai`, `grouper`), and a test fails if
  upstream fixes them. The pseudo-circular index is built once per catalogue
  (`umbrella_ciriquant_index`, with CIRIquant's own `generate_index` and `build_index`)
  instead of once per run; counts are identical (native test). Counts are exported for every catalogue circRNA (CIRIquant's own
  GTF omits BSJ = 0), with `not_indexed` and `no_support` kept distinct from zero.
- **Shared STAR.** One STAR 2.7.11b run per run (chimeric junctions only, no BAM,
  DCC's joint-mapping settings) serves both CIRCexplorer2 and DCC; report their
  agreement as shared-evidence support, not two independent discoveries. The STAR
  index is built once from the umbrella genome and annotation (`circrna_star_sjdb_overhang`,
  default 100). CIRCexplorer2 calls are computed from the kept STAR cache, so they can
  be made in `cached_only`; DCC runs per comparison on the kept STAR caches, with its
  CircSkip step switched off.
- **Targets.** `prepare_<tool>` and `quant_<tool>` for `ciriquant`, `circexplorer2` and
  `dcc`; `quant_umbrella_modules` includes them when selected.

Tests: `tests/test_circrna_catalogue.py`, `tests/test_run_ciriquant.py`,
`tests/test_circrna_discovery.py`, `tests/test_circrna_workflow_dag.py`,
`tests/test_umbrella_module_scope.py` and `tests/test_circrna_workflow_native.py` (a real
Snakemake run on a synthetic genome). The native tests need conda prefixes made from the
two environment files:

```bash
CIRCRNA_CIRIQUANT_ENV=/path/to/ciriquant-env CIRCRNA_DISCOVERY_ENV=/path/to/discovery-env MICROEXONATOR_SNAKEMAKE=$(which snakemake) python -m unittest tests.test_run_ciriquant tests.test_circrna_discovery tests.test_circrna_workflow_native -v
```

Real STAR runs only on Linux in these tests (STAR's macOS builds cannot spawn
`--readFilesCommand`; the module streams gzipped reads through named pipes instead).

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

### vast-tools speed work (pending a pin bump)

The refurbished vast-tools branches (`perf/*` in `vast-tools-refurbished`) give
byte-identical `to_combine` content with `align` about 18% faster at `-c 8`
(concurrent lanes), plus an optional Rust trimmer. They are local and
unpushed, so `VASTTOOLS_COMMIT` in `src/umbrella_modules.py` is deliberately
unchanged: a bump must point at a fetchable commit, and it changes
`module_reference_id`, so every vast-tools cache is rebuilt. Do it as one
deliberate step once the commit is published. `vasttools_turbo` (above) is the
workflow-level speed-up and runs on the pinned commit.

- The cache gzip in `src/run_vasttools.py` only compresses the small
  `to_combine` files, so switching it to pigz would gain nothing measurable.
# Optional circRNA differential comparisons

After circRNA counting, enable comparison-only inference with:

```yaml
circrna_differential: true
circrna_hisat_delta: true
circrna_delta_min_samples: 2
circrna_delta_min_reads: 5
circrna_delta_draws: 10000
circrna_delta_threshold: 0.1
```

The existing selected circRNA targets then include `umbrella_circrna_differential`.
Add that rule to any `--allowed-rules` list. The conda environment is
`envs/umbrella-circrna-diff.yaml`; the PBS example includes a one-core, 8 GB entry.
Outputs live beneath each tool/comparison's `differential/<settings-id>/` directory.
Changing these settings does not change counting-cache identities.

`abundance.tsv.gz` uses edgeR quasi-likelihood inference on integer BSJ counts,
TMMwsp normalization and BH adjustment within this tool/comparison. This follows
the recommended replicated count-model family, but does not reproduce native
CIRI_DE_replicate's gene-derived normalization. The tested universe consists of
the retained native calls; discovery/annotation filtering limits it. Missing
counts are not converted to zeros. With inadequate support, `status.tsv` reports
`insufficient_data` rather than manufacturing significant events: fewer than two
replicates per side, fewer than `circrna_abundance_min_events` (default 10) events
before or after `filterByExpr`, or a replicate without BSJ reads. An edgeR fit error
is `fit_failed` with its message; R warnings are kept as `warning` lines. Read
`events_tested` and `normalization.tsv` with the results: when filtering leaves few
events, most of them changed, TMMwsp and the dispersion estimate absorb the change.

`hisat_delta.tsv.gz` is a separate **experimental** Whippet-style Bayesian
comparison of BSJ/(BSJ+linear support), independently at the donor and acceptor.
Linear support comes from unique-fragment HISAT2 junction counts sharing the
same chromosome, strand and exact splice boundary. It is not native CIRIquant's
2*BSJ/(2*BSJ+FSJ), nor an estimate of complete circular transcript abundance.
The two boundaries are not independent replication and must not count as two
tool votes. `hisat_evidence.tsv.gz` retains the per-replicate numerator and
denominator. Cross-aligner count differences and possible shared fragment
evidence mean the posterior confidence is not calibrated FDR.
In count simulations (`python -m tests.simulate_circrna_delta`; 120 events, 3 vs 3,
P(|delta| > 0.1) >= 0.95) the pooled moment-fit posterior, like MicroExonator Delta,
made no calls under a homogeneous null, but 7/120 under a null with heterogeneous
replicates, and found 14-18/120 moderate changes (0.3 vs 0.5). Treat calls from
overdispersed comparisons with caution.

The Bayesian route requires **already retained** junction shards at
`junctions/<reference>/<project>/<group>/<batch>.junctions.tsv.gz`. Missing
shards stop the DAG; enabling this option never implicitly regenerates HISAT2
alignments. For a new cohort, finish retaining those shards first, then run the
comparison from caches. Set `circrna_hisat_delta: false` to run abundance only.
The current module-only rerun does not itself create these group/batch shards.

No one-to-one inference: at least two eligible biological replicates per group
are required. Technical units are collapsed first. Opposite-strand circles,
unknown library protocols and detector-only novel events need separate review.
Assess tool agreement within the jointly measurable event universe, separately
for abundance and usage. Shared STAR/HISAT2 evidence is correlated; consensus
alone does not establish specificity.

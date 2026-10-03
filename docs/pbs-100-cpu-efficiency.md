# Running `quant_umbrella` efficiently on PBS with a 100-CPU cap

Notes for running the umbrella on a PBS cluster where one user can hold at most
100 CPUs at a time (the Donnelly Centre limit), with a Snakemake profile that
reads `cluster.PBS.json`. The command these notes assume:

```bash
snakemake -s MicroExonator.smk quant_umbrella --profile base \
  --use-conda --conda-frontend conda --conda-prefix $C \
  -k --rerun-incomplete -j 50 --resources get_data=8
```

The rule set uses `retries:`, so this is Snakemake 7 or later. Line numbers
refer to the `refurbishment/umbrella` branch.

## The short version

1. Tell Snakemake about the CPU cap with a `cpus` resource. `-j` counts jobs,
   not cores.
2. Lower `ppn`/`threads` on the staging jobs. They mostly wait on the network
   but each one holds 4 cores.
3. In `cluster.PBS.json`, fix the entries that request more or fewer CPUs than
   the rule uses, and give the small legacy rules realistic memory and walltime.
4. In the profile, add a `cluster-status` script and a `latency-wait` setting
   so jobs PBS kills don't hang the run.
5. Time one representative run at different thread counts before scaling up.

## 1. `-j 50` limits jobs, not CPUs

In cluster mode, `-j 50` means "at most 50 submitted jobs". Each job requests
its own `ppn` from `cluster.PBS.json`:

| Per-run job (quant_umbrella) | ppn | mem | Most at once in 100 CPUs |
|---|---|---|---|
| `umbrella_hisat2` | 8 | 16 GB | 12 |
| `Round2_bowtie_to_tags` (legacy ME quant) | 8 | 2 GB | 12 |
| `Round1_bwa_mem_to_tags` | 5 | 10 GB | 20 |
| `umbrella_stage_reads` / `umbrella_stage_single` | 4 | 3 / 8 GB | 25 |
| `umbrella_legacy_fastq` | 4 | 2 GB | 25 |
| `umbrella_featurecounts` | 4 | 3 GB | 25 |
| `umbrella_rmats_prep` | 4 | 6 GB | 25 |
| `umbrella_salmon_quant` | 4 | 24 GB | 25 |
| `umbrella_whippet_quant` | 1 | 12 GB | 100 |
| `umbrella_junctions`, `umbrella_coverage_run`, per-run legacy scripts | 1 | 2–10 GB | 100 |

So 50 jobs can ask for anywhere from 50 to 400 CPUs. PBS enforces the 100-CPU
limit and holds the rest in the queue, which causes three problems:

- **Snakemake loses control of the order.** Queued jobs start in PBS order, not
  Snakemake's (`priority:` stops having an effect). This matters most for the
  temporary BAM and FASTQ files. A run's reads and its `aligned.bam` are
  deleted only after every consumer has finished (junctions, featureCounts,
  rMATS prep, coverage, Salmon, Whippet, MicroExonator). If the consumers sit
  behind other runs' alignments, temporaries pile up on scratch.
- **Waiting jobs use up `-j` slots.** If 40 of the 50 slots are 8-CPU jobs
  waiting in PBS, Snakemake can't submit the cheap 1-CPU jobs that would fill
  the remaining cores.
- **Queue limits.** Some PBS setups also cap how many jobs one user can have
  queued. Check whether Donnelly does.

**Fix: count CPUs as a Snakemake resource.** Snakemake applies global
`--resources` limits in cluster mode too. In `base/config.yaml`:

```yaml
default-resources:
  - cpus=1
set-resources:
  # keep in step with ppn in cluster.PBS.json
  - umbrella_hisat2:cpus=8
  - Round2_bowtie_to_tags:cpus=8
  - Round1_bwa_mem_to_tags:cpus=5
  - umbrella_stage_reads:cpus=4
  - umbrella_stage_single:cpus=4
  - umbrella_legacy_fastq:cpus=4
  - umbrella_featurecounts:cpus=4
  - umbrella_rmats_prep:cpus=4
  - umbrella_salmon_quant:cpus=4
  - umbrella_fastqc:cpus=2
  - umbrella_hisat2_index:cpus=8
  - umbrella_salmon_index:cpus=8
  - umbrella_rmats_post:cpus=4
  - umbrella_leafcutter:cpus=4
  # add the module rules you enable (umbrella_majiq_*, umbrella_qapa_*, ...)
```

Then put the cap on the command line. A command-line `--resources` replaces the
profile's `resources:` entry rather than merging with it, so list both
resources there:

```bash
... -k --rerun-incomplete -j 150 --resources get_data=8 cpus=96
```

With `cpus` doing the limiting, `-j` can go above 100 so that 1-CPU jobs always
have a slot. Keep the `cpus` cap slightly below 100 if anything else of yours
(an interactive session, or Snakemake itself running inside a job) uses part of
the allowance. Check the setup with `snakemake ... -n` first, and confirm that
`qstat -u $USER` never shows more than about 100 cores running.

The `cpus` values must match `ppn` in `cluster.PBS.json`. If they drift apart,
Snakemake's count and PBS's count disagree.

## 2. Staging (downloads) holds more cores than it uses

The rule's own comment (`rules/umbrella_inputs.smk:59`) notes that while
downloading, staging uses 24–59 % of **one** core. Only fasterq-dump/pigz at
the end use all 4. With `--resources get_data=8`, that's 32 cores (a third of
the cap) held mostly idle by downloads. Local (non-SRA) runs also take 4 cores
each, without a `get_data` slot.

Options, in order of effort:

- Set `ppn: 2` for `umbrella_stage_reads` and `umbrella_stage_single` in
  `cluster.PBS.json`, and `threads: 2` in the rule. Expect slower extraction
  and compression in return for more alignments running.
- Keep `get_data` around 6–8. More concurrent downloads mostly fill scratch
  with FASTQs waiting for alignment slots.
- If Donnelly lets jobs run on the login or a transfer node, make staging a
  `localrule` so it doesn't use the PBS allowance at all (it still respects
  `get_data`).

## 3. Thread counts on the per-run tools

With a hard CPU cap, total throughput across samples matters more than how
fast one sample finishes. A thread is only worth giving if the tool actually
speeds up with it.

- **`Round2_bowtie_to_tags`** (`rules/Round2.smk:374`) is
  `gzip -dc | awk | bowtie -p 8 | awk | python3 select_tag_alignments.py`.
  Bowtie against the small tag index is fast, so the single-threaded
  gzip/awk/Python stages probably limit it. If 3 threads give about the same
  speed as 8, you can run 33 at once instead of 12. Change `threads:` in the
  rule, `ppn` in the JSON and the `cpus` override together.
- **`Round1_bwa_mem_to_tags`** (5 threads, piped to awk) is in the same
  situation, if `quant_umbrella` pulls in discovery.
- **`umbrella_hisat2`** scales well up to about 8 threads, and `samtools sort -@`
  reuses them, so 8 is reasonable. It is the job that occupies the most
  core-hours, though, so measure 4 vs 8 as well.
- **Whippet, junctions, coverage and the per-run scripts** are single-threaded.
  They are cheap in CPUs, so per-node memory, not the CPU cap, decides how many
  can run (`umbrella_whippet_quant` asks for 12 GB).

How to measure: pick one typical run, use `--forcerun <rule>` with a few
`threads` values (or add a `benchmark:` directive), and compare wall time.
Pick the thread count where doubling threads no longer roughly halves the
time.

## 4. `cluster.PBS.json` fixes

Compared against the `threads:` the rules actually use:

**Request more CPUs than they use:**

| Rule | ppn | Rule threads | Fix |
|---|---|---|---|
| `bowtie_genome_index` | 8 | 1 | `ppn: 1`, or add `--threads` to the rule |
| `Output` | 2 | 1 | `ppn: 1` |

**Use more CPUs than they request** (no JSON entry, so they get `ppn: 1`; only
matters if you use these aligners): `total_STAR_to_Genome` (5),
`total_olego_to_Genome` (10), `total_tophat_to_Genome` (5), `mv_STAR` (5).
These oversubscribe the node, or get killed where CPU use is enforced.

**Default memory and walltime too large.** About 30 legacy entries fall back to
`__default__` (10 GB, 24 h), including:

- `*_alingment_pre_processing`
- `ME_psi_to_quant`
- `SamToBam`, `BamIndex`
- `correct_quant`, `get_PSI_sparse_quants_*`
- the detection filters and the delta rules

Most finish in minutes. PBS backfills short jobs into gaps sooner, and a 10 GB
request blocks scheduling wherever memory is also limited. The umbrella entries
already use realistic values (1–8 h). Something like 2–4 GB and 2–4 h for the
small legacy steps would bring them in line. Check a few `logs/*.err` or
`qstat -f` records for real peak memory before tightening.

## 5. Profile settings

Things to have in `base/config.yaml` (it is the same file as above):

- **`cluster-status: <script>`.** Without one, Snakemake doesn't notice when PBS
  kills a job for exceeding walltime or memory. The job just never finishes and
  `-k` keeps waiting. A short script that maps `qstat -f <jobid>` or
  `tracejob <jobid>` to `running` / `success` / `failed` is enough. The submit
  command must print the job ID (plain `qsub` does).
- **`latency-wait: 60`.** Shared filesystems can be slow to show new files, and
  without this you get false "missing output" failures.
- **`max-jobs-per-second: 5` and `max-status-checks-per-second: 1`.** These keep
  the PBS server responsive when many short jobs are submitted at once.
- **Check the resource line of `cluster:`.**
  - Torque uses `-l nodes={cluster.nodes}:ppn={cluster.ppn},mem={cluster.mem},walltime={cluster.walltime}`.
  - PBS Pro uses `-l select=1:ncpus={cluster.ppn}:mem={cluster.mem} -l walltime={cluster.walltime}`.

  A job whose `ppn` never reaches PBS silently runs with one core.
- **`--rerun-incomplete`** (already in the command) is right here: a job PBS
  killed leaves incomplete outputs behind.

## 6. Where Snakemake runs, and scratch

- Running Snakemake itself inside a batch job uses one CPU of the 100 for the
  whole run. If the login node allows long-lived processes, run it in `tmux`
  there instead.
- Set `umbrella_stage_tmpdir` to node-local scratch if the nodes have it.
  Uncompressed reads are written there during staging, and 25 jobs doing this
  on the shared filesystem at once slow each other down.
- Watch scratch during the first batch: staged FASTQs and `aligned.bam` per run
  accumulate when consumers lag (section 1). `snakemake --delete-temp-output -n`
  lists what is still being held.

## 7. Plan the one-off big jobs

The reference builds (`umbrella_hisat2_index`, `umbrella_salmon_index` at 32 GB,
`umbrella_whippet_index`, and the MAJIQ/QAPA references if enabled) run once
per reference but need the most memory. Build the reference first with its own
target while nothing else is running, then start `quant_umbrella`. That way the
per-run jobs don't compete with them for CPUs and memory.

## 8. Slow runs hold up the joins

Shards, joins and comparisons wait for every run in their group and batch, so
one very deep run can leave most of the 100 cores idle at the end of a batch.
Order the manifest, or split batches, so the deepest runs start early. Where
it's practical, keep batches of similar depth together.

## Before the big run

1. `snakemake ... -n` with the new profile, to check the `cpus` overrides are
   accepted.
2. Run one batch of 5–10 runs. Check `qstat -u $USER` core totals, real memory
   per rule, and scratch growth.
3. Adjust `ppn`/`threads`, then scale `-j` and `get_data`.

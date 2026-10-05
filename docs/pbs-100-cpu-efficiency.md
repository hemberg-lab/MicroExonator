# Configuring `quant_umbrella` for a 100-CPU PBS limit

Configuration changes for running the umbrella on a PBS cluster where one user
can hold at most 100 CPUs at a time, with a Snakemake profile (`base`) that
reads `cluster.PBS.json`. Current command:

```bash
snakemake -s MicroExonator.smk quant_umbrella --profile base \
  --use-conda --conda-frontend conda --conda-prefix $C \
  -k --rerun-incomplete -j 50 --resources get_data=8
```

The changes touch three places, and none of them edits the rules:

- the profile, `base/config.yaml`
- `cluster.PBS.json`
- the command line

The rule set uses `retries:`, so this assumes Snakemake 7 or later. Every
option below exists there.

Three values must agree for each rule:

- `threads` (from the rule, or overridden with `set-threads`)
- `ppn` in `cluster.PBS.json`
- the `cpus` resource from `set-resources`

If they drift apart, the tool, PBS and Snakemake each count a different number
of cores.

## 1. Profile: count CPUs, not jobs

`-j 50` limits Snakemake to 50 submitted jobs, but each job asks PBS for its own
`ppn` (1 to 8), so 50 jobs can mean anywhere from 50 to 400 CPUs. PBS holds the
extras in its own queue. Then:

- Snakemake no longer controls what runs first, so temporary FASTQs and BAMs pile up on scratch.
- Jobs waiting in PBS take up `-j` slots that cheap 1-CPU jobs could have used.

Giving every job a `cpus` resource lets Snakemake enforce the limit itself:

```yaml
# base/config.yaml
default-resources:
  - cpus=1
set-resources:
  - umbrella_hisat2:cpus=8
  - Round2_bowtie_to_tags:cpus=4
  - Round1_bwa_mem_to_tags:cpus=5
  - umbrella_stage_reads:cpus=2
  - umbrella_stage_single:cpus=2
  - umbrella_legacy_fastq:cpus=4
  - umbrella_featurecounts:cpus=4
  - umbrella_rmats_prep:cpus=4
  - umbrella_salmon_quant:cpus=4
  - umbrella_fastqc:cpus=2
  - umbrella_hisat2_index:cpus=8
  - umbrella_salmon_index:cpus=8
  - umbrella_rmats_post:cpus=4
  - umbrella_leafcutter:cpus=4
  # add any enabled module rules with threads > 1 (umbrella_majiq_*, umbrella_qapa_*, umbrella_dapars2_compare)
```

## 2. Profile: thread overrides

`set-threads` changes a rule's `threads` without editing the rule. Two rules
benefit:

```yaml
set-threads:
  # network-bound while downloading (the rule's own note: 24-59 % of one core);
  # at 4, get_data=8 holds 32 of the 100 CPUs mostly idle
  - umbrella_stage_reads=2
  - umbrella_stage_single=2
  # the single-threaded gzip | awk | python stages around bowtie set the pace;
  # 4 instead of 8 runs twice as many samples at once
  - Round2_bowtie_to_tags=4
```

Time one run at both values before you commit to them. Keep
`umbrella_hisat2` at 8 threads: it scales well, and `samtools sort -@` uses the
same threads.

## 3. Profile: job monitoring and submission

```yaml
cluster-status: "base/pbs_status.py"   # see below
latency-wait: 60
max-jobs-per-second: 5
max-status-checks-per-second: 1
```

- **`cluster-status`.** Without one, Snakemake doesn't notice when PBS kills a
  job for exceeding walltime or memory, and `-k` waits forever. The script gets
  the job ID and must print `running`, `success` or `failed`. It can work it out
  from `qstat -fx <jobid>` (PBS Pro) or `qstat -f <jobid>` / `tracejob <jobid>`
  (Torque), using `job_state` and `exit_status`.
- **`latency-wait`.** Prevents false "missing output" failures on a slow shared
  filesystem.
- **Submission rates.** Keep the PBS server responsive when many short jobs are
  submitted at once.

Check that the `cluster:` line passes the CPU request through, or every job runs
on one core:

- Torque:
  `qsub -l nodes={cluster.nodes}:ppn={cluster.ppn},mem={cluster.mem},walltime={cluster.walltime} -q {cluster.queue} -N {cluster.name} -o {cluster.output} -e {cluster.error}`
- PBS Pro:
  `qsub -l select=1:ncpus={cluster.ppn}:mem={cluster.mem} -l walltime={cluster.walltime} ...`

## 4. `cluster.PBS.json`

**Match the new thread counts:**

```json
"umbrella_stage_reads":  { "ppn": "2", "threads": "2" },
"umbrella_stage_single": { "ppn": "2", "threads": "2" },
"Round2_bowtie_to_tags": { "ppn": "4", "threads": "4" }
```

(Keep the other fields these entries already have.)

**Entries that request more CPUs than the rule uses:**

| Rule | ppn now | Rule threads | Set |
|---|---|---|---|
| `bowtie_genome_index` | 8 | 1 | `ppn: 1` |
| `Output` | 2 | 1 | `ppn: 1` |

**Entries that are missing.** These rules fall back to `ppn: 1` but run with more
threads, which oversubscribes the node. Add them only if you use these aligners:

| Rule | Rule threads |
|---|---|
| `total_STAR_to_Genome` | 5 |
| `mv_STAR` | 5 |
| `total_tophat_to_Genome` | 5 |
| `total_olego_to_Genome` | 10 |

**Default memory and walltime too large.** About 30 short legacy rules use
`__default__` (10 GB, 24 h), including:

- `Round1_alingment_pre_processing`
- `ME_psi_to_quant`
- `SamToBam`, `BamIndex`
- `correct_quant`, `get_PSI_sparse_quants_*`
- `detection_filter_*`
- the `*delta*` rules

PBS backfills short jobs into idle gaps sooner, and smaller memory requests fit
on more nodes. A tighter default fixes all of them at once:

```json
"__default__": {
    "mem_mb": "4000",
    "mem": "4000mb",
    "walltime": "04:00:00",
    "runtime": 14400,
    ...
}
```

Keep explicit, larger values for the rules that need them. The umbrella
entries already set their own. Among the legacy rules, give an explicit entry
to any that relied on the old default, such as `download_fastq`,
`Round1_filter` and `ME_reads`. Use real peak memory from `qstat -f` or the
logs as the guide.

## 5. Command line

```bash
snakemake -s MicroExonator.smk quant_umbrella --profile base \
  --use-conda --conda-frontend conda --conda-prefix $C \
  -k --rerun-incomplete -j 150 --resources get_data=8 cpus=96
```

- **`cpus=96`.** This is now the real limit. It sits a little under 100 to leave
  room for anything else of yours running at the same time.
- **List both resources here.** A command-line `--resources` replaces the
  profile's `resources:` entry instead of adding to it.
- **`-j 150`.** Raised so that 1-CPU jobs always find a free slot. `cpus`
  enforces the CPU limit.
- **`get_data=8`.** Still sensible. At 2 threads, 8 downloads now hold 16 CPUs
  instead of 32. Going higher mainly fills scratch with reads waiting to be
  aligned.

Before the full run, do a dry run (`-n`) to check that the profile is accepted.
Once it's running, check that `qstat -u $USER` shows about 100 CPUs in use.

## 6. Pipeline config

- **`umbrella_stage_tmpdir`.** Point it at node-local scratch if the nodes have
  it. Staging writes uncompressed reads there, and many concurrent jobs doing
  that on the shared filesystem slow each other down.

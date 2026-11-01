# Continuing an existing umbrella run with circular Whippet output

Deploy executable code by Git commits. Copy only run-specific manifests,
approval records, scheduler configuration and other data outside Git. The
helpers below are versioned under `src/`; no code archive or copied helper
script is needed on the cluster.

Stop existing drivers/jobs before renaming a run, changing its code or
archiving its outputs. Retain its config and `.snakemake` metadata. Use a
compatibility symlink for the former root while historical absolute paths
remain in the config or staging provenance. Do not clear metadata or bypass
protected outputs to force reuse.

## Compose an unpaired continuation

`src/umbrella_continuation.py` takes an existing run and separately approved
additional inputs. It verifies the old cohort against explicit approved
manifests/comparisons, preserves old processing and shard identities, retains
existing settings, and produces a new plan directory. Existing files are not
rewritten. Unexpected old inputs or an existing plan directory are refused.

```bash
"$RUN_PYTHON" src/umbrella_continuation.py \
  --root "$PWD" --out "$PWD/continuation" \
  --previous-root /previous/run/path --reference-id EXPECTED_REFERENCE_ID \
  --additional-manifest umbrella_continuation_inputs/additional-manifest.tsv \
  --additional-comparisons umbrella_continuation_inputs/additional-comparisons.yaml \
  --approved-manifest umbrella_continuation_inputs/approved-existing-manifest.tsv \
  --approved-comparisons umbrella_continuation_inputs/approved-existing-comparisons.yaml
```

This helper deliberately emits unpaired comparisons. It keeps any deferred
`design` declarations in `suspended-paired-designs.json`; it does not silently
enable another model. Only source inputs already approved for that unpaired
adaptation should be supplied. Metadata admission/preflight does not validate
experimental independence or native-tool sequencing QC.

The plan includes baseline and combined manifests, comparisons, prior and
proposed configs, a path-only migration config, explicit Whippet backfill
targets and a recoverable archive plan. Old rows keep the extra
`whippet_read_source: backfill` metadata on later resumes. It changes only
Whippet's private read staging, not ME/Salmon processing hashes. New inputs
use ordinary shared staging.

Review and, if applicable, apply `src/staging_provenance.py` using the plan's
`prior.config.yaml` and `migration-proposed.yaml`. That migration accepts only
manifest/comparison path changes and saves exact metadata backups. Do not mix
new Whippet flags into the path-only migration config.

## Archive and backfill Whippet

```bash
"$RUN_PYTHON" src/whippet_backfill.py preview --plan continuation
"$RUN_PYTHON" src/whippet_backfill.py archive --plan continuation
"$RUN_PYTHON" src/umbrella_launch.py backfill \
  --plan continuation --profile umbrella_continuation_inputs/profile
```

The archive action moves only the planned old Whippet directories and linear
Whippet shards to a recoverable local backup. It records the unchanged file
hashes and timestamps of ME/Salmon results and original shared validation
markers. It refuses repeated or interrupted archival until the receipt is
reviewed. Only after inspecting the real production dry run and resources,
add `--submit` to launch the backfill.

The launcher uses the current interpreter's sibling Snakemake executable,
requires Python 3.11 and Snakemake 7, places explicit targets before options,
and defaults to a dry run. It defaults to 12 scheduler jobs and acquisition
concurrency 2; neither is a disk-space cap. The public `quant_whippet` target
also requests MicroExonator work and is not used for this isolated phase.

## Verify and continue

```bash
"$RUN_PYTHON" src/whippet_backfill.py verify --plan continuation
cp continuation/config.continue.yaml config.yaml
"$RUN_PYTHON" src/umbrella_launch.py continue \
  --plan continuation --profile umbrella_continuation_inputs/profile
```

Verification requires certified circular core outputs, unchanged non-BS
linear shard contents, and unchanged ME/Salmon/shared-validation hashes and
timestamps. A continuation requires the verified manifest/comparison hashes
to remain unchanged. Review the dry-run log with
`src/umbrella_resume_check.py LOG --plan continuation --require-backfill`.
The `continue --submit` launcher repeats that dry run and rejects old-library
ME/Salmon producer jobs before submitting. Review other jobs, quota and
scheduler resources too; this is not a full cluster audit.

`src/slurm_status.py` is the versioned Snakemake 7 Slurm status callback.
Point the copied scheduler profile at its absolute path in the Git checkout.
No code should be copied into the profile directory.

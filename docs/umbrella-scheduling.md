# Disk-saving umbrella priorities (Snakemake 7)

Opt in with the following key in the active configuration:

```yaml
umbrella_disk_saving_priority: true
```

Or append `--config umbrella_disk_saving_priority=true` to the normal Snakemake
invocation. Keep the same target, manifest, comparisons and enabled tools.
This works for both MicroExonator + Whippet selection and the full umbrella.
Absent/false preserves existing rule priorities; non-umbrella workflows are
unaffected. Invalid boolean values fail early.

The policy assigns priorities after all rules are loaded, without changing
producer commands, parameters, inputs, outputs, threads or temporary flags:

| Priority | Work |
|---|---|
| 110 | Validation, legacy bridge, read-length extraction, ME preprocessing and ME read extraction |
| 100 | Parallel FASTQ readers/conversion and expensive per-library BAM consumers |
| 90 | Whippet and other direct BAM readers |
| 80 | ME filtering/coverage steps consuming temporary ME intermediates |
| 0 | New native FASTQ acquisition |

The exact rule map is in `src/disk_saving_priority.py`. Other rules retain their
original priorities. The map is static, not a last-consumer detector; no tool
is enabled by assigning it a priority. References and group comparisons do not
gain new dependencies. Lane merging remains within the existing library job.

Snakemake still deletes temporary files only after all required consumers
finish. The preference is not strict depth-first execution or a disk quota.
Jobs must be ready and fit resource limits; cluster queues control the start
order after submission. The acquisition resource limits running/submitted
acquisition jobs, not accumulated libraries or retained bytes. A failed
consumer can leave required reads on disk. Do not manually delete intermediates
that unfinished consumers need.

## Resuming a stopped direct run

Wait until no old driver or submitted jobs remain, make enough free disk space
available, and update code without overwriting local run configuration. Keep
the original tool selection and all normal provenance checks. Run the original
command with the additional scheduling option and `--dry-run`, then inspect
unexpected old-library staging, protected-output errors or reference rebuilds
before removing `--dry-run`. The scheduling option alone needs no provenance
migration; a manifest-path migration is a separate issue documented in
`umbrella-staging-provenance.md`.

Changing this setting alone must not rerun completed outputs in a direct
Snakemake workspace. It does not remove incomplete or disk-full job records;
retain the usual `--rerun-incomplete`. No automatic unlocking, output deletion,
permission changes or metadata cleanup is part of this feature.

The experimental execution controller separately freezes workflow/config
content. Existing controller packets cannot silently adopt changed code or
settings: follow its explicit revision contract. Do not switch an existing
direct run to that controller to use this scheduling option.

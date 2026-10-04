# Guarded two-pass umbrella executions (experimental)

This opt-in launcher separates MicroExonator + Whippet discovery from later
supporting analyses. Existing direct Snakemake launches are unchanged unless
`umbrella_execution.enabled` is true. Do not enroll an existing production run:
it has no controller completion ledger. Full native-tool and PBS transition
validation is still required before production use.

## Configuration

Keep your master manifest, comparisons and reference configuration unchanged.
Add this block to a **separate** configuration file:

```yaml
umbrella_execution:
  enabled: true
  packet_id: pilot-v1
  execution_id: discovery-01
  mode: discovery
```

Discovery always requests MicroExonator and ordinary Whippet. For MegaSearch,
the generated configuration uses standalone `delta_method: microexonator` and
disables `umbrella_optional.microexonator_whippet_delta`. The optional legacy
hybrid implementation remains available outside these execution profiles.
MicroExonator's own PSI correction is unaffected.

After discovery completes, change only the execution block in that separate
file, keeping the same packet ID:

```yaml
umbrella_execution:
  enabled: true
  packet_id: pilot-v1
  execution_id: complement-01
  mode: complement
  comparisons: [comparison-of-interest]
  tools: [hisat2, salmon, rmats, coverage, suppa2, majiq, dapars2, qapa, multiqc]
```

Complement excludes MicroExonator, Whippet and LeafCutter. Tools can be narrowed
to an admitted subset. MAJIQ/DaPars2/QAPA use `ingest`, not a cache-only launch.
Omitting selection uses all included master libraries and comparisons. Selecting
only comparisons selects their participating libraries. Optional `libraries`
accepts processing-library IDs or member run IDs and retains every lane in each
selected library. Explicit library selection must still satisfy both sides of
each selected comparison. Use `comparisons: []` with explicit libraries for
quantification without differential comparisons. No group or batch is relabelled.

## Launch and resume

From the checkout, with Python 3, PyYAML and Snakemake 7 available:

```bash
python3 src/umbrella_execution.py --configfile config.execution.yaml --dry-run -- --profile base --use-conda --conda-frontend conda --conda-prefix "$C" -k --rerun-incomplete -j 50 --resources get_data=8
```

Remove the launcher's `--dry-run` to execute. A dry run creates an attempt record
but cannot mark tools complete. Repeat the same command, execution ID and
selection after an interruption. Changed selection needs a new execution ID;
changed immutable inputs or producer configuration needs a new explicit packet
revision. Completed discovery cannot be regenerated under another execution ID.

The launcher accepts a bounded set of scheduler/resource/conda options. Arbitrary
targets, force/rebuild options, alternate work directories, config overrides,
and passthrough dry-run flags are rejected. Profiles are resolved to absolute
paths and audited; an inherited `SNAKEMAKE_PROFILE` is ignored. The root
`cluster.PBS.json`, if present, is copied into the execution workspace and hashed.
No PBS integration test has been performed.

## Outputs and guarantees

Outputs live under `umbrella/executions/<packet_id>/<execution_id>/work/`, with
their own generated manifest/config, Snakemake metadata and scientific results.
The master files remain untouched. A complement may reacquire transient reads.
Its planned producers are checked before launch, and cannot write through a
reused reference symlink or schedule discovery producers.

`packet.json` records selections, attempts, requested/completed tools, identities,
Snakemake version and retained file hashes. Completed products and reports are
made read-only and verified before reuse. Reference/index contents are fully
hashed, including Whippet's exon sidecar and seeded assets; hashing large inputs
has an I/O cost. No genome reindexing is performed by the controller itself.

Each execution keeps its own synthesis/QC. `cumulative_inventory.json` and
`cumulative_reports.json` collect paths, hashes and cohort/execution scope across
completed passes. These are provenance-aware collections, **not a merged
cross-tool biological consensus table**. Completion means validated output
inventory, not evidence that differential inference succeeded; inspect each
tool's status and results. Actual tool/conda-package version capture is not yet
complete. Keep this workflow experimental until native transitions, scheduler
behavior and any required cross-pass biological synthesis are validated.

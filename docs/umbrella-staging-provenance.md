# Reusing libraries when the manifest moves

Staging tracks the per-library processing hash and scratch/restaging options,
not the location of the master manifest. A new comparison or manifest path must
not reacquire deleted temporary reads for unchanged completed libraries.
Protected validation markers and immutable shard checks remain in place.
Changes to sources, layout, library membership or reference still invalidate
processing; this fix does not make those changes reusable.

## Existing Snakemake 7 workspaces

Older staging metadata records the manifest path as a parameter. Updating the
rule alone therefore does not remove the first false rerun. After every driver
and submitted job has finished, preview the explicit migration:

```bash
python3 src/staging_provenance.py --prior-config previous/config.yaml --configfile correction/config.yaml
```

Then, only if the reviewed paths/settings are correct:

```bash
python3 src/staging_provenance.py --prior-config previous/config.yaml --configfile correction/config.yaml --apply
```

It accepts only manifest/comparison path changes in the config, checks all old
processing/shard identities, refuses workspace locks and changed settings,
and validates each historical manifest named in the recorded staging command.
That supports both original and later intake manifests in one workspace.
A missing historical manifest, mismatched hash, incomplete staging record or
missing completed validation marker stops the migration.

Only recognised `umbrella_stage_reads` / `umbrella_stage_single` records change.
Exact originals and a receipt are saved under
`umbrella/provenance-migrations/staging-manifest-path-v1-*/` before updates.
It removes just the obsolete path parameter and updates the code fingerprint
from the verified legacy template to the current template. Future code changes
remain detectable. Other parameters, inputs, environment and
timestamps remain. It never deletes results, changes permissions, enables
restaging, removes protection or cleans all metadata.

Run a normal dry run afterward with all provenance triggers enabled. Only
new controls should need staging/per-library analysis. Stop for unexpected old
samples or reference jobs. Do not use blanket mtime-only triggers, metadata
cleanup or touching results to force reuse.

Migration is not automatic during workflow parsing; repeating it skips already
migrated records. If interrupted, keep its backup receipt and inspect before
continuing. Exact backup record files can be restored to `.snakemake/metadata/`
while no driver/jobs are active; never copy `receipt.json` into metadata.

## Regression coverage

Tiny native Snakemake tests use four-base local FASTQs and the real staging /
validation rules. No genome, index or network is used. Tests cover relocation
after temporary-read deletion, new-library-only scheduling, source-change
refusal and real legacy paired-end provenance migration. Unit tests cover
preview, backup, idempotence, mixed historical paths, locks and hash guards.
Actual native-tool/PBS execution remains subject to the cluster dry run.

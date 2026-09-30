"""Keep existing umbrella shards immutable when the manifest is edited.

Shards are filed by {reference_id}/{project_id}/{group}/{batch_id}, and the
per-run files they were built from are temporary, so a shard can only be
rebuilt by downloading its reads again. Two checks protect them:

check_existing_shards()
    Compares every shard already on disk with the manifest. A run that moved
    to another group or batch, or a shard whose runs changed (added, removed or
    excluded), is refused with a message saying what to do instead: select runs
    on a metadata column in the comparisons file, put new runs in a new
    batch_id, or leave a run out of a comparison with `exclude`.

migrate_identity_params()
    One-time step for workflows started before the run identity was narrowed
    to the processing fields (umbrella_manifest.RUN_IDENTITY). Snakemake
    recorded the old whole-row hashes as params of the staging and shard
    rules; with the new hashes those would look changed and rerun, which means
    downloading everything again. This drops the recorded params of those
    rules' existing outputs (start and end times, code and inputs are kept),
    once, and leaves a marker so later runs skip it.
"""

import glob
import gzip
import json
import os
from pathlib import Path


# shard file (relative pattern), number of fixed columns, columns per run
SHARD_FILES = (
    ("splicing/{}.microexonator.tsv.gz", 1, 5),
    ("junctions/{}.junctions.tsv.gz", 8, 1),
    ("genes/{}.featurecounts.tsv.gz", 2, 1),
    ("genes/{}.salmon_counts.tsv.gz", 1, 1),
    ("splicing/{}.whippet.tsv.gz", 5, 6),
)

# rules whose params carried a manifest hash before the identity was narrowed
IDENTITY_RULES = frozenset((
    "umbrella_stage_reads", "umbrella_stage_single", "umbrella_validate_reads",
    "umbrella_junction_shard", "umbrella_capture_shard", "umbrella_featurecounts_shard",
    "umbrella_qc_shard", "umbrella_rmats_inventory", "umbrella_coverage_shard",
    "umbrella_salmon_shard", "umbrella_microexonator_shard", "umbrella_whippet_shard",
))
MIGRATION_MARKER = "umbrella/.run_identity_v2"


def shard_runs(path, fixed, width):
    """Run IDs recorded in a shard's header."""
    with gzip.open(path, "rt") as stream:
        header = stream.readline().rstrip("\n").split("\t")
    columns = header[fixed:]
    if width == 1:
        return set(columns)
    return {column.rsplit(".", 1)[0] for column in columns[::width]}


def existing_shards(reference_ids, root="."):
    """{(reference_id, project_id, group, batch_id): set of run IDs} for shards on disk."""
    root = Path(root)
    found = {}
    for pattern, fixed, width in SHARD_FILES:
        prefix, suffix = pattern.split("{}")
        for path in glob.glob(str(root / pattern.format("*/*/*/*"))):
            relative = Path(path).relative_to(root).as_posix()
            parts = relative[len(prefix):-len(suffix)].split("/")
            if len(parts) != 4 or parts[0] not in reference_ids:
                continue
            key = tuple(parts)
            if key not in found:
                found[key] = shard_runs(path, fixed, width)
    return found


def check_existing_shards(manifest, root="."):
    """Messages for every way the manifest contradicts shards already written."""
    included = manifest.included_runs()
    reference_ids = {run.reference_id for run in included}
    shards = existing_shards(reference_ids, root)
    if not shards:
        return []
    stored_in = {}
    for key, runs in shards.items():
        for run in runs:
            stored_in.setdefault(run, set()).add(key)
    problems = []
    for run in included:
        here = (run.reference_id, run.project_id, run.group, run.batch_id)
        elsewhere = sorted(key for key in stored_in.get(run.run_id, ()) if key != here
                           and key[0] == run.reference_id)
        if elsewhere:
            problems.append(
                "run {} is stored in shard {} but the manifest now puts it in {}/{}. "
                "Keep its original group and batch_id; to compare it under a new condition, "
                "add a metadata column and select on it in the comparisons file".format(
                    run.run_id, "/".join(elsewhere[0]), run.group, run.batch_id))
    moved = {problem.split()[1] for problem in problems}
    for key, recorded in sorted(shards.items()):
        current = {run.run_id for run in included
                   if (run.reference_id, run.project_id, run.group, run.batch_id) == key}
        added = sorted(current - recorded - moved)
        removed = sorted(recorded - current - moved)
        if added:
            problems.append(
                "shard {} already exists without run(s) {}. Shards are immutable: "
                "give new runs a new batch_id".format("/".join(key), ", ".join(added)))
        if removed:
            problems.append(
                "shard {} holds run(s) {} that the manifest no longer includes. Keep them "
                "included and leave them out of a comparison with `exclude` instead".format(
                    "/".join(key), ", ".join(removed)))
    return problems


def migrate_identity_params(metadata_dir=".snakemake/metadata", marker=MIGRATION_MARKER):
    """Drop recorded params of IDENTITY_RULES outputs, once. Returns how many records changed."""
    marker = Path(marker)
    if marker.exists():
        return 0
    changed = 0
    for directory, _, files in os.walk(metadata_dir):
        for name in files:
            path = os.path.join(directory, name)
            try:
                with open(path) as stream:
                    record = json.load(stream)
            except (OSError, ValueError):
                continue
            if record.get("rule") not in IDENTITY_RULES or "params" not in record:
                continue
            del record["params"]
            temporary = "{}.{}.tmp".format(path, os.getpid())
            with open(temporary, "w") as stream:
                json.dump(record, stream)
            os.replace(temporary, path)
            changed += 1
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text("run identity v2 ({}): dropped recorded params of {} outputs\n".format(
        ", ".join(sorted(IDENTITY_RULES)), changed))
    return changed

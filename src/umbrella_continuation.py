"""Compose an in-place, unpaired continuation without modifying the live run."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import yaml



def read_tsv(path):
    with Path(path).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path, columns, rows):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def local_path(root, value, previous_root, destination_root):
    for prefix in filter(None, (str(previous_root or ""), str(destination_root))):
        if value.startswith(prefix + "/"):
            return root / value[len(prefix) + 1:]
    path = Path(value)
    return path if path.is_absolute() else root / path


def compose(existing_root, out, *, additional_manifest, additional_comparisons,
            approved_manifest, approved_comparisons, reference_id, previous_root=None, destination_root=None):
    existing_root, out = Path(existing_root), Path(out)
    if out.exists():
        raise FileExistsError(out)
    destination_root = Path(destination_root or existing_root.resolve())
    if __package__:
        from src.umbrella_manifest import load_umbrella_manifest, REQUIRED
        from src.comparison_preflight import load_comparisons, preflight
        from src.umbrella_identity import check_existing_shards
        from src.whippet_flags import quant_flags
    else:
        from umbrella_manifest import load_umbrella_manifest, REQUIRED
        from comparison_preflight import load_comparisons, preflight
        from umbrella_identity import check_existing_shards
        from whippet_flags import quant_flags
    prior_path = existing_root / "config.yaml"
    original = prior_path.read_bytes()
    prior = yaml.safe_load(original)
    if prior.get("umbrella_allow_restage") not in (None, False, "F", "false"):
        raise ValueError("Inspect prior restaging authorization before continuing")
    old_manifest_path = local_path(existing_root, prior["umbrella_manifest"], previous_root, destination_root)
    old_comp_path = local_path(existing_root, prior["umbrella_comparisons"], previous_root, destination_root)
    old_columns, old_rows = read_tsv(old_manifest_path)
    new_columns, new_rows = read_tsv(additional_manifest)
    baseline_manifest = Path(approved_manifest)
    baseline_comps = Path(approved_comparisons)
    approved_old = {r["run_id"]: r for r in read_tsv(baseline_manifest)[1]}
    for row in old_rows:
        approved = approved_old.get(row["run_id"])
        if approved is None:
            raise ValueError("Existing run is outside the approved baseline cohort: " + row["run_id"])
        for field in REQUIRED:
            if field == "reference_id":
                if row[field] not in ("auto", reference_id):
                    raise ValueError("Existing reference differs")
            elif row[field] != approved[field]:
                raise ValueError(f"Existing approved metadata differs: {row['run_id']} {field}")
    old_comps = yaml.safe_load(old_comp_path.read_text())["comparisons"]
    allowed = {c["comparison_id"]: c for c in yaml.safe_load(baseline_comps.read_text())["comparisons"]}
    for c in old_comps:
        if c != allowed.get(c["comparison_id"]):
            raise ValueError("Existing comparison is not in the confirmed approved baseline definitions: " + c["comparison_id"])
    # This extra column affects only Whippet's read source. Shared staging,
    # ME and Salmon processing hashes remain unchanged.
    if "whippet_read_source" not in old_columns:
        old_columns.append("whippet_read_source")
    for row in old_rows:
        row["whippet_read_source"] = "backfill"
    for row in new_rows:
        row["whippet_read_source"] = "shared"
    columns = old_columns + [c for c in new_columns if c not in old_columns]
    old_ids = {r["run_id"] for r in old_rows}
    if old_ids & {r["run_id"] for r in new_rows}:
        raise ValueError("Overlapping intake runs need an explicit reuse review")
    combined = [{c: row.get(c, "") for c in columns} for row in old_rows + new_rows]
    new_comps = yaml.safe_load(Path(additional_comparisons).read_text())["comparisons"]
    suspended = {c["comparison_id"]: c.pop("design") for c in new_comps if "design" in c}
    comparisons = old_comps + new_comps
    if len({c["comparison_id"] for c in comparisons}) != len(comparisons):
        raise ValueError("Duplicate comparison ID")
    out.mkdir(parents=True)
    write_tsv(out / "manifest.tsv", columns, combined)
    # Keep a local baseline copy for the narrow backfill target and audit.
    write_tsv(out / "baseline-manifest.tsv", old_columns, old_rows)
    (out / "comparisons.yaml").write_text(yaml.safe_dump({"comparisons": comparisons}, sort_keys=False))
    (out / "baseline-comparisons.yaml").write_text(yaml.safe_dump({"comparisons": old_comps}, sort_keys=False))
    (out / "suspended-paired-designs.json").write_text(json.dumps(suspended, indent=2) + "\n")
    manifest = load_umbrella_manifest(out / "manifest.tsv").with_reference_id(reference_id)
    baseline = load_umbrella_manifest(out / "baseline-manifest.tsv").with_reference_id(reference_id)
    for run in baseline.included_runs():
        if baseline.run_sha256(run.run_id) != manifest.run_sha256(run.run_id):
            raise ValueError("Old processing identity changed: " + run.run_id)
        args = run.reference_id, run.project_id, run.group, run.batch_id
        if baseline.shard_sha256(*args) != manifest.shard_sha256(*args):
            raise ValueError("Old shard changed: " + run.run_id)
    problems = check_existing_shards(manifest, existing_root)
    if problems:
        raise ValueError("\n".join(problems))
    preflights = [preflight(manifest, c) for c in load_comparisons(out / "comparisons.yaml")]
    if not all(p["inference_supported"] for p in preflights):
        raise ValueError("A comparison failed native metadata preflight")
    (out / "preflights.json").write_text(json.dumps(preflights, indent=2) + "\n")
    try:
        relative_plan = out.resolve().relative_to(existing_root.resolve())
    except ValueError:
        relative_plan = Path(out.name)
    remote = str(destination_root / relative_plan)
    cfg = dict(prior, working_directory=str(destination_root) + "/", whippet_circ_source="core")
    cfg["whippet_flags"] = quant_flags(prior.get("whippet_flags", ""))
    cfg.update(umbrella_manifest=remote + "/manifest.tsv", umbrella_comparisons=remote + "/comparisons.yaml")
    (out / "prior.config.yaml").write_bytes(original)
    (out / "config.continue.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    backfill = dict(cfg, umbrella_manifest=remote + "/baseline-manifest.tsv",
                    umbrella_comparisons=remote + "/baseline-comparisons.yaml")
    (out / "config.backfill.yaml").write_text(yaml.safe_dump(backfill, sort_keys=False))
    # Migration accepts only path changes; do not mix flags/settings into it.
    migration = dict(prior, umbrella_manifest=remote + "/manifest.tsv", umbrella_comparisons=remote + "/comparisons.yaml")
    (out / "migration-proposed.yaml").write_text(yaml.safe_dump(migration, sort_keys=False))
    targets, archives = [], []
    for r in baseline.included_runs():
        targets += [r.work_dir + "/whippet/quant.provenance.json"]
        folder = r.work_dir + "/whippet"
        if (existing_root / folder).exists():
            archives.append(folder)
    groups = sorted({(r.reference_id, r.project_id, r.group, r.batch_id) for r in baseline.included_runs()})
    for group in groups:
        prefix = "splicing/" + "/".join(group)
        for suffix in (".whippet.tsv.gz", ".whippet_checksums.json"):
            path = prefix + suffix
            targets.append(path)
            if (existing_root / path).exists():
                archives.append(path)
    (out / "backfill-targets.txt").write_text("\n".join(targets) + "\n")
    (out / "backfill-archive.json").write_text(json.dumps(archives, indent=2) + "\n")
    summary = dict(reference_id=reference_id, old_raw_runs=len(old_rows), new_raw_runs=len(new_rows), raw_runs=len(combined),
                   processing_inputs=len(manifest.included_runs()), comparisons=len(comparisons),
                   old_comparisons=len(old_comps), new_comparisons=len(new_comps), paired_mode=False,
                   old_processing_identities_unchanged=len(baseline.included_runs()),
                   source_config_sha256=hashlib.sha256(original).hexdigest(),
                   source_manifest_sha256=hashlib.sha256(old_manifest_path.read_bytes()).hexdigest(),
                   source_comparisons_sha256=hashlib.sha256(old_comp_path.read_bytes()).hexdigest(),
                   max_read_len=cfg.get("max_read_len"), min_number_files_detected=cfg.get("min_number_files_detected"))
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--out", type=Path)
    parser.add_argument("--additional-manifest", type=Path, required=True)
    parser.add_argument("--additional-comparisons", type=Path, required=True)
    parser.add_argument("--approved-manifest", type=Path, required=True)
    parser.add_argument("--approved-comparisons", type=Path, required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--previous-root")
    parser.add_argument("--destination-root", type=Path)
    args = parser.parse_args()
    print(json.dumps(compose(args.root, args.out or args.root / "continuation",
                            additional_manifest=args.additional_manifest, additional_comparisons=args.additional_comparisons,
                            approved_manifest=args.approved_manifest, approved_comparisons=args.approved_comparisons,
                            reference_id=args.reference_id, previous_root=args.previous_root,
                            destination_root=args.destination_root), indent=2))

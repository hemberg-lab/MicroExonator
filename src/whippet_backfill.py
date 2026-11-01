"""Back up only old Whippet outputs, then verify circular backfill and reused results."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sys


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def no_locks(root):
    locks = root / ".snakemake/locks"
    if locks.exists() and any(locks.iterdir()):
        raise ValueError("Wait until all Snakemake drivers/jobs have stopped; workspace has locks")


def archive(root, plan, apply=False):
    no_locks(root)
    receipt = plan / "backfill-receipt.json"
    if receipt.exists():
        raise ValueError("Backfill receipt already exists; use verify or inspect the saved recovery plan")
    paths = json.loads((plan / "backfill-archive.json").read_text())
    witnesses = {}
    for prefix, pattern in (("splicing", "*.microexonator*"), ("genes", "*.salmon*")):
        for path in (root / prefix).rglob(pattern):
            if path.is_file():
                witnesses[str(path.relative_to(root))] = {"sha256": sha(path), "mtime_ns": path.stat().st_mtime_ns}
    for path in (root / "umbrella/work").rglob("reads.valid"):
        witnesses[str(path.relative_to(root))] = {"sha256": sha(path), "mtime_ns": path.stat().st_mtime_ns}
    summary = {"archive_paths": paths, "preserved_result_files": len(witnesses)}
    if not apply:
        return summary
    backup = root / "continuation-backups" / ("whippet-pre-circ-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    backup.mkdir(parents=True)
    state = dict(summary, backup=str(backup), witnesses=witnesses, moved=[], complete=False)
    receipt.write_text(json.dumps(state, indent=2) + "\n")
    for relative in paths:
        source = root / relative
        if not source.exists():
            raise ValueError("Planned Whippet output disappeared: " + str(source))
        target = backup / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise ValueError("Refusing to overwrite backup " + str(target))
        source.rename(target)  # same filesystem; originals remain fully recoverable
        state["moved"].append(relative)
        receipt.write_text(json.dumps(state, indent=2) + "\n")
    state["complete"] = True
    receipt.write_text(json.dumps(state, indent=2) + "\n")
    return state


def verify(root, plan):
    no_locks(root)
    sys.path.insert(0, str(root))
    from src.umbrella_manifest import load_umbrella_manifest
    from src.whippet_core import validate_core
    state = json.loads((plan / "backfill-receipt.json").read_text())
    if not state["complete"]:
        raise ValueError("Whippet backup was interrupted; inspect recovery receipt")
    for relative, expected in state["witnesses"].items():
        if sha(root / relative) != expected["sha256"] or (root / relative).stat().st_mtime_ns != expected["mtime_ns"]:
            raise ValueError("An existing MicroExonator/Salmon/validation file changed: " + relative)
    baseline = load_umbrella_manifest(plan / "baseline-manifest.tsv").with_reference_id(json.loads((plan / "summary.json").read_text())["reference_id"])
    for run in baseline.included_runs():
        validate_core(root / run.work_dir / "whippet", run.run_id,
                      {"reference_id": run.reference_id, "processing_sha256": baseline.run_sha256(run.run_id)})
    checked = 0
    for relative in state["moved"]:
        if relative.endswith(".whippet.tsv.gz"):
            with gzip.open(Path(state["backup"]) / relative, "rb") as old, gzip.open(root / relative, "rb") as new:
                while True:
                    a, b = old.read(1024 * 1024), new.read(1024 * 1024)
                    if a != b:
                        raise ValueError("Linear Whippet shard changed after --circ backfill: " + relative)
                    if not a:
                        break
            checked += 1
    result = dict(verified=True, circular_inputs=len(baseline.included_runs()), unchanged_linear_shards=checked,
                  unchanged_me_salmon_and_validation_files=len(state["witnesses"]),
                  continuation_manifest_sha256=sha(plan / "manifest.tsv"),
                  continuation_comparisons_sha256=sha(plan / "comparisons.yaml"))
    (plan / "backfill-verified.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("preview", "archive", "verify"))
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--plan", type=Path, default=Path("continuation"))
    args = parser.parse_args()
    plan = args.plan if args.plan.is_absolute() else args.root / args.plan
    result = verify(args.root, plan) if args.action == "verify" else archive(args.root, plan, apply=args.action == "archive")
    print(json.dumps(result, indent=2))

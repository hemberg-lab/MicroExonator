"""Reject old-library ME/Salmon producer reruns in a continuation dry-run log."""
import argparse
import json
from pathlib import Path
import hashlib
import re
import sys


def audit(log, old_ids):
    failures = []
    consumers = {"umbrella_salmon_quant", "bowtie_to_genome", "Round2_bowtie_to_tags", "ME_reads",
                 "ME_SJ_coverage", "coverage_to_PSI_report", "correct_quant", "get_PSI_sparse_quants_se"}
    chunks = re.split(r"(?m)^(?:local)?rule ", log)
    for chunk in chunks[1:]:
        rule = chunk.split(":", 1)[0]
        if rule not in consumers:
            continue
        for line in chunk.splitlines():
            if line.strip().startswith("wildcards:"):
                ids = re.findall(r"(?:sample|run_id)=([^,\s]+)", line)
                failures += [(rule, run) for run in ids if run in old_ids]
    return failures


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--plan", type=Path, default=Path("continuation"))
    parser.add_argument("--require-backfill", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(Path.cwd()))
    from src.umbrella_manifest import load_umbrella_manifest
    if args.require_backfill:
        verified = json.loads((args.plan / "backfill-verified.json").read_text())
        for name, key in (("manifest.tsv", "continuation_manifest_sha256"),
                          ("comparisons.yaml", "continuation_comparisons_sha256")):
            if hashlib.sha256((args.plan / name).read_bytes()).hexdigest() != verified[key]:
                raise SystemExit("STOP: verified continuation selection changed: " + name)
    old = {r.run_id for r in load_umbrella_manifest(args.plan / "baseline-manifest.tsv").included_runs()}
    problems = audit(args.log.read_text(), old)
    if problems:
        raise SystemExit("STOP: unexpected old ME/Salmon producer jobs: " + repr(problems))
    print("No old per-library MicroExonator/Salmon producer jobs found. Inspect other jobs and resources before submission.")

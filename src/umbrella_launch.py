"""Dry-run-first launcher for a prepared circular backfill or continuation.

Code is deployed by Git. The plan and scheduler profile are local run data.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

if __package__:
    from src.umbrella_resume_check import audit
    from src.umbrella_manifest import load_umbrella_manifest
else:
    from umbrella_resume_check import audit
    from umbrella_manifest import load_umbrella_manifest


def build_command(root, plan, profile, phase, submit=False, jobs=12, acquisitions=2, snakemake=None):
    root, plan = Path(root), Path(plan)
    if jobs < 1 or acquisitions < 1:
        raise ValueError("jobs and acquisitions must be positive")
    if phase == "backfill":
        receipt = json.loads((plan / "backfill-receipt.json").read_text())
        if not receipt.get("complete"):
            raise ValueError("Backfill archive was interrupted; inspect the recovery receipt")
        targets = (plan / "backfill-targets.txt").read_text().splitlines()
        if not targets or any(not t.startswith(("umbrella/work/", "splicing/")) for t in targets):
            raise ValueError("Unexpected backfill targets")
    elif phase == "continue":
        verified = json.loads((plan / "backfill-verified.json").read_text())
        for name, key in (("manifest.tsv", "continuation_manifest_sha256"),
                          ("comparisons.yaml", "continuation_comparisons_sha256")):
            if hashlib.sha256((plan / name).read_bytes()).hexdigest() != verified[key]:
                raise ValueError("Verified continuation selection changed: " + name)
        targets = ["quant_umbrella"]
    else:
        raise ValueError("phase must be backfill or continue")
    executable = str(snakemake or Path(sys.executable).with_name("snakemake"))
    config = plan / ("config.backfill.yaml" if phase == "backfill" else "config.continue.yaml")
    return [executable, "-s", "MicroExonator.smk", *targets, "--configfile", str(config),
            "--profile", str(profile), "-j", str(jobs), "--resources", "get_data=" + str(acquisitions),
            "--rerun-incomplete"] + ([] if submit else ["-n"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("backfill", "continue"))
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--plan", type=Path, default=Path("continuation"))
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--acquisitions", type=int, default=2)
    parser.add_argument("--snakemake", type=Path)
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 11):
        raise SystemExit("Use the verified Snakemake 7 / Python 3.11 environment")
    root = args.root.resolve()
    plan = args.plan.resolve() if args.plan.is_absolute() else root / args.plan
    profile = args.profile.resolve() if args.profile.is_absolute() else root / args.profile
    command = build_command(root, plan, profile, args.phase, args.submit,
                            args.jobs, args.acquisitions, args.snakemake)
    version = subprocess.check_output([command[0], "--version"], text=True).strip()
    if not version.startswith("7."):
        raise SystemExit("Snakemake 7 is required")
    if args.submit and args.phase == "continue":
        preview = subprocess.run(build_command(root, plan, profile, args.phase, False,
                                              args.jobs, args.acquisitions, args.snakemake),
                                 cwd=root, capture_output=True, text=True)
        text = preview.stdout + preview.stderr
        path = plan / "continue-submit-preview.log"
        path.write_text(text)
        if preview.returncode:
            raise SystemExit("Continuation dry run failed; inspect " + str(path))
        old = {r.run_id for r in load_umbrella_manifest(plan / "baseline-manifest.tsv").included_runs()}
        failures = audit(text, old)
        if failures:
            raise SystemExit("Unexpected old ME/Salmon producer jobs: " + repr(failures))
    os.chdir(root)
    os.execv(command[0], command)


if __name__ == "__main__":
    main()

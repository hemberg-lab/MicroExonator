"""Check a finished umbrella smoke run (see run_umbrella_smoke.py).

  - every kept shard matches its checksum guard and the run's reference_id
  - the junction shards recover the skipping junction of at least 80% of the
    simulated events (both flanking exons are known from truth/events.tsv)
  - preflights: both comparisons support inference; rMATS is refused for the
    mixed-layout A_vs_B and allowed for A_vs_B_single_end
  - each comparison has a synthesis report, and rMATS shows as unsupported in
    the mixed-layout one
"""

import argparse
import csv
import gzip
import json
import pathlib
import sys


FIXTURE = pathlib.Path(__file__).resolve().parent
REPOSITORY = FIXTURE.parents[2]
sys.path.insert(0, str(REPOSITORY))

from src.join_umbrella_shards import verify_shard  # noqa: E402


def skipping_junction(event):
    """1-based closed intron between the two flanking exons of a truth event."""
    exons = []
    for field in ("upstream_exon", "downstream_exon"):
        chrom, span, strand = event[field].rsplit(":", 2)
        start, end = (int(value) for value in span.split("-"))
        exons.append((start, end))
    exons.sort()
    return (chrom, exons[0][1] + 1, exons[1][0] - 1, strand)


def observed_junctions(workdir):
    found = set()
    for path in sorted(pathlib.Path(workdir, "junctions").rglob("*.junctions.tsv.gz")):
        with gzip.open(path, "rt") as stream:
            header = next(stream).rstrip("\n").split("\t")
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if any(int(value) for value in fields[8:]):
                    found.add((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return found


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True, type=pathlib.Path)
    parser.add_argument("--min-recovery", type=float, default=0.8)
    args = parser.parse_args(argv)
    workdir = args.workdir.resolve()
    with open(workdir / "config.yaml") as stream:
        config = json.load(stream)
    with open(config["umbrella_manifest"]) as stream:
        reference_id = next(csv.DictReader(stream, delimiter="\t"))["reference_id"]
    failures = []

    shards = [path for kind in ("junctions", "genes", "splicing")
              for path in (workdir / kind).rglob("*.tsv.gz")
              if not path.name.endswith(".capture.tsv.gz")]
    for path in shards:
        try:
            verify_shard(path, reference_id)
        except ValueError as error:
            failures.append(str(error))
    if not shards:
        failures.append("no shards found")

    with open(FIXTURE / "truth" / "events.tsv") as stream:
        events = list(csv.DictReader(stream, delimiter="\t"))
    found = observed_junctions(workdir)
    recovered = sum(skipping_junction(event) in found for event in events)
    if recovered < args.min_recovery * len(events):
        failures.append("skipping junctions recovered for {}/{} events".format(recovered, len(events)))

    root = workdir / "comparisons" / reference_id / "mm10sim"
    expectations = {"A_vs_B": False, "A_vs_B_single_end": True}
    for comparison, rmats_supported in expectations.items():
        preflight = json.loads((root / comparison / "preflight.json").read_text())
        if not preflight["inference_supported"]:
            failures.append("{}: inference unexpectedly unsupported: {}".format(
                comparison, preflight["reasons"]))
        if preflight["tools"]["rmats"]["supported"] != rmats_supported:
            failures.append("{}: rMATS support should be {}".format(comparison, rmats_supported))
        synthesis = (root / comparison / "synthesis.tsv").read_text()
        rmats_line = [line for line in synthesis.splitlines() if line.startswith("rmats\t")]
        if not rmats_supported and (not rmats_line or "\tunsupported\t" not in rmats_line[0]):
            failures.append("{}: synthesis should show rMATS as unsupported".format(comparison))

    print("shards checked: {}; skipping junctions recovered: {}/{}".format(
        len(shards), recovered, len(events)))
    if failures:
        for failure in failures:
            print("FAIL: " + failure, file=sys.stderr)
        sys.exit(1)
    print("umbrella smoke run passed")


if __name__ == "__main__":
    main()

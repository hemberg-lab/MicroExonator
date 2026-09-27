"""One row of raw QC values per run, as an immutable group shard.

No outlier calls are made here: a batch shard is too small to judge them.
That happens after all shards of a project and group are joined.
"""

import argparse
import gzip
import json
import re

if __package__:
    from src.reduce_featurecounts import read_summary
    from src.shard_guard import write_immutable_bundle
else:
    from reduce_featurecounts import read_summary
    from shard_guard import write_immutable_bundle


COLUMNS = ("run_id", "layout", "declared_strandedness", "inferred_strandedness",
           "sense_fraction", "hisat2_overall_alignment", "spliced_fragments",
           "junction_fragments", "ambiguous_strand", "featurecounts_assigned",
           "featurecounts_total", "featurecounts_assigned_fraction")


def hisat2_rate(path):
    """Overall alignment rate from a HISAT2 --new-summary file, as a fraction."""
    with open(path) as stream:
        for line in stream:
            match = re.search(r"Overall alignment rate:\s*([0-9.]+)%", line)
            if match:
                return round(float(match.group(1)) / 100, 4)
    raise ValueError("no overall alignment rate in {}".format(path))


def row(run):
    with open(run["junction_summary"]) as stream:
        junctions = json.load(stream)
    assigned = read_summary(run["featurecounts_summary"])
    total = sum(assigned.values())
    values = {
        "run_id": run["run_id"], "layout": run["layout"],
        "declared_strandedness": run["strandedness"],
        "inferred_strandedness": junctions["inferred_strandedness"],
        "sense_fraction": junctions["sense_fraction"],
        "hisat2_overall_alignment": hisat2_rate(run["hisat2_summary"]),
        "spliced_fragments": junctions["spliced_fragments"],
        "junction_fragments": junctions["junction_fragments"],
        "ambiguous_strand": junctions["ambiguous_strand"],
        "featurecounts_assigned": assigned.get("Assigned", 0),
        "featurecounts_total": total,
        "featurecounts_assigned_fraction": round(assigned.get("Assigned", 0) / total, 4) if total else None,
    }
    return "\t".join("NA" if values[column] is None else str(values[column])
                     for column in COLUMNS) + "\n"


def reduce_group(runs, prefix, reference_id, manifest_sha256):
    runs = sorted(runs, key=lambda run: run["run_id"])
    content = "\t".join(COLUMNS) + "\n" + "".join(row(run) for run in runs)
    write_immutable_bundle(
        {prefix + ".qc.tsv.gz": gzip.compress(content.encode(), mtime=0)},
        prefix + ".qc_checksums.json",
        {"reference_id": reference_id, "manifest_sha256": manifest_sha256,
         "runs": [run["run_id"] for run in runs]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--runs-json", required=True,
                        help="JSON list of {run_id, layout, strandedness, hisat2_summary, "
                             "junction_summary, featurecounts_summary}")
    args = parser.parse_args()
    with open(args.runs_json) as stream:
        runs = json.load(stream)
    reduce_group(runs, args.prefix, args.reference_id, args.manifest_sha256)


if __name__ == "__main__":
    main()

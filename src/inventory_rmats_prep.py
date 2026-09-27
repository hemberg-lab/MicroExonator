"""Keep one run's rMATS-turbo prep output as a durable, checksummed file.

rMATS `--task prep` writes `<tmp>/<timestamp>_<n>.rmats`. That file records
the BAM path it was made from; `--task post` must later be given the same BAM
path, even though the BAM itself is gone by then. This script copies the single
.rmats file to a stable run-named path and writes an inventory row with that
recorded BAM path, the read length and the layout.
"""

import argparse
import hashlib
import json
from pathlib import Path

if __package__:
    from src.shard_guard import write_immutable, write_immutable_bundle
else:
    from shard_guard import write_immutable, write_immutable_bundle


def single_prep_file(tmp_dir):
    files = sorted(Path(tmp_dir).glob("*.rmats"))
    if len(files) != 1:
        raise ValueError("expected one .rmats file in {}, found {}".format(tmp_dir, len(files)))
    return files[0]


def keep_run(tmp_dir, destination, record):
    """Copy the prep file and write its inventory JSON next to it."""
    content = single_prep_file(tmp_dir).read_bytes()
    record = dict(record, sha256=hashlib.sha256(content).hexdigest())
    write_immutable(destination, content)
    write_immutable(str(destination) + ".json",
                    (json.dumps(record, sort_keys=True, indent=2) + "\n").encode())


def inventory_text(record_paths):
    rows = []
    for path in sorted(record_paths):
        with open(path) as stream:
            rows.append(json.load(stream))
    columns = ("run_id", "bam", "layout", "read_length", "prep_file", "sha256")
    lines = ["\t".join(columns) + "\n"]
    for row in sorted(rows, key=lambda row: row["run_id"]):
        lines.append("\t".join(str(row[column]) for column in columns) + "\n")
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    keep = sub.add_parser("keep")
    keep.add_argument("--tmp", required=True)
    keep.add_argument("--destination", required=True)
    keep.add_argument("--run-id", required=True)
    keep.add_argument("--bam", required=True)
    keep.add_argument("--layout", required=True)
    keep.add_argument("--read-length", required=True, type=int)
    inventory = sub.add_parser("inventory")
    inventory.add_argument("--prefix", required=True)
    inventory.add_argument("--reference-id", required=True)
    inventory.add_argument("--manifest-sha256", required=True)
    inventory.add_argument("records", nargs="+")
    args = parser.parse_args()
    if args.command == "keep":
        keep_run(args.tmp, args.destination, {
            "run_id": args.run_id, "bam": args.bam, "layout": args.layout,
            "read_length": args.read_length, "prep_file": args.destination})
    else:
        write_immutable_bundle(
            {args.prefix + ".rmats_inventory.tsv": inventory_text(args.records).encode()},
            args.prefix + ".rmats_checksums.json",
            {"reference_id": args.reference_id, "manifest_sha256": args.manifest_sha256})


if __name__ == "__main__":
    main()

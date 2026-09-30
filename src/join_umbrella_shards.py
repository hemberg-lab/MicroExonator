"""Join group-by-batch shards of one kind into one deterministic matrix.

Each shard is checked against its checksum guard (content hash and
reference_id) before joining, so a shard edited after it was written, or one
from another reference bundle, is refused. Runs must not repeat across shards.

Kinds, their key columns and how rows missing from a shard are filled:

  junctions      chrom start end strand + label; max_anchor is the maximum,
                 multi_total and short_anchor_total are summed; runs filled 0
  featurecounts  gene_id length; runs filled 0
  salmon         gene_id or transcript_id; runs filled 0
  microexonator  ME; five columns per run, filled empty (sparse)
  whippet        Gene Node Coord Strand Type, in row order; shards must have
                 identical rows (same index), else the join is refused

With `runs` (from a comparison preflight, `--select`), only those runs'
columns are kept: a comparison can take some of a shard's runs, so selections
on metadata work without rebuilding shards. Sparse rows (MicroExonator) empty
in every kept run are dropped; junction summary columns still describe all
runs of the shard.

collapse_technical_runs() sums count columns of runs that belong to one
biological replicate, so technical runs never become independent replicates.
"""

import argparse
import gzip
import hashlib
import json
from collections import OrderedDict, defaultdict
from pathlib import Path


KINDS = {
    "junctions": {"key": 5, "summary": 3, "width": 1, "fill": "0"},
    "featurecounts": {"key": 2, "summary": 0, "width": 1, "fill": "0"},
    "salmon": {"key": 1, "summary": 0, "width": 1, "fill": "0"},
    "microexonator": {"key": 1, "summary": 0, "width": 5, "fill": ""},
    "whippet": {"key": 5, "summary": 0, "width": 6, "fill": None},
}


def _guard_for(shard):
    shard = str(shard)
    for suffix in (".tsv.gz", ".tsv"):
        if shard.endswith(suffix):
            stem = shard[:-len(suffix)]
            break
    else:
        raise ValueError("unexpected shard name: {}".format(shard))
    # X.junctions.tsv.gz -> X.junctions_checksums.json; Salmon's matrices share
    # one guard: X.salmon_counts.tsv.gz -> X.salmon_checksums.json.
    prefix, name = stem.rsplit(".", 1)
    candidates = [stem + "_checksums.json",
                  "{}.{}_checksums.json".format(prefix, name.split("_")[0])]
    for candidate in candidates:
        if Path(candidate).exists():
            return candidate
    raise ValueError("no checksum guard next to {}".format(shard))


def verify_shard(shard, reference_id):
    with open(_guard_for(shard)) as stream:
        guard = json.load(stream)
    if guard["metadata"].get("reference_id") != reference_id:
        raise ValueError("shard {} belongs to reference {}, not {}".format(
            shard, guard["metadata"].get("reference_id"), reference_id))
    digest = hashlib.sha256(Path(shard).read_bytes()).hexdigest()
    recorded = {Path(path).name: value for path, value in guard["sha256"].items()}
    if recorded.get(Path(shard).name) != digest:
        raise ValueError("shard {} does not match its checksum guard".format(shard))
    return digest


def _read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        header = next(stream).rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in stream]
    return header, rows


def _run_of(column, width):
    # multi-column kinds name columns <run>.<field>; field names have no dots
    return column if width == 1 else column.rsplit(".", 1)[0]


def _select(header, rows, fixed, width, runs):
    """Keep the fixed columns and the column blocks of `runs`."""
    columns = header[fixed:]
    if len(columns) % width:
        raise ValueError("shard columns do not divide into runs")
    keep = list(range(fixed))
    for start in range(fixed, len(header), width):
        if _run_of(header[start], width) in runs:
            keep.extend(range(start, start + width))
    return [header[i] for i in keep], [[row[i] for i in keep] for row in rows]


def join(shards, kind, reference_id, runs=None):
    spec = KINDS[kind]
    fixed = spec["key"] + spec["summary"]
    headers, tables, digests = [], [], {}
    for shard in sorted(shards):
        digests[str(shard)] = verify_shard(shard, reference_id)
        header, rows = _read(shard)
        if runs is not None:
            header, rows = _select(header, rows, fixed, spec["width"], set(runs))
            if len(header) == fixed:
                raise ValueError("shard {} holds none of the selected runs".format(shard))
            if spec["fill"] == "":
                rows = [row for row in rows if any(row[fixed:])]
        headers.append(header)
        tables.append(rows)
    if len({tuple(header[:fixed]) for header in headers}) != 1:
        raise ValueError("shards have different fixed columns")
    run_columns = [header[fixed:] for header in headers]
    flat = [column for columns in run_columns for column in columns]
    if len(flat) != len(set(flat)):
        raise ValueError("a run appears in more than one shard")

    if kind == "whippet":
        events = [[row[:fixed] for row in rows] for rows in tables]
        if any(event != events[0] for event in events[1:]):
            raise ValueError("Whippet shards have different rows: not the same index")
        order = sorted(range(len(tables)), key=lambda i: run_columns[i])
        lines = ["\t".join(headers[0][:fixed] + [c for i in order for c in run_columns[i]]) + "\n"]
        for r, event in enumerate(events[0]):
            lines.append("\t".join(event + [v for i in order for v in tables[i][r][fixed:]]) + "\n")
        return "".join(lines), digests

    merged = OrderedDict()
    for i, rows in enumerate(tables):
        for row in rows:
            key = tuple(row[:spec["key"]])
            entry = merged.setdefault(key, {"summary": None, "runs": {}})
            summary = row[spec["key"]:fixed]
            if spec["summary"]:
                if entry["summary"] is None:
                    entry["summary"] = [int(value) for value in summary]
                else:
                    entry["summary"] = [max(entry["summary"][0], int(summary[0]))] + [
                        a + int(b) for a, b in zip(entry["summary"][1:], summary[1:])]
            entry["runs"][i] = row[fixed:]
    order = sorted(range(len(tables)), key=lambda i: run_columns[i])
    lines = ["\t".join(headers[0][:fixed] + [c for i in order for c in run_columns[i]]) + "\n"]

    def sort_key(key):
        return (key[0], int(key[1]), int(key[2]), key[3]) if kind == "junctions" else key

    for key in sorted(merged, key=sort_key):
        entry = merged[key]
        values = list(key) + [str(value) for value in (entry["summary"] or [])]
        for i in order:
            values += entry["runs"].get(i, [spec["fill"]] * len(run_columns[i]))
        lines.append("\t".join(values) + "\n")
    return "".join(lines), digests


def collapse_technical_runs(text, kind, replicate_of):
    """Sum count columns by biological replicate (single-width count kinds only)."""
    spec = KINDS[kind]
    if spec["width"] != 1 or spec["fill"] != "0":
        raise ValueError("only count matrices can be collapsed by summing")
    fixed = spec["key"] + spec["summary"]
    lines = text.splitlines()
    header = lines[0].split("\t")
    runs = header[fixed:]
    replicates = sorted({replicate_of[run] for run in runs})
    members = defaultdict(list)
    for index, run in enumerate(runs):
        members[replicate_of[run]].append(fixed + index)
    out = ["\t".join(header[:fixed] + replicates) + "\n"]
    for line in lines[1:]:
        fields = line.split("\t")
        sums = []
        for replicate in replicates:
            total = sum(float(fields[i]) for i in members[replicate])
            sums.append(str(int(total)) if total.is_integer() else format(total, ".12g"))
        out.append("\t".join(fields[:fixed] + sums) + "\n")
    return "".join(out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kind", required=True, choices=sorted(KINDS))
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--collapse", help="JSON {run_id: biological_replicate_id}, or a "
                        "comparison preflight JSON, whose 'collapse' map is used")
    parser.add_argument("--select", help="comparison preflight JSON: keep only its runs")
    parser.add_argument("shards", nargs="+")
    args = parser.parse_args(argv)
    runs = None
    if args.select:
        with open(args.select) as stream:
            runs = sorted(json.load(stream)["collapse"])
    text, digests = join(args.shards, args.kind, args.reference_id, runs)
    if args.collapse:
        with open(args.collapse) as stream:
            mapping = json.load(stream)
        text = collapse_technical_runs(text, args.kind, mapping.get("collapse", mapping))
    with open(args.output, "wb") as stream:
        stream.write(gzip.compress(text.encode(), mtime=0))
    with open(args.provenance, "w") as stream:
        json.dump({"kind": args.kind, "reference_id": args.reference_id,
                   "shards": digests, "runs": runs, "output_sha256": hashlib.sha256(text.encode()).hexdigest()},
                  stream, sort_keys=True, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()

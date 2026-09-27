"""Group-by-batch shards of per-run MicroExonator and Whippet splicing outputs.

The shards keep every field the delta modules read, verbatim (no numeric
reformatting), so each run's delta input can be rebuilt from them:

  microexonator  one row per microexon; per run ME_coverages, excluding_covs,
                 PSI, CI_Lo and CI_Hi from `*.corrected.PSI.gz`
                 (read by src/me_delta.py). A microexon absent from a run's
                 sparse table is left empty for that run.
  whippet        Whippet's rows in their original order (whippet-delta pairs
                 files line by line); the event columns Gene, Node, Coord,
                 Strand, Type once, and per run Psi, CI_Width, CI_Lo,Hi,
                 Total_Reads, Complexity and Entropy (columns 6-11, the ones
                 whippet-delta reads). The path columns (Inc_Paths, Exc_Paths,
                 Edges) are not kept; rebuilt files carry NA there.

Per-run inputs are not deleted here.
"""

import argparse
import gzip

if __package__:
    from src.shard_guard import write_immutable_bundle
else:
    from shard_guard import write_immutable_bundle


ME_FIELDS = ("ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi")
ME_HEADER = ("sample", "ME") + ME_FIELDS
WHIPPET_EVENT = ("Gene", "Node", "Coord", "Strand", "Type")
WHIPPET_FIELDS = ("Psi", "CI_Width", "CI_Lo,Hi", "Total_Reads", "Complexity", "Entropy")
WHIPPET_PATHS = ("Inc_Paths", "Exc_Paths", "Edges")


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _gzip_text(content):
    return gzip.compress(content.encode(), mtime=0)


# ------------------------------------------------------------ MicroExonator

def read_me_table(path):
    rows = {}
    with _open(path) as stream:
        header = tuple(next(stream).rstrip("\n").split("\t"))
        if header != ME_HEADER:
            raise ValueError("unexpected corrected PSI header in {}: {}".format(path, header))
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if fields[1] in rows:
                raise ValueError("duplicate microexon {} in {}".format(fields[1], path))
            rows[fields[1]] = fields[2:]
    return rows


def me_shard_text(tables):
    runs = sorted(tables)
    parsed = {run: read_me_table(tables[run]) for run in runs}
    microexons = sorted(set().union(*(set(rows) for rows in parsed.values())))
    header = ["ME"] + ["{}.{}".format(run, field) for run in runs for field in ME_FIELDS]
    lines = ["\t".join(header) + "\n"]
    empty = [""] * len(ME_FIELDS)
    for me in microexons:
        values = [me]
        for run in runs:
            values += parsed[run].get(me, empty)
        lines.append("\t".join(values) + "\n")
    return "".join(lines)


def _runs_from_header(header, fields, first):
    width = len(fields)
    columns = header[first:]
    if len(columns) % width:
        raise ValueError("shard columns do not divide into runs")
    runs = []
    for i in range(0, len(columns), width):
        run = columns[i].rsplit("." + fields[0], 1)[0]
        expected = ["{}.{}".format(run, field) for field in fields]
        if columns[i:i + width] != expected:
            raise ValueError("malformed shard columns for run {}".format(run))
        runs.append(run)
    return runs


def rebuild_me_table(shard_path, run):
    """The run's corrected PSI table, rows sorted by microexon."""
    with _open(shard_path) as stream:
        header = next(stream).rstrip("\n").split("\t")
        runs = _runs_from_header(header, ME_FIELDS, 1)
        if run not in runs:
            raise ValueError("run {} not in {}".format(run, shard_path))
        offset = 1 + runs.index(run) * len(ME_FIELDS)
        lines = ["\t".join(ME_HEADER) + "\n"]
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            values = fields[offset:offset + len(ME_FIELDS)]
            if any(values):
                lines.append("\t".join([run, fields[0]] + values) + "\n")
    return "".join(lines)


# ------------------------------------------------------------ Whippet

def read_whippet(path):
    events, values = [], []
    with _open(path) as stream:
        header = next(stream).rstrip("\n").split("\t")
        if tuple(header[:11]) != WHIPPET_EVENT + WHIPPET_FIELDS:
            raise ValueError("unexpected Whippet psi header in {}".format(path))
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            events.append(fields[:5])
            values.append(fields[5:11] + [""] * max(0, 11 - len(fields)))
    return events, values


def whippet_shard_text(psi_files):
    runs = sorted(psi_files)
    parsed = {run: read_whippet(psi_files[run]) for run in runs}
    events = parsed[runs[0]][0]
    for run in runs[1:]:
        if parsed[run][0] != events:
            raise ValueError("Whippet rows differ for run {}: not the same index".format(run))
    header = list(WHIPPET_EVENT) + ["{}.{}".format(run, field) for run in runs for field in WHIPPET_FIELDS]
    lines = ["\t".join(header) + "\n"]
    for i, event in enumerate(events):
        row = list(event)
        for run in runs:
            row += parsed[run][1][i][:len(WHIPPET_FIELDS)]
        lines.append("\t".join(row) + "\n")
    return "".join(lines)


def rebuild_whippet_psi(shard_path, run):
    with _open(shard_path) as stream:
        header = next(stream).rstrip("\n").split("\t")
        runs = _runs_from_header(header, WHIPPET_FIELDS, len(WHIPPET_EVENT))
        if run not in runs:
            raise ValueError("run {} not in {}".format(run, shard_path))
        offset = len(WHIPPET_EVENT) + runs.index(run) * len(WHIPPET_FIELDS)
        lines = ["\t".join(WHIPPET_EVENT + WHIPPET_FIELDS + WHIPPET_PATHS) + "\n"]
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            lines.append("\t".join(fields[:5] + fields[offset:offset + len(WHIPPET_FIELDS)]
                                   + ["NA"] * len(WHIPPET_PATHS)) + "\n")
    return "".join(lines)


# ------------------------------------------------------------ CLI

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("microexonator", "whippet"):
        shard = sub.add_parser(name)
        shard.add_argument("--prefix", required=True)
        shard.add_argument("--reference-id", required=True)
        shard.add_argument("--manifest-sha256", required=True)
        shard.add_argument("--input", action="append", required=True, help="run_id=path")
    rebuild = sub.add_parser("rebuild")
    rebuild.add_argument("--kind", choices=("microexonator", "whippet"), required=True)
    rebuild.add_argument("--shard", required=True)
    rebuild.add_argument("--run", required=True)
    rebuild.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    if args.command == "rebuild":
        text = (rebuild_me_table if args.kind == "microexonator" else rebuild_whippet_psi)(
            args.shard, args.run)
        with _open(args.output, "wt") as stream:
            stream.write(text)
        return
    inputs = dict(value.split("=", 1) for value in args.input)
    if len(inputs) != len(args.input):
        parser.error("duplicate run ID")
    text = (me_shard_text if args.command == "microexonator" else whippet_shard_text)(inputs)
    suffix = ".microexonator.tsv.gz" if args.command == "microexonator" else ".whippet.tsv.gz"
    write_immutable_bundle({args.prefix + suffix: _gzip_text(text)},
                           args.prefix + suffix.replace(".tsv.gz", "_checksums.json"),
                           {"reference_id": args.reference_id,
                            "manifest_sha256": args.manifest_sha256, "runs": sorted(inputs)})


if __name__ == "__main__":
    main()

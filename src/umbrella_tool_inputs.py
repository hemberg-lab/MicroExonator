"""Write the input files that each comparison tool expects, from joined matrices.

Subcommands
  delta      rebuild per-run MicroExonator or Whippet delta inputs from a joined
             splicing matrix (same layout as a shard)
  leafcutter per-replicate LeafCutter junction files from the joined junction
             matrix (already summed per replicate), plus the groups file
  suppa      SUPPA2 TPM tables for group A and B from joined transcript TPM
             (runs averaged per biological replicate)

LeafCutter junction files are six columns: chrom, start, end, '.', count,
strand, in the convention of clustering/leafcutter_cluster.py, which reads
(A, B) and keys the intron as (A, B + 1): A is the last base of the upstream
exon and B the last intron base (1-based). The intron key is then (upstream
exon end, downstream exon start), matching GTF exons, as leafcutter_ds.R
expects. Junction shards store the intron itself (first and last intron base),
so A = intron start - 1. LeafCutter names samples after the junction files,
so the groups file uses the file names (replicate + '.junc').
"""

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

if __package__:
    from src.reduce_splicing import rebuild_me_table, rebuild_whippet_psi
else:
    from reduce_splicing import rebuild_me_table, rebuild_whippet_psi


def _read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        header = next(stream).rstrip("\n").split("\t")
        return header, [line.rstrip("\n").split("\t") for line in stream]


def write_delta_inputs(kind, joined, runs, outputs):
    rebuild = rebuild_me_table if kind == "microexonator" else rebuild_whippet_psi
    for run, output in zip(runs, outputs):
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(output, "wt") as stream:
            stream.write(rebuild(joined, run))


def leafcutter_files(joined, preflight, directory):
    """Returns the list of junction files written; also writes groups.txt."""
    header, rows = _read(joined)
    replicates = preflight["replicates"]["a"] + preflight["replicates"]["b"]
    index = {name: i for i, name in enumerate(header)}
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    for replicate in replicates:
        path = directory / (replicate + ".junc")
        with open(path, "w") as stream:
            for row in rows:
                count = int(row[index[replicate]])
                if count:
                    stream.write("{}\t{}\t{}\t.\t{}\t{}\n".format(
                        row[0], int(row[1]) - 1, row[2], count, row[3]))
        written.append(str(path))
    # absolute paths: LeafCutter is run from inside the directory
    (directory / "juncfiles.txt").write_text("".join(
        str(Path(path).resolve()) + "\n" for path in written))
    (directory / "groups.txt").write_text("".join(
        "{}.junc\t{}\n".format(replicate, side) for side in ("a", "b")
        for replicate in preflight["replicates"][side]))
    return written


def suppa_tables(joined, preflight):
    """{'a': text, 'b': text}: SUPPA2 TPM tables, one column per replicate."""
    header, rows = _read(joined)
    runs = header[1:]
    members = defaultdict(list)
    for column, run in enumerate(runs, start=1):
        members[preflight["collapse"][run]].append(column)
    tables = {}
    for side in ("a", "b"):
        replicates = preflight["replicates"][side]
        lines = ["\t".join(replicates) + "\n"]      # SUPPA: no header for the ID column
        for row in rows:
            values = [sum(float(row[c]) for c in members[r]) / len(members[r]) for r in replicates]
            lines.append(row[0] + "\t" + "\t".join(format(v, ".6g") for v in values) + "\n")
        tables[side] = "".join(lines)
    return tables


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    delta = sub.add_parser("delta")
    delta.add_argument("--kind", choices=("microexonator", "whippet"), required=True)
    delta.add_argument("--joined", required=True)
    delta.add_argument("--run", action="append", required=True, help="run_id=output path")
    leafcutter = sub.add_parser("leafcutter")
    leafcutter.add_argument("--joined", required=True)
    leafcutter.add_argument("--preflight", required=True)
    leafcutter.add_argument("--directory", required=True)
    suppa = sub.add_parser("suppa")
    suppa.add_argument("--joined", required=True)
    suppa.add_argument("--preflight", required=True)
    suppa.add_argument("--output-a", required=True)
    suppa.add_argument("--output-b", required=True)
    args = parser.parse_args(argv)
    if args.command == "delta":
        pairs = [value.split("=", 1) for value in args.run]
        write_delta_inputs(args.kind, args.joined, [p[0] for p in pairs], [p[1] for p in pairs])
        return
    with open(args.preflight) as stream:
        preflight = json.load(stream)
    if args.command == "leafcutter":
        leafcutter_files(args.joined, preflight, args.directory)
    else:
        tables = suppa_tables(args.joined, preflight)
        Path(args.output_a).write_text(tables["a"])
        Path(args.output_b).write_text(tables["b"])


if __name__ == "__main__":
    main()

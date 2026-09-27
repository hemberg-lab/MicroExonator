"""Merge per-run featureCounts gene counts into one immutable group shard."""

import argparse
import gzip

if __package__:
    from src.shard_guard import write_immutable_bundle
else:
    from shard_guard import write_immutable_bundle


def read_counts(path):
    """featureCounts output: a comment line, a header, then Geneid ... Length count."""
    genes = []
    with open(path) as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Geneid":
                if len(fields) != 7:
                    raise ValueError("expected one BAM per featureCounts table: {}".format(path))
                continue
            genes.append((fields[0], fields[5], fields[6]))
    if not genes:
        raise ValueError("empty featureCounts table: {}".format(path))
    return genes


def read_summary(path):
    values = {}
    with open(path) as stream:
        next(stream)
        for line in stream:
            status, count = line.rstrip("\n").split("\t")
            values[status] = int(count)
    return values


def shard_text(counts_by_run):
    runs = sorted(counts_by_run)
    tables = {run: read_counts(counts_by_run[run]) for run in runs}
    reference = [(gene, length) for gene, length, _ in tables[runs[0]]]
    for run in runs[1:]:
        if [(gene, length) for gene, length, _ in tables[run]] != reference:
            raise ValueError("featureCounts gene order or length differs for {}".format(run))
    lines = ["gene_id\tlength\t" + "\t".join(runs) + "\n"]
    for i, (gene, length) in enumerate(reference):
        lines.append("{}\t{}\t{}\n".format(gene, length, "\t".join(
            tables[run][i][2] for run in runs)))
    return "".join(lines)


def reduce_group(counts_by_run, prefix, reference_id, manifest_sha256):
    if not counts_by_run:
        raise ValueError("featureCounts shard requires at least one run")
    write_immutable_bundle(
        {prefix + ".featurecounts.tsv.gz": gzip.compress(shard_text(counts_by_run).encode(), mtime=0)},
        prefix + ".featurecounts_checksums.json",
        {"reference_id": reference_id, "manifest_sha256": manifest_sha256,
         "runs": sorted(counts_by_run)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--counts", action="append", required=True, help="run_id=path")
    args = parser.parse_args()
    counts = dict(value.split("=", 1) for value in args.counts)
    if len(counts) != len(args.counts):
        parser.error("duplicate run ID")
    reduce_group(counts, args.prefix, args.reference_id, args.manifest_sha256)


if __name__ == "__main__":
    main()

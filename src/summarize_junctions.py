"""Junction evidence from an independent spliced alignment, and index capture.

Subcommands
  extract       SAM on stdin (samtools view of one run's BAM) -> per-run junction
                table and summary JSON.
  strandedness  print the library strandedness argument for a tool.
  catalog       reference annotations -> labelled intron catalog.
  shard         per-run tables of one group and batch -> immutable junction shard.
  capture       junction shard + catalog -> immutable per-site capture shard.

Coordinates: introns are 1-based and closed (first and last intron base).

Counting rules
  - Only primary, mapped, QC-passed alignments count (flags 0x4, 0x100, 0x200
    and 0x800 are skipped).
  - A fragment counts once per junction: when both mates span the same
    junction it is one count, with the larger of the two anchors.
  - Anchor = the smaller number of aligned bases on either side of the intron,
    up to the neighbouring intron or the read end. Fragments with an anchor of
    at least --min-anchor (default 8) count as unique (NH 1) or multimapping
    (NH > 1); shorter ones count only as short_anchor.
  - The junction strand is the XS tag. Spliced reads without XS are counted as
    ambiguous in the summary and left out of the table.
  - Strandedness: for spliced fragments with XS, whether the fragment's
    orientation (read 2 flipped) agrees with XS.

Capture: Whippet splices any annotated donor to any annotated acceptor of the
same gene, so a junction is captured by the index when both of its ends are
intron ends in the Whippet annotation (label index_intron or index_sites).
Site sets are kept per chromosome and strand, not per gene.
"""

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict

if __package__:
    from src.shard_guard import write_immutable_bundle
else:
    from shard_guard import write_immutable_bundle


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x800
LABELS = ("index_intron", "index_sites", "annotation", "hint", "novel")
CAPTURED = {"index_intron", "index_sites"}
MIN_STRAND_EVIDENCE = 1000
TOOL_STRAND = {
    "featurecounts": {"unstranded": "0", "secondstrand": "1", "firststrand": "2"},
}
SHARD_FIXED = ("chrom", "start", "end", "strand", "label", "max_anchor",
               "multi_total", "short_anchor_total")


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _gzip_text(content):
    return gzip.compress(content.encode(), mtime=0)


def read_junctions(pos, cigar):
    """Return [(start, end, anchor)] for every N operation of one alignment."""
    ref = pos
    blocks, introns, aligned = [], [], 0
    for length, op in CIGAR.findall(cigar):
        length = int(length)
        if op in "M=X":
            aligned += length
            ref += length
        elif op == "D":
            ref += length
        elif op == "N":
            introns.append((ref, ref + length - 1))
            blocks.append(aligned)
            aligned = 0
            ref += length
    blocks.append(aligned)
    return [(start, end, min(blocks[i], blocks[i + 1]))
            for i, (start, end) in enumerate(introns)]


def infer_strandedness(sense, antisense):
    total = sense + antisense
    if total < MIN_STRAND_EVIDENCE:
        return "undetermined"
    fraction = sense / total
    if fraction >= 0.8:
        return "secondstrand"
    if fraction <= 0.2:
        return "firststrand"
    return "unstranded"


class JunctionCounter:
    def __init__(self, min_anchor=8):
        self.min_anchor = min_anchor
        self.table = defaultdict(lambda: [0, 0, 0, 0])  # unique, multi, short, max_anchor
        self.pending = {}
        self.summary = {"spliced_fragments": 0, "junction_fragments": 0,
                        "ambiguous_strand": 0, "sense": 0, "antisense": 0}

    def _record(self, chrom, junctions, nh, orientation):
        """junctions: {(start, end, strand): anchor} for one fragment."""
        if not junctions:
            return
        self.summary["spliced_fragments"] += 1
        for (start, end, strand), anchor in sorted(junctions.items()):
            if strand not in ("+", "-"):
                self.summary["ambiguous_strand"] += 1
                continue
            self.summary["junction_fragments"] += 1
            self.summary["sense" if orientation == strand else "antisense"] += 1
            values = self.table[(chrom, start, end, strand)]
            if anchor >= self.min_anchor:
                values[0 if nh == 1 else 1] += 1
            else:
                values[2] += 1
            values[3] = max(values[3], anchor)

    def add(self, line):
        if line.startswith("@"):
            return
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            return
        flag = int(fields[1])
        if flag & SKIP_FLAGS:
            return
        chrom, pos, cigar = fields[2], int(fields[3]), fields[5]
        tags = {tag[:2]: tag[5:] for tag in fields[11:] if len(tag) > 5}
        nh = int(tags.get("NH", "1"))
        strand = tags.get("XS", ".")
        junctions = {}
        for start, end, anchor in read_junctions(pos, cigar):
            key = (start, end, strand)
            junctions[key] = max(anchor, junctions.get(key, 0))
        orientation = "-" if flag & 0x10 else "+"
        if flag & 0x80:
            orientation = "+" if orientation == "-" else "-"
        qname = fields[0]
        paired = bool(flag & 0x1) and not flag & 0x8 and fields[6] in ("=", chrom)
        if paired and qname in self.pending:
            mate_chrom, mate_junctions, mate_nh, mate_orientation = self.pending.pop(qname)
            for key, anchor in junctions.items():
                mate_junctions[key] = max(anchor, mate_junctions.get(key, 0))
            self._record(mate_chrom, mate_junctions, mate_nh, mate_orientation)
        elif paired and int(fields[7]) >= pos:
            # The mate comes later in coordinate order; hold this read until then.
            self.pending[qname] = (chrom, junctions, nh, orientation)
        else:
            self._record(chrom, junctions, nh, orientation)

    def finish(self):
        # Mates that never appeared (filtered as secondary or QC-failed).
        for qname in sorted(self.pending):
            self._record(*self.pending[qname])
        self.pending.clear()
        informative = self.summary["sense"] + self.summary["antisense"]
        self.summary["sense_fraction"] = (
            round(self.summary["sense"] / informative, 4) if informative else None)
        self.summary["inferred_strandedness"] = infer_strandedness(
            self.summary["sense"], self.summary["antisense"])
        return self

    def table_text(self):
        lines = ["chrom\tstart\tend\tstrand\tunique\tmulti\tshort_anchor\tmax_anchor\n"]
        for key in sorted(self.table):
            lines.append("\t".join(map(str, key + tuple(self.table[key]))) + "\n")
        return "".join(lines)


def strandedness_argument(tool, declared, summary_path, layout="SE"):
    """Declared strandedness wins unless it is auto; then use the inferred one."""
    if declared == "auto":
        with open(summary_path) as stream:
            declared = json.load(stream)["inferred_strandedness"]
        if declared == "undetermined":
            declared = "unstranded"
    if tool == "hisat2":
        return {"unstranded": "", "secondstrand": "F", "firststrand": "R"}[declared] * (
            2 if layout == "PE" else 1)
    return TOOL_STRAND[tool][declared]


# ---------------------------------------------------------------- catalog

def _gtf_introns(path):
    exons = defaultdict(list)
    with _open(path) as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if match:
                exons[(match.group(1), fields[0], fields[6])].append(
                    (int(fields[3]), int(fields[4])))
    introns = set()
    for (transcript, chrom, strand), blocks in exons.items():
        blocks.sort()
        for (_, left_end), (right_start, _) in zip(blocks, blocks[1:]):
            if right_start - left_end > 1:
                introns.add((chrom, left_end + 1, right_start - 1, strand))
    return introns


def _hint_introns(path):
    """HISAT2 splice-site file: chrom, 0-based last exon base, 0-based first exon base, strand."""
    introns = set()
    with _open(path) as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 4 and not line.startswith("#"):
                introns.add((fields[0], int(fields[1]) + 2, int(fields[2]), fields[3]))
    return introns


def build_catalog(whippet_gtf, annotation_gtf, splice_sites):
    index = _gtf_introns(whippet_gtf)
    annotation = _gtf_introns(annotation_gtf)
    hints = _hint_introns(splice_sites)
    lines = ["chrom\tstart\tend\tstrand\tlabel\n"]
    for key in sorted(index | annotation | hints):
        label = "index_intron" if key in index else (
            "annotation" if key in annotation else "hint")
        lines.append("{}\t{}\t{}\t{}\t{}\n".format(*key, label))
    return "".join(lines)


class Catalog:
    def __init__(self, path):
        self.labels = {}
        self.left = defaultdict(set)
        self.right = defaultdict(set)
        with _open(path) as stream:
            next(stream)
            for line in stream:
                chrom, start, end, strand, label = line.rstrip("\n").split("\t")
                key = (chrom, int(start), int(end), strand)
                self.labels[key] = label
                if label == "index_intron":
                    self.left[(chrom, strand)].add(key[1])
                    self.right[(chrom, strand)].add(key[2])

    def label(self, key):
        label = self.labels.get(key)
        if label == "index_intron":
            return label
        chrom, start, end, strand = key
        if start in self.left[(chrom, strand)] and end in self.right[(chrom, strand)]:
            return "index_sites"
        return label or "novel"


# ---------------------------------------------------------------- shard

def _read_run_table(path):
    rows = {}
    with _open(path) as stream:
        header = next(stream).rstrip("\n").split("\t")
        if header[:8] != ["chrom", "start", "end", "strand", "unique", "multi",
                          "short_anchor", "max_anchor"]:
            raise ValueError("invalid junction table: {}".format(path))
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            rows[(fields[0], int(fields[1]), int(fields[2]), fields[3])] = tuple(
                int(value) for value in fields[4:8])
    return rows


def shard_text(tables, catalog):
    """tables: {run_id: path}. Unique counts per run; other counts summed."""
    runs = sorted(tables)
    parsed = {run: _read_run_table(tables[run]) for run in runs}
    keys = sorted(set().union(*(set(rows) for rows in parsed.values())))
    lines = ["\t".join(SHARD_FIXED + tuple(runs)) + "\n"]
    for key in keys:
        values = [parsed[run].get(key, (0, 0, 0, 0)) for run in runs]
        fixed = [key[0], key[1], key[2], key[3], catalog.label(key),
                 max(value[3] for value in values),
                 sum(value[1] for value in values),
                 sum(value[2] for value in values)]
        lines.append("\t".join(map(str, fixed + [value[0] for value in values])) + "\n")
    return "".join(lines)


def read_shard(path):
    with _open(path) as stream:
        header = next(stream).rstrip("\n").split("\t")
        if tuple(header[:len(SHARD_FIXED)]) != SHARD_FIXED:
            raise ValueError("invalid junction shard: {}".format(path))
        runs = header[len(SHARD_FIXED):]
        rows = []
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], int(fields[1]), int(fields[2]), fields[3])
            rows.append((key, fields[4], [int(value) for value in fields[len(SHARD_FIXED):]]))
    return runs, rows


# ---------------------------------------------------------------- capture

def capture_text(shard_path, catalog):
    """Per index splice site and run: captured and uncaptured junction reads.

    Only sites with at least one uncaptured read in some run are written, with
    the strongest uncaptured (competing) junction. A site with no junction
    reads in a run is 'unassessable' for that run.
    """
    runs, rows = read_shard(shard_path)
    sites = defaultdict(lambda: {"captured": [0] * len(runs), "uncaptured": [0] * len(runs),
                                 "competitor": None, "competitor_reads": 0})
    for key, label, counts in rows:
        chrom, start, end, strand = key
        captured = label in CAPTURED
        for side, position, index_set in (("left", start, catalog.left),
                                          ("right", end, catalog.right)):
            if position not in index_set[(chrom, strand)]:
                continue
            site = sites[(chrom, position, strand, side)]
            bucket = site["captured" if captured else "uncaptured"]
            for i, count in enumerate(counts):
                bucket[i] += count
            total = sum(counts)
            if not captured and total > site["competitor_reads"]:
                site["competitor"] = "{}:{}-{}:{}".format(*key)
                site["competitor_reads"] = total
    lines = ["chrom\tposition\tstrand\tside\tcompetitor\tcompetitor_reads\t" + "\t".join(
        "{0}.captured\t{0}.uncaptured\t{0}.status".format(run) for run in runs) + "\n"]
    for key in sorted(sites):
        site = sites[key]
        if not any(site["uncaptured"]):
            continue
        per_run = []
        for captured, uncaptured in zip(site["captured"], site["uncaptured"]):
            total = captured + uncaptured
            status = "unassessable" if total == 0 else "{:.4g}".format(captured / total)
            per_run += [str(captured), str(uncaptured), status]
        lines.append("\t".join(map(str, key)) + "\t{}\t{}\t".format(
            site["competitor"], site["competitor_reads"]) + "\t".join(per_run) + "\n")
    return "".join(lines)


def run_capture_rates(shard_path):
    """Share of unique junction reads captured by the index, per run."""
    runs, rows = read_shard(shard_path)
    captured, total = [0] * len(runs), [0] * len(runs)
    for _, label, counts in rows:
        for i, count in enumerate(counts):
            total[i] += count
            if label in CAPTURED:
                captured[i] += count
    return {run: (round(c / t, 4) if t else None) for run, c, t in zip(runs, captured, total)}


# ---------------------------------------------------------------- CLI

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    extract = sub.add_parser("extract")
    extract.add_argument("--table", required=True)
    extract.add_argument("--summary", required=True)
    extract.add_argument("--min-anchor", type=int, default=8)

    strand = sub.add_parser("strandedness")
    strand.add_argument("--tool", required=True, choices=("featurecounts", "hisat2"))
    strand.add_argument("--declared", required=True)
    strand.add_argument("--summary")
    strand.add_argument("--layout", default="SE")

    catalog = sub.add_parser("catalog")
    catalog.add_argument("--whippet-gtf", required=True)
    catalog.add_argument("--annotation-gtf", required=True)
    catalog.add_argument("--splice-sites", required=True)
    catalog.add_argument("--output", required=True)

    shard = sub.add_parser("shard")
    shard.add_argument("--catalog", required=True)
    shard.add_argument("--prefix", required=True)
    shard.add_argument("--reference-id", required=True)
    shard.add_argument("--manifest-sha256", required=True)
    shard.add_argument("--table", action="append", required=True, help="run_id=path")

    capture = sub.add_parser("capture")
    capture.add_argument("--catalog", required=True)
    capture.add_argument("--shard", required=True)
    capture.add_argument("--prefix", required=True)
    capture.add_argument("--reference-id", required=True)
    capture.add_argument("--manifest-sha256", required=True)

    args = parser.parse_args(argv)
    if args.command == "extract":
        counter = JunctionCounter(args.min_anchor)
        for line in sys.stdin:
            counter.add(line)
        counter.finish()
        with _open(args.table, "wt") as stream:
            stream.write(counter.table_text())
        with open(args.summary, "w") as stream:
            json.dump(counter.summary, stream, sort_keys=True, indent=2)
            stream.write("\n")
    elif args.command == "strandedness":
        print(strandedness_argument(args.tool, args.declared, args.summary, args.layout))
    elif args.command == "catalog":
        with open(args.output, "wb") as stream:
            stream.write(_gzip_text(build_catalog(
                args.whippet_gtf, args.annotation_gtf, args.splice_sites)))
    elif args.command == "shard":
        tables = dict(value.split("=", 1) for value in args.table)
        if len(tables) != len(args.table):
            parser.error("duplicate run ID")
        write_immutable_bundle(
            {args.prefix + ".junctions.tsv.gz": _gzip_text(shard_text(tables, Catalog(args.catalog)))},
            args.prefix + ".junctions_checksums.json",
            {"reference_id": args.reference_id, "manifest_sha256": args.manifest_sha256,
             "runs": sorted(tables)})
    elif args.command == "capture":
        write_immutable_bundle(
            {args.prefix + ".capture.tsv.gz": _gzip_text(capture_text(args.shard, Catalog(args.catalog))),
             args.prefix + ".capture_rates.json": (json.dumps(
                 run_capture_rates(args.shard), sort_keys=True, indent=2) + "\n").encode()},
            args.prefix + ".capture_checksums.json",
            {"reference_id": args.reference_id, "manifest_sha256": args.manifest_sha256})


if __name__ == "__main__":
    main()

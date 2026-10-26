"""Frozen circRNA catalogue for the circRNA modules (CIRIquant and the support tables).

Everything runs on the umbrella assembly (hg38). Sources on another assembly (circBase is
hg19) are lifted first with `lift`, boundary by boundary through a UCSC chain file.

Input: BED6 lists of back-spliced junctions (BSJ), standard BED coordinates: column 2 is
the 0-based first base of the circRNA, column 3 its last base (1-based, inclusive),
column 6 the strand. Several sources can be concatenated (e.g. circAtlas hg38 and circBase
lifted to hg38); names of rows with the same BSJ are merged.

Output (`build`):
  catalogue.bed       CIRIquant's BED: chrom, 1-based start, end, circ_id, '.', strand.
                      CIRIquant reads column 2 as 1-based (BedParser: length = end - start + 1)
                      and needs circ_id to begin with "chrom:" (it splits on ':'), so
                      circ_id = chrom:start|end, its own convention.
  annotation.tsv.gz   circ_id, chrom, start, end, strand, span, source names, splice motif
  report.json         counts, drops by reason, pseudo-reference size, catalogue_id

Rows are dropped, with a reason, when the contig is not in the genome, coordinates are
outside it or empty, the strand is not + or -, the same BSJ is listed on both strands
(CIRIquant keys circRNAs by chrom:start|end only), the span exceeds max_span (0 = no cap),
or more than half of the pseudo-circular sequence would be N (CIRIquant skips those).
CIRIquant's pseudo-reference is the genomic span written twice, so the build refuses
catalogues whose sum of 2 x span exceeds max_pseudo_bases.
"""

import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


class Fasta:
    """Random access to an uncompressed FASTA through its samtools .fai."""

    def __init__(self, fasta, fai=None):
        self.fasta = str(fasta)
        self.index = {}
        with open(fai or self.fasta + ".fai") as stream:
            for line in stream:
                name, length, offset, line_bases, line_bytes = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(line_bases), int(line_bytes))
        self._handle = None

    def length(self, chrom):
        return self.index[chrom][0]

    def fetch(self, chrom, start0, end):
        """Sequence of [start0, end) on the + strand, upper case."""
        length, offset, line_bases, line_bytes = self.index[chrom]
        start0, end = max(0, start0), min(length, end)
        if end <= start0:
            return ""
        if self._handle is None:
            self._handle = open(self.fasta, "rb")
        first = offset + (start0 // line_bases) * line_bytes + start0 % line_bases
        last = offset + ((end - 1) // line_bases) * line_bytes + (end - 1) % line_bases
        self._handle.seek(first)
        raw = self._handle.read(last - first + 1).decode()
        return raw.replace("\n", "").replace("\r", "").upper()


def splice_motif(fasta, chrom, start1, end, strand):
    """Intron dinucleotides flanking the BSJ, read in the transcript's direction.

    The BSJ joins the donor after `end` to the acceptor before `start1`. On + the
    canonical pair is GT (end+1..end+2) and AG (start1-2..start1-1)."""
    donor = fasta.fetch(chrom, end, end + 2)          # end+1, end+2 (1-based)
    acceptor = fasta.fetch(chrom, start1 - 3, start1 - 1)  # start1-2, start1-1
    if strand == "-":
        donor, acceptor = acceptor.translate(COMPLEMENT)[::-1], donor.translate(COMPLEMENT)[::-1]
    return "{}-{}".format(donor, acceptor)


def circ_id(chrom, start1, end):
    return "{}:{}|{}".format(chrom, start1, end)


def build(source, fasta, out_dir, max_span=0, max_pseudo_bases=2_000_000_000, max_n_fraction=0.5,
          settings_extra=None, annotate=()):
    """`source` rows are quantified; `annotate` rows only add their names to the same BSJ
    (same coordinates and strand) when a quantified source already has it."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genome = Fasta(fasta)
    drops = {}
    rows = {}                  # (chrom, start1, end) -> {strand: set(names)}
    total = 0

    def drop(reason):
        drops[reason] = drops.get(reason, 0) + 1

    sources = [source] if isinstance(source, (str, Path)) else list(source)

    def lines():
        for path in sources:
            with _open(path) as stream:
                yield from stream

    for line in lines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        total += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            drop("fewer_than_6_columns")
            continue
        chrom, strand = fields[0], fields[5]
        try:
            start0, end = int(fields[1]), int(fields[2])
        except ValueError:
            drop("non_integer_coordinates")
            continue
        if strand not in ("+", "-"):
            drop("strand_not_plus_or_minus")
            continue
        if chrom not in genome.index:
            drop("contig_not_in_genome")
            continue
        if start0 < 0 or end <= start0 or end > genome.length(chrom):
            drop("outside_contig_or_empty")
            continue
        if max_span and end - start0 > max_span:
            drop("span_above_max_span")
            continue
        rows.setdefault((chrom, start0 + 1, end), {}).setdefault(strand, set()).add(fields[3])

    annotations = 0
    for path in ([annotate] if isinstance(annotate, (str, Path)) else list(annotate)):
        with _open(path) as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or line.startswith("#"):
                    continue
                key = (fields[0], int(fields[1]) + 1, int(fields[2]))
                if key in rows and fields[5] in rows[key]:
                    rows[key][fields[5]].add(fields[3])
                    annotations += 1

    kept = []
    for (chrom, start1, end), strands in sorted(rows.items()):
        if len(strands) > 1:
            drop("same_bsj_on_both_strands")
            continue
        strand, names = next(iter(strands.items()))
        sequence = genome.fetch(chrom, start1 - 1, end)
        if sequence.count("N") > len(sequence) * max_n_fraction:
            drop("mostly_n")
            continue
        kept.append((chrom, start1, end, strand, sorted(names)))

    pseudo_bases = sum(2 * (end - start1 + 1) for _, start1, end, _, _ in kept)
    settings = dict(settings_extra or {}, max_span=max_span, max_n_fraction=max_n_fraction,
                    coordinates="BED6 input (0-based start); CIRIquant BED output (1-based start)")
    bed = "".join("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start1, end, circ_id(chrom, start1, end), strand)
                  for chrom, start1, end, strand, _ in kept)
    catalogue_id = hashlib.sha256(json.dumps(settings, sort_keys=True).encode() + bed.encode()).hexdigest()[:16]
    motifs = {}
    report = {"catalogue_id": catalogue_id, "source_rows": total, "kept": len(kept), "annotation_names": annotations,
              "dropped": dict(sorted(drops.items())), "pseudo_reference_bases": pseudo_bases,
              "max_pseudo_bases": max_pseudo_bases, "settings": settings}
    if pseudo_bases > max_pseudo_bases:
        report["refused"] = ("the pseudo-circular reference would hold {} bases (2 x span summed), above "
                             "circrna_max_pseudo_bases {}; filter the catalogue or raise the limit".format(
                                 pseudo_bases, max_pseudo_bases))
        (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        raise SystemExit(report["refused"])
    (out / "catalogue.bed").write_text(bed)
    with gzip.open(out / "annotation.tsv.gz", "wt") as stream:
        stream.write("circ_id\tchrom\tstart\tend\tstrand\tspan\tsource_names\tsplice_motif\n")
        for chrom, start1, end, strand, names in kept:
            motif = splice_motif(genome, chrom, start1, end, strand)
            motifs[motif] = motifs.get(motif, 0) + 1
            stream.write("\t".join([circ_id(chrom, start1, end), chrom, str(start1), str(end), strand,
                                    str(end - start1 + 1), ",".join(names), motif]) + "\n")
    report["splice_motifs"] = dict(sorted(motifs.items(), key=lambda item: -item[1])[:10])
    report["canonical_gt_ag"] = motifs.get("GT-AG", 0)
    (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


FORMATS = ("bed6", "circatlas", "circbase", "circpedia")


def convert(source, out, fmt, tissue=None):
    """Database file -> BED6 (0-based start), checked against each format's convention.

    bed6       already BED6
    circatlas  circAtlas 3.0 human_bed_v3.0.txt: header, Chro Start End Strand circAltas_ID;
               Start is the 1-based first base (canonical GT-AG only under that reading)
    circbase   circBase hsa_hg19_circRNA.txt: '#' header, chrom start end strand ID ...;
               start is 0-based (BED)
    circpedia  CIRCpedia v3 Homo-sapiens.txt: circID, Location chr:start-end(strand) with a
               0-based start; `tissue` keeps rows whose long-read (TGS) evidence names it,
               e.g. "Central nervous system"
    """
    if fmt not in FORMATS:
        raise SystemExit("unknown catalogue format {}; choose from {}".format(fmt, ", ".join(FORMATS)))
    kept = skipped = 0
    location = re.compile(r"^([^:]+):(\d+)-(\d+)\(([+-])\)$")
    with _open(source) as stream, open(out, "w") as sink:
        header = None
        for line in stream:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fmt == "circpedia" and header is None:
                header = fields
                tgs = header.index("TGS")
                continue
            if fmt == "circatlas" and fields[0] == "Chro":
                continue
            if line.startswith("#"):
                continue
            if fmt == "bed6":
                chrom, start0, end, name, strand = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5]
            elif fmt == "circatlas":
                chrom, start0, end, strand, name = fields[0], int(fields[1]) - 1, int(fields[2]), fields[3], fields[4]
            elif fmt == "circbase":
                chrom, start0, end, strand, name = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]
            else:
                match = location.match(fields[1])
                if not match or (tissue and tissue not in fields[tgs]):
                    skipped += 1
                    continue
                chrom, start0, end, strand = match.group(1), int(match.group(2)), int(match.group(3)), match.group(4)
                name = fields[0]
            sink.write("\t".join([chrom, str(start0), str(end), name, "0", strand]) + "\n")
            kept += 1
    return {"format": fmt, "tissue": tissue or "", "kept": kept, "skipped": skipped}


class Chains:
    """UCSC chain file (e.g. hg19ToHg38.over.chain.gz): point lift from target to query.

    Only bases inside aligned blocks lift; a base inside blocks of two chains is
    ambiguous and does not lift (no best-chain guessing)."""

    def __init__(self, path):
        import bisect
        self._bisect = bisect
        self.blocks = {}          # tName -> sorted [(tStart, tEnd, qName, qStrand, qSize, qStart, chain_id)]
        with _open(path) as stream:
            chain = None
            for line in stream:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "chain":
                    t_name, t_start = fields[2], int(fields[5])
                    q_name, q_size, q_strand, q_start = fields[7], int(fields[8]), fields[9], int(fields[10])
                    chain = [t_name, t_start, q_name, q_size, q_strand, q_start, fields[12] if len(fields) > 12 else ""]
                    continue
                size = int(fields[0])
                t_name, t_pos, q_name, q_size, q_strand, q_pos, chain_id = chain
                self.blocks.setdefault(t_name, []).append(
                    (t_pos, t_pos + size, q_name, q_strand, q_size, q_pos, chain_id))
                if len(fields) == 3:
                    chain[1] = t_pos + size + int(fields[1])
                    chain[5] = q_pos + size + int(fields[2])
        for blocks in self.blocks.values():
            blocks.sort()
        self.starts = {name: [block[0] for block in blocks] for name, blocks in self.blocks.items()}

    def lift(self, chrom, pos0):
        """(qName, 0-based + strand position, strand_flipped) or a failure reason string."""
        blocks = self.blocks.get(chrom)
        if not blocks:
            return "contig_not_in_chain"
        hits = []
        index = self._bisect.bisect_right(self.starts[chrom], pos0) - 1
        # walk back over blocks that start earlier but may still cover pos0 (overlapping chains)
        while index >= 0 and len(hits) < 2:
            t_start, t_end, q_name, q_strand, q_size, q_start, chain_id = blocks[index]
            if t_start <= pos0 < t_end:
                q = q_start + pos0 - t_start
                hits.append((q_name, q if q_strand == "+" else q_size - 1 - q, q_strand == "-"))
            if pos0 - t_start > 50_000_000:
                break
            index -= 1
        if not hits:
            return "boundary_not_in_aligned_block"
        if len(hits) > 1:
            return "boundary_in_two_chains"
        return hits[0]


def lift(source, chain_path, out_bed, failed_path):
    """Lift a BED6 BSJ list between assemblies, boundary by boundary.

    Both the first and the last base of each circRNA must lift, to the same chromosome
    and the same orientation; the strand flips with the chain's. Rows that fail are kept
    in failed_path with their reason. The span may change (indels inside the circle);
    the change is reported, and the splice motif is checked later by `build`."""
    chains = Chains(chain_path)
    kept, failed, span_changed = 0, {}, 0
    with _open(source) as stream, open(out_bed, "w") as out, open(failed_path, "w") as bad:
        bad.write("chrom\tstart\tend\tname\tstrand\treason\n")
        for line in stream:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start0, end, name, strand = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5]
            first, last = chains.lift(chrom, start0), chains.lift(chrom, end - 1)
            reason = next((value for value in (first, last) if isinstance(value, str)), None)
            if reason is None and (first[0] != last[0] or first[2] != last[2]):
                reason = "boundaries_lift_to_different_contigs_or_orientations"
            if reason is None:
                low, high = sorted((first[1], last[1]))
                new_strand = {"+": "-", "-": "+"}[strand] if first[2] else strand
                if (first[1] > last[1]) != first[2]:
                    reason = "boundaries_swap_order"
            if reason:
                failed[reason] = failed.get(reason, 0) + 1
                bad.write("\t".join([chrom, str(start0), str(end), name, strand, reason]) + "\n")
                continue
            if high + 1 - low != end - start0:
                span_changed += 1
            out.write("\t".join([first[0], str(low), str(high + 1), name, fields[4] if len(fields) > 4 else "0",
                                 new_strand]) + "\n")
            kept += 1
    return {"lifted": kept, "failed": dict(sorted(failed.items())), "span_changed": span_changed}


def read_catalogue(path):
    """{circ_id: (chrom, start1, end, strand)} from a built catalogue.bed."""
    entries = {}
    with open(path) as stream:
        for line in stream:
            chrom, start1, end, name, _, strand = line.rstrip("\n").split("\t")[:6]
            entries[name] = (chrom, int(start1), int(end), strand)
    return entries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    b = commands.add_parser("build")
    b.add_argument("--input", required=True, nargs="+", help="BED6 files of back-spliced junctions (may be .gz)")
    b.add_argument("--genome", required=True, help="uncompressed FASTA with a samtools .fai next to it")
    b.add_argument("--out-dir", required=True)
    b.add_argument("--max-span", type=int, default=0, help="drop circRNAs spanning more bases (0 = no cap)")
    b.add_argument("--max-pseudo-bases", type=int, default=2_000_000_000)
    b.add_argument("--annotate", nargs="*", default=[], help="BED6 files that only add names to kept BSJs")
    c = commands.add_parser("convert", help="database file -> BED6 (0-based start)")
    c.add_argument("--input", required=True)
    c.add_argument("--format", required=True, choices=FORMATS)
    c.add_argument("--tissue", default="")
    c.add_argument("--out", required=True)
    c.add_argument("--report", required=True)
    l = commands.add_parser("lift", help="lift a BED6 BSJ list between assemblies (e.g. circBase hg19 -> hg38)")
    l.add_argument("--input", required=True)
    l.add_argument("--chain", required=True, help="UCSC chain file, e.g. hg19ToHg38.over.chain.gz")
    l.add_argument("--out", required=True)
    l.add_argument("--failed", required=True)
    l.add_argument("--report", required=True)
    args = parser.parse_args(argv)
    if args.command == "convert":
        report = convert(args.input, args.out, args.format, args.tissue or None)
        Path(args.report).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(json.dumps(report))
        return
    if args.command == "lift":
        report = lift(args.input, args.chain, args.out, args.failed)
        Path(args.report).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(json.dumps(report))
        return
    report = build(args.input, args.genome, args.out_dir, args.max_span, args.max_pseudo_bases,
                   annotate=args.annotate)
    print(json.dumps({key: report[key] for key in ("catalogue_id", "kept", "dropped", "pseudo_reference_bases")}))


if __name__ == "__main__":
    sys.exit(main())

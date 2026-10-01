"""3' UTR resources for the APA modules, derived from the configured GTF.

qapa-db       QAPA's gene identifier table (Gene stable ID, Transcript stable
              ID, Gene type, Transcript type, Gene name), with Ensembl version
              suffixes removed as QAPA expects. QAPA keeps protein-coding genes
              and transcripts only.
dapars-utr    DaPars2's annotated 3' UTRs (BED6, name "transcript|gene|chrom|
              strand"): the last exon's part downstream of the coding end, for
              protein-coding transcripts whose coding end lies in the last exon
              (3' UTRs containing introns are skipped, as in QAPA). Identical
              regions are kept once. Also writes the merged, unstranded UTR
              windows used to keep per-run coverage small.
qapa-gtf      The GTF lines QAPA's genePred is made from: every feature of the
              selected transcripts (all, or only those tagged `basic`, as in
              QAPA's own GENCODE basic input).
polya-bed     GENCODE's poly(A) feature GTF as the BED6 QAPA's -g expects
              (name `polyA_site`, other features dropped).
check-bed     Fails when a built 3' UTR library has fewer than --min entries.

Coordinates in BED are 0-based, half open. Transcripts without a CDS, and
synthetic models without gene_type/transcript_type, are skipped and counted
in the report. Both APA references need coding ends: an annotation with no
CDS features at all (an exon-only build) is refused instead of silently
giving an empty or tiny library.
"""

import argparse
import gzip
import json
import re
from collections import defaultdict


VERSION = re.compile(r"\.\d+(_PAR_Y)?$")
TAG = re.compile(r'\btag "([^"]*)"')
TRANSCRIPT_SETS = ("basic", "all")
NO_CDS = ("the APA annotation has no CDS features, so 3' UTRs cannot be located (an exon-only "
          "GTF?); set apa_annotation_gtf to a GTF with CDS lines, e.g. GENCODE comprehensive")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _attributes(text):
    return dict(re.findall(r'(\S+) "([^"]*)"', text))


def selected(tags, transcripts="basic"):
    if transcripts not in TRANSCRIPT_SETS:
        raise ValueError("apa_transcripts must be basic or all, not {}".format(transcripts))
    return transcripts == "all" or "basic" in tags


def read_transcripts(gtf):
    transcripts = {}
    with _open(gtf) as stream:
        for line in stream:
            if not line.startswith("#"):
                _add_feature(transcripts, line.rstrip("\n").split("\t"))
    return transcripts


def _add_feature(transcripts, fields):
    if len(fields) < 9 or fields[2] not in ("transcript", "exon", "CDS", "stop_codon"):
        return
    attrs = _attributes(fields[8])
    tid = attrs.get("transcript_id")
    if not tid:
        return
    record = transcripts.setdefault(tid, {
        "chrom": fields[0], "strand": fields[6], "exons": [], "cds": [],
        "gene_id": attrs.get("gene_id", ""), "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")),
        "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
        "transcript_type": attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
        "tags": set()})
    record["tags"].update(TAG.findall(fields[8]))
    start, end = int(fields[3]) - 1, int(fields[4])
    if fields[2] == "exon":
        record["exons"].append((start, end))
    elif fields[2] in ("CDS", "stop_codon"):
        record["cds"].append((start, end))


def strip_version(identifier):
    return VERSION.sub("", identifier)


def qapa_db(transcripts, out):
    rows = set()
    for tid, record in transcripts.items():
        if not record["gene_type"] or not record["transcript_type"]:
            continue
        rows.add((strip_version(record["gene_id"]), strip_version(tid), record["gene_type"],
                  record["transcript_type"], record["gene_name"]))
    with open(out, "w") as stream:
        stream.write("Gene stable ID\tTranscript stable ID\tGene type\tTranscript type\tGene name\n")
        for row in sorted(rows, key=lambda row: (row[0], row[1])):
            stream.write("\t".join(row) + "\n")
    return len(rows)


def three_prime_utr(record):
    """(start, end) of the last exon's part past the coding end, or a skip reason."""
    if not record["cds"]:
        return "no_cds"
    exons = sorted(record["exons"])
    if not exons:
        return "no_exons"
    if record["strand"] == "+":
        coding_end = max(end for _, end in record["cds"])
        last = exons[-1]
        if coding_end < last[0]:
            return "utr_with_intron"
        if coding_end >= last[1]:
            return "no_utr"
        return (coding_end, last[1])
    coding_start = min(start for start, _ in record["cds"])
    last = exons[0]
    if coding_start > last[1]:
        return "utr_with_intron"
    if coding_start <= last[0]:
        return "no_utr"
    return (last[0], coding_start)


def coding_transcripts(transcripts):
    """Protein-coding transcripts that carry CDS features."""
    return sum(1 for record in transcripts.values()
               if record["transcript_type"] == "protein_coding" and record["cds"])


def dapars_utr(transcripts, out_bed, out_windows, transcript_set="all", min_utrs=0):
    if not coding_transcripts(transcripts):
        raise ValueError(NO_CDS)
    skipped = defaultdict(int)
    seen = {}
    for tid in sorted(transcripts):
        record = transcripts[tid]
        if record["gene_type"] != "protein_coding" or record["transcript_type"] != "protein_coding":
            skipped["not_protein_coding"] += 1
            continue
        if not selected(record["tags"], transcript_set):
            skipped["not_basic"] += 1
            continue
        utr = three_prime_utr(record)
        if isinstance(utr, str):
            skipped[utr] += 1
            continue
        key = (record["chrom"], utr[0], utr[1], record["strand"])
        if key in seen:
            skipped["duplicate_region"] += 1
            continue
        seen[key] = "{}|{}|{}|{}".format(strip_version(tid), record["gene_name"],
                                        record["chrom"], record["strand"])
    rows = sorted(seen.items())
    if len(rows) < min_utrs:
        raise ValueError("only {} DaPars2 3' UTRs (minimum apa_min_utrs = {}); skipped: {}".format(
            len(rows), min_utrs, dict(skipped)))
    with open(out_bed, "w") as stream:
        for (chrom, start, end, strand), name in rows:
            stream.write("{}\t{}\t{}\t{}\t0\t{}\n".format(chrom, start, end, name, strand))
    merged = []
    for (chrom, start, end, _), _name in sorted(rows, key=lambda item: item[0][:3]):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    with open(out_windows, "w") as stream:
        for chrom, start, end in merged:
            stream.write("{}\t{}\t{}\n".format(chrom, start, end))
    return len(rows), dict(skipped)


def qapa_gtf(gtf, out, transcript_set="basic"):
    """Write the selected transcripts' feature lines; refuse an annotation without CDS."""
    kept, cds = set(), 0
    with _open(gtf) as source, open(out, "w") as target:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or 'transcript_id "' not in fields[8]:
                continue
            if not selected(TAG.findall(fields[8]), transcript_set):
                continue
            target.write(line)
            kept.add(_attributes(fields[8])["transcript_id"])
            cds += fields[2] == "CDS"
    if not cds:
        raise ValueError(NO_CDS)
    return len(kept)


def polya_bed(gtf, out):
    """GENCODE poly(A) features to the BED6 QAPA's -g reads (polyA_site only)."""
    rows = 0
    with _open(gtf) as source, open(out, "w") as target:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) < 9 or fields[2] != "polyA_site":
                continue
            target.write("{}\t{}\t{}\tpolyA_site\t0\t{}\n".format(
                fields[0], int(fields[3]) - 1, fields[4], fields[6]))
            rows += 1
    if not rows:
        raise ValueError("{} has no polyA_site features".format(gtf))
    return rows


def check_bed(bed, minimum, label):
    with open(bed) as stream:
        rows = sum(1 for line in stream if line.strip() and not line.startswith(("#", "track")))
    if rows < minimum:
        raise ValueError("{} has only {} entries (minimum apa_min_utrs = {}): check the APA "
                         "annotation (CDS features, apa_transcripts)".format(label, rows, minimum))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    db = commands.add_parser("qapa-db")
    db.add_argument("--gtf", required=True)
    db.add_argument("--output", required=True)
    utr = commands.add_parser("dapars-utr")
    utr.add_argument("--gtf", required=True)
    utr.add_argument("--bed", required=True)
    utr.add_argument("--windows", required=True)
    utr.add_argument("--report", required=True)
    utr.add_argument("--transcripts", choices=TRANSCRIPT_SETS, default="all")
    utr.add_argument("--min-utrs", type=int, default=0)
    genes = commands.add_parser("qapa-gtf")
    genes.add_argument("--gtf", required=True)
    genes.add_argument("--output", required=True)
    genes.add_argument("--transcripts", choices=TRANSCRIPT_SETS, default="basic")
    polya = commands.add_parser("polya-bed")
    polya.add_argument("--gtf", required=True)
    polya.add_argument("--output", required=True)
    check = commands.add_parser("check-bed")
    check.add_argument("--bed", required=True)
    check.add_argument("--min", type=int, required=True)
    check.add_argument("--label", default="3' UTR library")
    args = parser.parse_args(argv)
    try:
        if args.command == "qapa-gtf":
            qapa_gtf(args.gtf, args.output, args.transcripts)
        elif args.command == "polya-bed":
            polya_bed(args.gtf, args.output)
        elif args.command == "check-bed":
            check_bed(args.bed, args.min, args.label)
        elif args.command == "qapa-db":
            qapa_db(read_transcripts(args.gtf), args.output)
        else:
            transcripts = read_transcripts(args.gtf)
            kept, skipped = dapars_utr(transcripts, args.bed, args.windows, args.transcripts, args.min_utrs)
            with open(args.report, "w") as stream:
                json.dump({"utrs": kept, "skipped": skipped, "transcripts": len(transcripts),
                           "transcript_set": args.transcripts}, stream, sort_keys=True, indent=2)
                stream.write("\n")
    except ValueError as error:
        parser.exit(1, "apa_annotation.py {}: error: {}\n".format(args.command, error))


if __name__ == "__main__":
    main()

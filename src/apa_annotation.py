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

Coordinates in BED are 0-based, half open. Transcripts without a CDS, and
synthetic models without gene_type/transcript_type, are skipped and counted
in the report.
"""

import argparse
import gzip
import json
import re
from collections import defaultdict


VERSION = re.compile(r"\.\d+(_PAR_Y)?$")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _attributes(text):
    return dict(re.findall(r'(\S+) "([^"]*)"', text))


def read_transcripts(gtf):
    transcripts = {}
    for line in _open(gtf):
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] not in ("transcript", "exon", "CDS", "stop_codon"):
            continue
        attrs = _attributes(fields[8])
        tid = attrs.get("transcript_id")
        if not tid:
            continue
        record = transcripts.setdefault(tid, {
            "chrom": fields[0], "strand": fields[6], "exons": [], "cds": [],
            "gene_id": attrs.get("gene_id", ""), "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")),
            "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
            "transcript_type": attrs.get("transcript_type", attrs.get("transcript_biotype", ""))})
        start, end = int(fields[3]) - 1, int(fields[4])
        if fields[2] == "exon":
            record["exons"].append((start, end))
        elif fields[2] in ("CDS", "stop_codon"):
            record["cds"].append((start, end))
    return transcripts


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


def dapars_utr(transcripts, out_bed, out_windows):
    skipped = defaultdict(int)
    seen = {}
    for tid in sorted(transcripts):
        record = transcripts[tid]
        if record["gene_type"] != "protein_coding" or record["transcript_type"] != "protein_coding":
            skipped["not_protein_coding"] += 1
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
    args = parser.parse_args(argv)
    transcripts = read_transcripts(args.gtf)
    if args.command == "qapa-db":
        qapa_db(transcripts, args.output)
    else:
        kept, skipped = dapars_utr(transcripts, args.bed, args.windows)
        with open(args.report, "w") as stream:
            json.dump({"utrs": kept, "skipped": skipped, "transcripts": len(transcripts)},
                      stream, sort_keys=True, indent=2)
            stream.write("\n")


if __name__ == "__main__":
    main()

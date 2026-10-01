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
apa-set       The MegaSearch APA set (apa_transcripts: megasearch): every line
              of the GENCODE transcripts that are protein-coding with CDS and
              their coding end in the last exon; free of the partial-model tags
              cds_start_NF, cds_end_NF, mRNA_start_NF, mRNA_end_NF; not a
              readthrough; and whose 3' end is supported: a Tier-1 model (MANE,
              Ensembl canonical, GENCODE Primary, or basic with TSL 1-2), or
              within --slop nt (50) of a PolyASite 2.0 cluster supported by
              >= 2 protocols or of a GENCODE polyA_site. The same evidence and
              tiers as the v2 Whippet annotation's end rule (megasearch-whippet-
              annotation-code, scripts/apa_transcript_set.py, which must select
              the same transcripts). Written as a reproducible gzip.
check-bed     Fails when a built 3' UTR library has fewer than --min entries.

Coordinates in BED are 0-based, half open. Transcripts without a CDS, those
whose 3' end is incomplete (cds_end_NF, mRNA_end_NF), and synthetic models
without gene_type/transcript_type, are skipped by dapars-utr and qapa-gtf and
counted in the report. Both APA references need coding ends: an annotation with no
CDS features at all (an exon-only build) is refused instead of silently
giving an empty or tiny library.
"""

import argparse
import bisect
import gzip
import json
import re
from collections import Counter, defaultdict


VERSION = re.compile(r"\.\d+(_PAR_Y)?$")
TAG = re.compile(r'\btag "([^"]*)"')
TRANSCRIPT_SETS = ("basic", "all")
# 3' end incomplete: no reliable 3' UTR (dropped by dapars-utr and qapa-gtf for any input)
INCOMPLETE_3PRIME = {"cds_end_NF", "mRNA_end_NF"}
PARTIAL = {"cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"}
TIER1_TAGS = {"MANE_Select", "MANE_Plus_Clinical", "Ensembl_canonical", "GENCODE_Primary"}
APA_SET_RULES = ("protein_coding", "has_cds", "coding_end_in_last_exon", "complete_model",
                 "supported_3prime_end", "not_readthrough")
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
        if record["tags"] & INCOMPLETE_3PRIME:
            skipped["incomplete_3prime"] += 1
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
    """Write the selected transcripts' feature lines; refuse an annotation without CDS.

    Transcripts without CDS lines, or with an incomplete 3' end, are left out."""
    coding = set()
    with _open(gtf) as source:
        for line in source:
            fields = line.split("\t")
            if not line.startswith("#") and len(fields) >= 9 and fields[2] == "CDS":
                coding.add(_attributes(fields[8]).get("transcript_id"))
    kept, cds = set(), 0
    with _open(gtf) as source, open(out, "w") as target:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or 'transcript_id "' not in fields[8]:
                continue
            tags = TAG.findall(fields[8])
            if not selected(tags, transcript_set) or INCOMPLETE_3PRIME.intersection(tags):
                continue
            if _attributes(fields[8])["transcript_id"] not in coding:
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


# ------------------------------------------------------------ MegaSearch APA set

def load_polya_evidence(gencode_polya, polyasite):
    """{(chrom, strand): sorted [(start, end)]}, 1-based closed: GENCODE polyA_site
    features and PolyASite 2.0 clusters supported by >= 2 protocols."""
    evidence = defaultdict(list)
    with _open(gencode_polya) as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) < 9 or fields[2] != "polyA_site":
                continue
            evidence[(fields[0], fields[6])].append((int(fields[3]), int(fields[4])))
    with _open(polyasite) as stream:
        for line in stream:
            if line.startswith(("#", "chromosome")):
                continue
            fields = line.rstrip("\n").split("\t")
            if int(fields[7]) < 2:
                continue
            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
            evidence[(chrom, fields[5])].append((int(fields[1]) + 1, int(fields[2])))
    for key in evidence:
        evidence[key].sort()
    return evidence


def near(evidence, chrom, strand, position, slop):
    """True when position lies within slop nt of an evidence interval."""
    intervals = evidence.get((chrom, strand))
    if not intervals:
        return False
    low = bisect.bisect_left(intervals, (position - slop - 100000, 0))
    high = bisect.bisect_right(intervals, (position + slop, 10 ** 12))
    return any(start - slop <= position <= end + slop for start, end in intervals[max(0, low):high])


def tier1(model):
    return bool(model["tags"] & TIER1_TAGS) or ("basic" in model["tags"] and model["tsl"] in ("1", "2"))


def apa_set_models(gtf):
    """{transcript_id: model}, 1-based closed exons and coding bases (CDS, stop_codon)."""
    models = {}
    with _open(gtf) as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon", "CDS", "stop_codon"):
                continue
            attrs = _attributes(fields[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            if fields[2] == "transcript":
                models[tid] = {"gene_id": attrs.get("gene_id", ""), "chrom": fields[0], "strand": fields[6],
                               "type": attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
                               "tsl": attrs.get("transcript_support_level", "missing"),
                               "tags": set(TAG.findall(fields[8])), "exons": [], "coding": []}
            elif tid in models:
                key = "exons" if fields[2] == "exon" else "coding"
                models[tid][key].append((int(fields[3]), int(fields[4])))
    for model in models.values():
        model["exons"].sort()
    return models


def apa_set_failures(model, evidence, slop=50):
    """The APA-set rules (APA_SET_RULES order) this transcript fails."""
    failed = []
    plus = model["strand"] == "+"
    if model["type"] != "protein_coding":
        failed.append("protein_coding")
    if not model["coding"]:
        failed.append("has_cds")
    else:
        coding_end = max(e for _, e in model["coding"]) if plus else min(s for s, _ in model["coding"])
        last = model["exons"][-1] if plus else model["exons"][0]
        if not last[0] <= coding_end <= last[1]:
            failed.append("coding_end_in_last_exon")
    if model["tags"] & PARTIAL:
        failed.append("complete_model")
    end = model["exons"][-1][1] if plus else model["exons"][0][0]
    if not (tier1(model) or near(evidence, model["chrom"], model["strand"], end, slop)):
        failed.append("supported_3prime_end")
    if "readthrough_transcript" in model["tags"]:
        failed.append("not_readthrough")
    return failed


def apa_set(gtf, gencode_polya, polyasite, out, slop=50):
    """Write the MegaSearch APA set GTF (reproducible gzip); return its report."""
    evidence = load_polya_evidence(gencode_polya, polyasite)
    models = apa_set_models(gtf)
    if not any(model["coding"] for model in models.values()):
        raise ValueError(NO_CDS)
    first, every = Counter(), Counter()
    chosen, genes = set(), set()
    for tid in sorted(models):
        failed = apa_set_failures(models[tid], evidence, slop)
        every.update(failed)
        if failed:
            first[failed[0]] += 1
        else:
            chosen.add(tid)
            genes.add(models[tid]["gene_id"])
    lines = 0
    with open(out, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as target, \
            _open(gtf) as source:
        target.write(b"##description: MegaSearch APA set (apa_annotation.py apa-set)\n")
        for line in source:
            if line.startswith("#"):
                target.write(line.encode())
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9:
                continue
            attrs = _attributes(fields[8])
            keep = attrs.get("gene_id") in genes if fields[2] == "gene" else attrs.get("transcript_id") in chosen
            if keep:
                target.write(line.encode())
                lines += 1
    basic = {tid for tid, model in models.items()
             if model["type"] == "protein_coding" and "basic" in model["tags"]}
    return {"transcripts_in_gtf": len(models), "selected_transcripts": len(chosen),
            "selected_genes": len(genes), "lines_written": lines, "slop": slop,
            "first_failed_rule": {rule: first[rule] for rule in APA_SET_RULES},
            "failing_each_rule": {rule: every[rule] for rule in APA_SET_RULES},
            "selected_tier1": sum(tier1(models[tid]) for tid in chosen),
            "basic_protein_coding": len(basic), "basic_protein_coding_not_selected": len(basic - chosen)}


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
    aset = commands.add_parser("apa-set")
    aset.add_argument("--gtf", required=True)
    aset.add_argument("--gencode-polya", required=True)
    aset.add_argument("--polyasite", required=True)
    aset.add_argument("--slop", type=int, default=50)
    aset.add_argument("--out", required=True)
    aset.add_argument("--report", required=True)
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
        elif args.command == "apa-set":
            report = apa_set(args.gtf, args.gencode_polya, args.polyasite, args.out, args.slop)
            with open(args.report, "w") as stream:
                json.dump(dict(report, inputs={"gtf": args.gtf, "gencode_polya": args.gencode_polya,
                                               "polyasite": args.polyasite}),
                          stream, sort_keys=True, indent=2)
                stream.write("\n")
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

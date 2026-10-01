"""QAPA module: per-run 3' UTR isoform quantification caches and PAU comparisons.

keep       After `salmon quant` against QAPA's UTR library: keep quant.sf
           (gzipped; Name, Length, EffectiveLength, TPM, NumReads),
           lib_format_counts.json and meta_info.json. Enough to recompute PAU
           and usage for any grouping without the reads.
pau        Per comparison, from the caches: one quant.sf per biological
           replicate (technical runs: NumReads summed, effective length
           read-weighted, TPM recomputed), `qapa quant` for PAU, and a count
           table per 3' UTR isoform (one row of QAPA's PAU table, keyed by its
           first transcript and UTR3 coordinates, grouped by the gene QAPA's
           identifier table gives that transcript) for the DEXSeq usage test.
normalize  Joins QAPA's PAU table with the DEXSeq results into the shared
           normalized columns: effect = mean PAU(A) - mean PAU(B) for that
           site (percent, as QAPA reports it); p and q from DEXSeq.
"""

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path

if __package__:
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from umbrella_modules import NORMALIZED, write_status

# QAPA's library names, as its create_merged_data.R splits them:
# TX_GENE[,TX_GENE...]_SPECIES_CHR_LASTEXONSTART_LASTEXONEND_STRAND_utr_START_END,
# with "::chr:start-end(strand)" appended by `qapa fasta` (bedtools getfasta)
QAPA_NAME = re.compile(r"^(?P<ids>[^_]+_[^_,]+(?:,[^_]+_[^_,]+)*)_[^_]+_[^_]+_\d+_\d+_[-+]_utr_"
                       r"(?P<start>\d+)_(?P<end>\d+)$")
VERSION = re.compile(r"\.\d+(_PAR_Y)?$")
PAU_EFFECT = "mean PAU(A) - mean PAU(B) of this poly(A) site within its gene, in percent"


def keep(salmon_dir, out_dir, run_id, extra=None):
    salmon_dir, out_dir = Path(salmon_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(salmon_dir / "quant.sf", "rb") as source, gzip.open(out_dir / "quant.sf.gz", "wb") as target:
        shutil.copyfileobj(source, target)
    shutil.copyfile(salmon_dir / "lib_format_counts.json", out_dir / "lib_format_counts.json")
    shutil.copyfile(salmon_dir / "aux_info" / "meta_info.json", out_dir / "meta_info.json")
    (out_dir / "cache.json").write_text(json.dumps(dict(extra or {}, run_id=run_id, tool="qapa"),
                                                   sort_keys=True, indent=2) + "\n")


def read_quant(path):
    rows = OrderedDict()
    with gzip.open(path, "rt") as stream:
        header = next(stream).rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            rows[fields[index["Name"]]] = (int(fields[index["Length"]]),
                                           float(fields[index["EffectiveLength"]]),
                                           float(fields[index["NumReads"]]))
    return rows


def combine_runs(paths):
    """One quant.sf's rows from several technical runs of one replicate."""
    tables = [read_quant(path) for path in paths]
    names = list(tables[0])
    if any(list(table) != names for table in tables[1:]):
        raise ValueError("technical runs were quantified against different UTR libraries")
    combined = OrderedDict()
    for name in names:
        reads = sum(table[name][2] for table in tables)
        if reads > 0:
            effective = sum(table[name][1] * table[name][2] for table in tables) / reads
        else:
            effective = sum(table[name][1] for table in tables) / len(tables)
        combined[name] = (tables[0][name][0], effective, reads)
    rate = {name: (reads / effective if effective > 0 else 0.0)
            for name, (_, effective, reads) in combined.items()}
    total = sum(rate.values()) or 1.0
    return [(name, length, effective, rate[name] / total * 1e6, reads)
            for name, (length, effective, reads) in combined.items()]


def site_key(name):
    """(first transcript, UTR3.Start, UTR3.End) of a library sequence: one PAU table row."""
    match = QAPA_NAME.match(name.split("::", 1)[0])
    if not match:
        raise ValueError("unexpected QAPA sequence name: {}".format(name))
    return VERSION.sub("", match.group("ids").split("_", 1)[0]), match.group("start"), match.group("end")


def pau_site(row):
    """The site_key of a PAU table row (QAPA keeps the name's transcripts and UTR3 ends)."""
    return VERSION.sub("", row["Transcript"].split(",")[0]), _int(row["UTR3.Start"]), _int(row["UTR3.End"])


def gene_of_transcript(db):
    """Transcript -> gene from QAPA's identifier table (protein-coding genes, as QAPA keeps)."""
    genes = {}
    with open(db) as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row.get("Gene type") == "protein_coding":
                genes[row["Transcript stable ID"]] = row["Gene stable ID"]
    return genes


def pau(preflight_path, caches, db, out_dir):
    preflight = json.loads(Path(preflight_path).read_text())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    replicates = preflight["replicates"]["a"] + preflight["replicates"]["b"]
    members = defaultdict(list)
    for run, replicate in preflight["collapse"].items():
        members[replicate].append(run)
    counts = defaultdict(lambda: defaultdict(float))
    genes = gene_of_transcript(db)
    with tempfile.TemporaryDirectory(prefix="qapa-", dir=str(out_dir)) as directory:
        work = Path(directory)
        paths = []
        for replicate in replicates:
            rows = combine_runs([Path(caches[run]) / "quant.sf.gz" for run in sorted(members[replicate])])
            (work / replicate).mkdir()
            with open(work / replicate / "quant.sf", "w") as stream:
                stream.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
                for name, length, effective, tpm, reads in rows:
                    stream.write("{}\t{}\t{:.3f}\t{:.6f}\t{:.3f}\n".format(name, length, effective, tpm, reads))
                    site = site_key(name)
                    if site[0] in genes:        # QAPA's PAU table drops the others too
                        counts[site][replicate] += reads
            # absolute: qapa runs inside the private folder, the workflow passes relative paths
            paths.append(str((work / replicate / "quant.sf").resolve()))
        with open(out_dir / "pau.tsv", "w") as stream, open(out_dir / "qapa.log", "w") as log:
            subprocess.run(["qapa", "quant", "--db", str(Path(db).resolve())] + paths, stdout=stream,
                           stderr=log, check=True, cwd=str(work))
    with gzip.open(out_dir / "site_counts.tsv.gz", "wt") as stream:
        stream.write("site_id\tgene_id\t" + "\t".join(replicates) + "\n")
        for site, by_replicate in sorted(counts.items()):
            values = "\t".join("{:.0f}".format(by_replicate[r]) for r in replicates)
            stream.write("{}\t{}\t{}\n".format("_".join(site), genes[site[0]], values))
    with open(out_dir / "samples.tsv", "w") as stream:
        stream.write("sample\tcondition\n")
        for side in ("a", "b"):
            for replicate in preflight["replicates"][side]:
                stream.write("{}\t{}\n".format(replicate, side))


def normalize(preflight_path, pau_table, dexseq_table, out_dir, ids=None):
    preflight = json.loads(Path(preflight_path).read_text())
    ids = ids or {}
    out_dir = Path(out_dir)
    base = dict(project_id=preflight["project_id"], comparison_id=preflight["comparison_id"], **ids)
    dexseq = {}
    if Path(dexseq_table).exists():
        with gzip.open(dexseq_table, "rt") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                dexseq[row["site_id"]] = row
    records = []
    with open(pau_table) as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            site = "_".join(pau_site(row))
            values = {side: [float(row[r + ".PAU"]) for r in preflight["replicates"][side]
                             if row.get(r + ".PAU") not in (None, "", "NA")] for side in ("a", "b")}
            effect = (sum(values["a"]) / len(values["a"]) - sum(values["b"]) / len(values["b"])
                      if values["a"] and values["b"] else None)
            test = dexseq.get(site, {})
            p, q = test.get("pvalue", ""), test.get("padj", "")
            p = "" if p in ("NA", None) else p
            q = "" if q in ("NA", None) else q
            status = ("single_site" if row.get("Num_Events") in ("1", "1.0")
                      else ("tested" if p else "untested"))
            records.append(dict(base, tool="qapa", tool_version=ids.get("tool_version", ""),
                                feature_id=row["APA_ID"], gene_id=row["Gene"],
                                gene_name=row.get("Gene_Name", ""), chrom=row["Chr"],
                                start=_int(row["UTR3.Start"]), end=_int(row["UTR3.End"]),
                                coordinates="QAPA UTR3.Start/UTR3.End; last exon {}-{}".format(
                                    row.get("LastExon.Start", ""), row.get("LastExon.End", "")),
                                strand=row["Strand"], event_class="3utr_polya_site",
                                effect="" if effect is None else "{:.6g}".format(effect),
                                effect_definition=PAU_EFFECT, n_a=len(values["a"]),
                                n_b=len(values["b"]), coverage="",
                                native_statistic=test.get("log2fold_a_b", ""),
                                statistic_type="DEXSeq site-usage log2 fold change (A/B)",
                                p_value=p, q_value=q, status=status))
    with gzip.open(out_dir / "results.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    tested = sum(bool(r["p_value"]) for r in records)
    status, reasons = "ok", []
    if not records:
        status, reasons = "empty", ["QAPA's PAU table has no sites"]
    elif not tested:
        status, reasons = "no_tests", ["no site has a DEXSeq p-value (genes with >= 2 sites, "
                                       ">= 2 replicates per side)"]
    write_status(out_dir / "status.json", "qapa", status, reasons, sites=len(records),
                 dexseq_rows=len(dexseq), tested=tested, **base)


def _int(value):
    try:
        return str(int(float(value)))
    except ValueError:
        return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    k = commands.add_parser("keep")
    k.add_argument("--salmon", required=True)
    k.add_argument("--out-dir", required=True)
    k.add_argument("--run-id", required=True)
    k.add_argument("--identity", default="{}")
    p = commands.add_parser("pau")
    p.add_argument("--preflight", required=True)
    p.add_argument("--cache", action="append", default=[], help="RUN=DIR")
    p.add_argument("--db", required=True)
    p.add_argument("--out-dir", required=True)
    n = commands.add_parser("normalize")
    n.add_argument("--preflight", required=True)
    n.add_argument("--pau", required=True)
    n.add_argument("--dexseq", required=True)
    n.add_argument("--out-dir", required=True)
    n.add_argument("--ids", default="{}")
    args = parser.parse_args(argv)
    if args.command == "keep":
        keep(args.salmon, args.out_dir, args.run_id, json.loads(args.identity))
    elif args.command == "pau":
        pau(args.preflight, dict(item.split("=", 1) for item in args.cache), args.db, args.out_dir)
    else:
        normalize(args.preflight, args.pau, args.dexseq, args.out_dir, json.loads(args.ids))


if __name__ == "__main__":
    main()

"""DaPars2 module: per-run 3' UTR coverage caches and per-comparison PDUI analysis.

coverage  Per run, from the temporary HISAT2 BAM: raw read coverage
          (`bedtools genomecov -split -bg`, spliced reads split, no scaling)
          clipped to the merged 3' UTR windows, gzipped bedGraph; intervals
          absent from it have zero coverage (DaPars2 reads gaps as zero).
          Also the primary mapped read count, DaPars2's depth input.
          Coverage is unstranded, as DaPars2 expects; loci where genes on
          opposite strands overlap are flagged in the results.
compare   Per comparison, from the caches only: sums technical runs into one
          bedGraph per biological replicate, runs DaPars2 by chromosome in a
          private folder, and writes
            native.tsv.gz   DaPars2's rows (Gene, fit_value, Predicted_Proximal_APA,
                            Loci, <replicate>_PDUI ...)
            results.tsv.gz  normalized rows: effect = mean PDUI(A) - mean PDUI(B);
                            two-sided Welch t-test where each side has >= 2 values
                            and some variance, BH over tested regions (exploratory:
                            DaPars2 itself reports no p-values)
            status.json
"""

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

if __package__:
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from umbrella_modules import NORMALIZED, write_status

# DaPars2 starts worker processes at module level; it needs the "fork" start
# method (the Linux default up to Python 3.13, not the macOS or 3.14+ default)
FORK_RUNNER = ("import multiprocessing, runpy, sys; multiprocessing.set_start_method('fork'); "
               "script = sys.argv[1]; sys.argv = sys.argv[1:]; runpy.run_path(script, run_name='__main__')")
PDUI_EFFECT = "mean PDUI(A) - mean PDUI(B); PDUI = distal poly(A) usage of the 3' UTR"


# ------------------------------------------------------------ per run

def coverage(bam, windows, out_dir, run_id, extra=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bedgraph = out_dir / "utr_coverage.bedgraph.gz"
    temporary = out_dir / "utr_coverage.bedgraph.gz.tmp"
    genomecov = subprocess.Popen(["bedtools", "genomecov", "-ibam", str(bam), "-bg", "-split"],
                                 stdout=subprocess.PIPE)
    clip = subprocess.Popen(["bedtools", "intersect", "-a", "stdin", "-b", str(windows)],
                            stdin=genomecov.stdout, stdout=subprocess.PIPE)
    genomecov.stdout.close()
    with gzip.open(temporary, "wb", compresslevel=6) as stream:
        shutil.copyfileobj(clip.stdout, stream)
    if clip.wait() or genomecov.wait():
        raise RuntimeError("coverage extraction failed for {}".format(run_id))
    os.replace(temporary, bedgraph)
    mapped = int(subprocess.run(["samtools", "view", "-c", "-F", "0x904", str(bam)],
                                check=True, capture_output=True, text=True).stdout.strip())
    (out_dir / "depth.json").write_text(json.dumps({"run_id": run_id, "mapped_reads": mapped,
                                                    "definition": "primary mapped reads (-F 0x904)"},
                                                   sort_keys=True) + "\n")
    (out_dir / "cache.json").write_text(json.dumps(dict(extra or {}, run_id=run_id, tool="dapars2"),
                                                   sort_keys=True, indent=2) + "\n")


# ------------------------------------------------------------ per comparison

def _intervals(path):
    with gzip.open(path, "rt") as stream:
        for line in stream:
            chrom, start, end, value = line.split("\t")[:4]
            yield chrom, int(start), int(end), float(value)


def sum_bedgraphs(paths, out, chrom_order):
    """Sum bedGraphs (sorted per chromosome) into one; gaps are zero."""
    rank = {chrom: index for index, chrom in enumerate(chrom_order)}
    by_chrom = {}
    for path in paths:
        for chrom, start, end, value in _intervals(path):
            by_chrom.setdefault(chrom, []).append((start, end, value))
    with open(out, "w") as stream:
        for chrom in sorted(by_chrom, key=lambda c: (rank.get(c, len(rank)), c)):
            events = []
            for start, end, value in by_chrom[chrom]:
                events.append((start, value))
                events.append((end, -value))
            events.sort()
            level, previous = 0.0, None
            for position, delta in events:
                if previous is not None and position > previous and abs(level) > 1e-9:
                    stream.write("{}\t{}\t{}\t{:g}\n".format(chrom, previous, position, level))
                level += delta
                previous = position


def welch(a, b):
    """Two-sided Welch t-test p-value, or None when not testable."""
    if len(a) < 2 or len(b) < 2:
        return None
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    se2 = va / len(a) + vb / len(b)
    if se2 <= 0:
        return None
    t = (ma - mb) / math.sqrt(se2)
    df = se2 ** 2 / ((va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
    return _incomplete_beta(df / 2.0, 0.5, df / (df + t * t))


def _incomplete_beta(a, b, x):
    """Regularized incomplete beta I_x(a, b) (continued fraction, Numerical Recipes).
    With a = df/2, b = 1/2 and x = df/(df + t^2) it is the two-sided t p-value."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                     + a * math.log(x) + b * math.log(1 - x))
    if x > (a + 1) / (a + b + 2):
        return 1.0 - _incomplete_beta(b, a, 1 - x)
    c, d = 1.0, 1.0 - (a + b) * x / (a + 1)
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    f = d
    for m in range(1, 300):
        for numerator in (m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)),
                          -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))):
            d = 1.0 + numerator * d
            d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
            c = 1.0 + numerator / c
            c = c if abs(c) > 1e-300 else 1e-300
            f *= c * d
        if abs(c * d - 1.0) < 1e-12:
            break
    return front * f / a


def benjamini_hochberg(pvalues):
    indexed = sorted((p, i) for i, p in enumerate(pvalues) if p is not None)
    q = [None] * len(pvalues)
    running = 1.0
    for rank in range(len(indexed), 0, -1):
        p, index = indexed[rank - 1]
        running = min(running, p * len(indexed) / rank)
        q[index] = running
    return q


def overlapping_opposite_strands(utr_bed):
    rows = []
    with open(utr_bed) as stream:
        for line in stream:
            chrom, start, end, name, _, strand = line.rstrip("\n").split("\t")[:6]
            rows.append((chrom, int(start), int(end), strand, "{}:{}-{}".format(chrom, start, end)))
    rows.sort()
    flagged, active = set(), []
    for chrom, start, end, strand, loci in rows:
        active = [item for item in active if item[0] == chrom and item[2] > start]
        for other in active:
            if other[3] != strand:
                flagged.update((loci, other[4]))
        active.append((chrom, start, end, strand, loci))
    return flagged


def compare(preflight_path, caches, utr_bed, dapars2_dir, out_dir, threads=1,
            coverage_threshold=10, ids=None):
    preflight = json.loads(Path(preflight_path).read_text())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ids = ids or {}
    status_path = out_dir / "status.json"
    base = dict(project_id=preflight["project_id"], comparison_id=preflight["comparison_id"], **ids)
    if not preflight["inference_supported"]:
        _empty_outputs(out_dir)
        write_status(status_path, "dapars2", "unsupported", preflight["reasons"], **base)
        return
    chroms = []
    with open(utr_bed) as stream:
        for line in stream:
            chrom = line.split("\t", 1)[0]
            if chrom not in chroms:
                chroms.append(chrom)
    replicates = {side: preflight["replicates"][side] for side in ("a", "b")}
    members = {}
    for run, replicate in preflight["collapse"].items():
        members.setdefault(replicate, []).append(run)
    with tempfile.TemporaryDirectory(prefix="dapars2-", dir=str(out_dir)) as directory:
        work = Path(directory)
        names, depths = [], []
        for replicate in replicates["a"] + replicates["b"]:
            runs = sorted(members[replicate])
            sum_bedgraphs([Path(caches[run]) / "utr_coverage.bedgraph.gz" for run in runs],
                          work / (replicate + ".bedgraph"), chroms)
            depth = sum(json.loads((Path(caches[run]) / "depth.json").read_text())["mapped_reads"] for run in runs)
            names.append(replicate + ".bedgraph")
            depths.append((replicate + ".bedgraph", depth))
        (work / "depth.txt").write_text("".join("{}\t{}\n".format(n, d) for n, d in depths))
        (work / "chroms.txt").write_text("".join(c + "\n" for c in chroms))
        shutil.copyfile(utr_bed, work / "utr.bed")
        (work / "config.txt").write_text(
            "Annotated_3UTR=utr.bed\nAligned_Wig_files={}\nOutput_directory=dapars2_out/\n"
            "Output_result_file=dapars2\nsequencing_depth_file=depth.txt\nNum_Threads={}\n"
            "Coverage_threshold={}\n".format(",".join(names), threads, coverage_threshold))
        with open(out_dir / "dapars2.log", "w") as log:
            subprocess.run([sys.executable, "-c", FORK_RUNNER,
                            str(Path(dapars2_dir) / "DaPars2_Multi_Sample_Multi_Chr.py"),
                            "config.txt", "chroms.txt"], cwd=work, stdout=log, stderr=subprocess.STDOUT,
                           check=True)
        header, rows = None, []
        for chrom in chroms:
            path = work / "dapars2_out_{}".format(chrom) / "dapars2_result_temp.{}.txt".format(chrom)
            if not path.exists():
                continue
            with open(path) as stream:
                lines = stream.read().splitlines()
            if not lines:
                continue
            header = header or lines[0].split("\t")
            rows.extend(line.split("\t") for line in lines[1:] if line)
    header = header or ["Gene", "fit_value", "Predicted_Proximal_APA", "Loci"] + [
        n.rsplit(".", 1)[0] + "_PDUI" for n in names]
    with gzip.open(out_dir / "native.tsv.gz", "wt") as stream:
        stream.write("\t".join(header) + "\n")
        for row in rows:
            stream.write("\t".join(row) + "\n")
    column = {name: index for index, name in enumerate(header)}
    flagged = overlapping_opposite_strands(utr_bed)
    normalized, pvalues = [], []
    for row in rows:
        values = {}
        for side in ("a", "b"):
            values[side] = []
            for replicate in replicates[side]:
                raw = row[column[replicate + "_PDUI"]]
                if raw not in ("", "NA", "nan", "None"):
                    values[side].append(float(raw))
        effect = (sum(values["a"]) / len(values["a"]) - sum(values["b"]) / len(values["b"])
                  if values["a"] and values["b"] else None)
        p = welch(values["a"], values["b"])
        if effect is None:
            status = "not_quantified"
        elif p is None:
            status = "untestable"
        else:
            status = "tested"
        if row[column["Loci"]] in flagged:
            status += ";opposite_strand_overlap"
        transcript, gene, chrom, strand = (row[column["Gene"]].split("|") + ["", "", "", ""])[:4]
        loci = row[column["Loci"]]
        start, end = loci.split(":")[1].split("-") if ":" in loci else ("", "")
        pvalues.append(p)
        normalized.append({
            "tool": "dapars2", "feature_id": row[column["Gene"]], "gene_id": transcript,
            "gene_name": gene, "chrom": chrom, "start": start, "end": end,
            "coordinates": "0-based half-open 3' UTR; proximal site " + row[column["Predicted_Proximal_APA"]],
            "strand": strand, "event_class": "tandem_3utr_apa",
            "effect": "" if effect is None else "{:.6g}".format(effect),
            "effect_definition": PDUI_EFFECT, "n_a": len(values["a"]), "n_b": len(values["b"]),
            "coverage": "", "native_statistic": row[column["fit_value"]],
            "statistic_type": "DaPars2 fit_value (mean squared error of the two-isoform fit)",
            "p_value": "" if p is None else "{:.6g}".format(p), "status": status})
    for record, q in zip(normalized, benjamini_hochberg(pvalues)):
        record["q_value"] = "" if q is None else "{:.6g}".format(q)
        record.update({key: base.get(key, "") for key in ("reference_id", "module_reference_id",
                                                          "analysis_id", "project_id", "comparison_id")})
        record["tool_version"] = ids.get("tool_version", "")
    with gzip.open(out_dir / "results.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(normalized)
    write_status(status_path, "dapars2", "ok", [], regions=len(rows),
                 tested=sum(p is not None for p in pvalues),
                 replicates={side: replicates[side] for side in ("a", "b")}, **base)


def _empty_outputs(out_dir):
    for name in ("native.tsv.gz", "results.tsv.gz"):
        with gzip.open(Path(out_dir) / name, "wt") as stream:
            stream.write("# inference not supported\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    cov = commands.add_parser("coverage")
    cov.add_argument("--bam", required=True)
    cov.add_argument("--windows", required=True)
    cov.add_argument("--out-dir", required=True)
    cov.add_argument("--run-id", required=True)
    cov.add_argument("--identity", default="{}", help="JSON recorded in cache.json")
    cmp_ = commands.add_parser("compare")
    cmp_.add_argument("--preflight", required=True)
    cmp_.add_argument("--cache", action="append", default=[], help="RUN=DIR")
    cmp_.add_argument("--utr-bed", required=True)
    cmp_.add_argument("--dapars2-dir", required=True)
    cmp_.add_argument("--out-dir", required=True)
    cmp_.add_argument("--threads", type=int, default=1)
    cmp_.add_argument("--coverage-threshold", type=int, default=10)
    cmp_.add_argument("--ids", default="{}", help="JSON of identity fields for the results")
    args = parser.parse_args(argv)
    if args.command == "coverage":
        coverage(args.bam, args.windows, args.out_dir, args.run_id, json.loads(args.identity))
    else:
        caches = dict(item.split("=", 1) for item in args.cache)
        compare(args.preflight, caches, args.utr_bed, args.dapars2_dir, args.out_dir,
                args.threads, args.coverage_threshold, json.loads(args.ids))


if __name__ == "__main__":
    main()

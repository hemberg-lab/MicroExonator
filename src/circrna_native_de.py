"""Native differential tests of the circRNA tools, run on their own per-comparison outputs.

software    Once: CIRIquant 1.1.3's libs/CIRI_DE.R and CircTest (dieterich-lab, pinned commit),
            fetched and checked against their SHA-256; CircTest installed into <root>/R.
ciriquant   CIRI_DE_replicate as CIRIquant runs it: Rscript libs/CIRI_DE.R --lib --bsj --gene
            (edgeR glmLRT on the BSJ matrix, library sizes and normalization factors from the
            gene matrix). Inputs are written in prep_CIRIquant's formats: the BSJ matrix from
            the comparison's kept CIRIquant counts (every circRNA with a BSJ read; not-indexed
            ones left out), the gene matrix from the core run's featureCounts shards (CIRIquant's
            manual uses StringTie's prepDE.py; any gene read count matrix serves its
            normalization). Group T is side A and C is side B, so logFC is log2(A / B).
dcc-linear  DCC's host-gene count (`DCC -G`, Genecount.comb_gen_count with
            countlinearsplicedreads=False, then CombineCounts as DCC's main): reads mapped at
            each circRNA's start and end base, averaged, from the run's genome BAM (here the
            umbrella HISAT2 BAM: the STAR caches keep no BAM). Units of one replicate are
            summed, as their chimeric junctions are pooled for DCC.
dcc-circtest CircTest on DCC's CircRNACount and that LinearCount (src/circrna_circtest.R).
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from collections import defaultdict
from pathlib import Path

CIRIQUANT_TAG = "v1.1.3"
CIRI_DE_URL = "https://raw.githubusercontent.com/bioinfo-biols/CIRIquant/{}/libs/CIRI_DE.R".format(CIRIQUANT_TAG)
CIRI_DE_SHA256 = "e3c63cff7f1511c6b128f9e561489a3769e4b6ce9cdad38b92bf7c56dbacd44b"
CIRCTEST_COMMIT = "0f5a86c030dc35f5659fb3a7e9be2f16349790ff"
CIRCTEST_URL = "https://github.com/dieterich-lab/CircTest/archive/{}.tar.gz".format(CIRCTEST_COMMIT)
CIRCTEST_SHA256 = {"R/Circ.test.R": "9bc4e924d271d40131b7afb9370ce2d34610fefbf07144336b7e66baf0bdcab4",
                   "R/Circ.filter.R": "f63463fb378a0c54d82fd7ca45b52d983f11ee78c25b8239f44639b5b018c641",
                   "DESCRIPTION": "bc6c2ed34a066f4f9e851aef21ceac3306d7644fa9862f9558a38568407f352e"}
SOFTWARE = {"ciriquant": "CIRI_DE_replicate (CIRIquant 1.1.3 libs/CIRI_DE.R)",
            "dcc": "CircTest dieterich-lab/CircTest@" + CIRCTEST_COMMIT[:8] + " on DCC 0.5.0 counts"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(url, target):
    with urllib.request.urlopen(url) as response, open(target, "wb") as stream:
        shutil.copyfileobj(response, stream)


def software(root, record):
    """CIRI_DE.R and CircTest at their pins; CircTest is installed into root/R."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    fetch(CIRI_DE_URL, root / "CIRI_DE.R")
    if sha256(root / "CIRI_DE.R") != CIRI_DE_SHA256:
        raise SystemExit("CIRI_DE.R is not CIRIquant {}'s".format(CIRIQUANT_TAG))
    with tempfile.TemporaryDirectory() as directory:
        archive = Path(directory) / "circtest.tar.gz"
        fetch(CIRCTEST_URL, archive)
        subprocess.run(["tar", "-xzf", str(archive), "-C", directory], check=True)
        source = Path(directory) / ("CircTest-" + CIRCTEST_COMMIT)
        bad = [name for name, digest in CIRCTEST_SHA256.items() if sha256(source / name) != digest]
        if bad:
            raise SystemExit("CircTest files differ from the pinned commit: " + ", ".join(bad))
        (root / "R").mkdir(exist_ok=True)
        subprocess.run(["R", "CMD", "INSTALL", "--no-test-load", "-l", str(root / "R"), str(source)], check=True)
    Path(record).write_text(json.dumps({"ciri_de": {"url": CIRI_DE_URL, "sha256": CIRI_DE_SHA256},
                                        "circtest": {"url": CIRCTEST_URL, "files": CIRCTEST_SHA256},
                                        "licences": {"CIRIquant": "MIT", "CircTest": "GPL (DESCRIPTION)"}},
                                       sort_keys=True, indent=2) + "\n")


def members(preflight):
    groups = defaultdict(list)
    for unit, replicate in preflight["collapse"].items():
        groups[replicate].append(unit)
    return {replicate: sorted(units) for replicate, units in groups.items()}


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_results(out, rows):
    from src.umbrella_modules import NORMALIZED
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("\t".join(NORMALIZED) + "\n")
        for row in rows:
            stream.write("\t".join(str(row.get(column, "")) for column in NORMALIZED) + "\n")


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def base_fields(preflight, ids, tool):
    return dict(ids or {}, project_id=preflight.get("project_id", ""),
                comparison_id=preflight.get("comparison_id", ""), tool=tool)


def status_fields(base):
    return {key: value for key, value in base.items() if key != "tool"}


# ------------------------------------------------------------ CIRIquant

def gene_matrix(shards, units_of):
    """{replicate: {gene: count}} summed over each replicate's units, from featureCounts shards."""
    wanted = {unit: replicate for replicate, units in units_of.items() for unit in units}
    counts, seen = defaultdict(lambda: defaultdict(int)), set()
    for path in shards:
        with gzip.open(path, "rt") as stream:
            header = stream.readline().rstrip("\n").split("\t")
            columns = [(i, wanted[name]) for i, name in enumerate(header) if name in wanted]
            seen.update(header[i] for i, _ in columns)
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                for i, replicate in columns:
                    counts[replicate][fields[0]] += int(fields[i])
    missing = sorted(set(wanted) - seen)
    if missing:
        raise SystemExit("featureCounts shards lack runs: " + ", ".join(missing))
    return counts


def ciriquant(preflight, source, shards, ciri_de, out_dir, ids=None, rscript="Rscript"):
    from src.umbrella_modules import write_status
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    base = base_fields(preflight, ids, "ciriquant")
    sides = {side: list(preflight["replicates"][side]) for side in ("a", "b")}
    source = Path(source)
    if min(len(sides["a"]), len(sides["b"])) < 2:
        write_results(out, [])
        write_status(out / "status.json", "ciriquant", "unsupported",
                     ["CIRI_DE_replicate needs at least 2 replicates per side"], **status_fields(base))
        return
    replicates = sorted(sides["a"] + sides["b"])     # prep_CIRIquant writes samples sorted
    group = {r: "T" for r in sides["a"]}
    group.update({r: "C" for r in sides["b"]})
    bsj_rows = read_table(source / "bsj_counts.tsv.gz")
    library = {row["replicate"]: row for row in read_table(source / "library.tsv")}
    info = {row["feature_id"]: row for row in read_table(source / "results.tsv.gz")}
    kept = []
    for row in bsj_rows:
        values = [row.get(r, "NA") for r in replicates]
        if "NA" in values:
            continue                                  # not indexed in some replicate: not measured
        counts = [int(v) for v in values]
        if sum(counts):
            kept.append((row["circ_id"], counts))
    with open(out / "library.csv", "w") as stream:
        stream.write("Sample,Total,Mapped,Circular,Group\n")
        for i, r in enumerate(replicates):
            circular = 2 * sum(counts[i] for _, counts in kept)
            stream.write(",".join([r, library[r]["total_reads"], library[r]["mapped_reads"], str(circular),
                                   group[r]]) + "\n")
    with open(out / "circRNA_bsj.csv", "w") as stream:
        stream.write(",".join(["circ_id"] + replicates) + "\n")
        for cid, counts in kept:
            stream.write(",".join([cid] + [str(c) for c in counts]) + "\n")
    genes = gene_matrix(shards, {r: members(preflight)[r] for r in replicates})
    names = sorted(set().union(*(genes[r] for r in replicates)))
    with open(out / "gene_count_matrix.csv", "w") as stream:
        stream.write(",".join(["gene_id"] + replicates) + "\n")
        for gene in names:
            stream.write(",".join([gene] + [str(genes[r].get(gene, 0)) for r in replicates]) + "\n")
    if not kept:
        write_results(out, [])
        write_status(out / "status.json", "ciriquant", "empty", ["no circRNA has a BSJ read"], **status_fields(base))
        return
    with open(out / "CIRI_DE.log", "w") as log:
        subprocess.run([rscript, str(ciri_de), "--lib=" + str(out / "library.csv"),
                        "--bsj=" + str(out / "circRNA_bsj.csv"), "--gene=" + str(out / "gene_count_matrix.csv"),
                        "--out=" + str(out / "circ_de.csv"), "--out2=" + str(out / "gene_de.csv")],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
    if not (out / "circ_de.csv").is_file():
        raise SystemExit("CIRI_DE.R wrote no circ_de.csv (see CIRI_DE.log)")
    rows = []
    with open(out / "circ_de.csv") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        for fields in reader:
            record = dict(zip(["circ_id"] + header[1:], fields))
            feature = info.get(record["circ_id"], {})
            rows.append(dict(base, tool_version=SOFTWARE["ciriquant"], feature_id=record["circ_id"],
                             chrom=feature.get("chrom", ""), start=feature.get("start", ""), end=feature.get("end", ""),
                             coordinates=record["circ_id"], strand=feature.get("strand", ""),
                             event_class="back_splice_junction", effect=record["logFC"],
                             effect_definition="edgeR log2 fold change of BSJ counts, A over B (CIRI_DE.R)",
                             n_a=len(sides["a"]), n_b=len(sides["b"]), coverage=feature.get("coverage", ""),
                             native_statistic=record["LR"], statistic_type="edgeR glmLRT likelihood ratio (CIRI_DE.R)",
                             p_value=record["PValue"], q_value=record["FDR"], status="tested"))
    rows.sort(key=lambda row: row["feature_id"])
    write_results(out, rows)
    write_status(out / "status.json", "ciriquant", "ok",
                 ["CIRI_DE_replicate (CIRIquant 1.1.3 libs/CIRI_DE.R) unmodified",
                  "gene matrix: core HISAT2 + featureCounts shards (CIRIquant's manual uses StringTie prepDE.py)"],
                 circrnas=len(kept), tested=len(rows), genes=len(names), **status_fields(base))


# ------------------------------------------------------------ DCC

def dcc_linear(preflight, dcc_dir, bams, genome, out_path, tmp_root=None):
    """DCC's own host-gene count (DCC -G) at the circRNAs DCC kept, one column per replicate."""
    from DCC import CombineCounts, genecount
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    replicates = [r for side in ("a", "b") for r in preflight["replicates"][side]]
    units_of = members(preflight)
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(Path(out_path).parent)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dcc-linear-", dir=tmp_root) as directory:
        work = Path(directory)
        with gzip.open(Path(dcc_dir) / "CircRNACount.tsv.gz", "rt") as stream:
            lines = [line for line in stream if line.strip()]
        header = lines[0].rstrip("\n").split("\t")
        samples = header[4:]
        if len(samples) != len(replicates):
            raise SystemExit("CircRNACount has {} samples, the comparison {} replicates".format(len(samples), len(replicates)))
        (work / "CircRNACount").write_text("".join(lines))
        per_replicate = []
        for index, replicate in enumerate(replicates):
            totals = None
            for unit in units_of[replicate]:
                counter = genecount.Genecount(str(work) + "/")
                output = str(work / "{}_{}.linear".format(index, unit))
                counter.comb_gen_count(str(work / "CircRNACount"), str(Path(bams[unit]).resolve()),
                                       str(Path(genome).resolve()), output, countlinearsplicedreads=False)
                values = [line.rstrip("\n").split("\t") for line in open(output)]
                totals = values if totals is None else [row[:5] + [str(int(row[5]) + int(old[5]))]
                                                        for row, old in zip(values, totals)]
            path = work / "s{:03d}.linear".format(index)
            path.write_text("".join("\t".join(row) + "\n" for row in totals))
            per_replicate.append(str(path))
        combiner = CombineCounts.Combine(str(work) + "/")
        result = combiner.combine(per_replicate, col=6, circ=False)
        combiner.writeouput_linear(str(work / "LinearCount"), result, "\t".join(samples), header=True)
        with open(work / "LinearCount") as source, gzip.open(out_path, "wt") as sink:
            shutil.copyfileobj(source, sink)


def dcc_circtest(preflight, dcc_dir, linear, circtest_lib, out_dir, settings=None, ids=None, rscript="Rscript"):
    from src.umbrella_modules import write_status
    try:
        from src.run_circrna_discovery import circ_id
    except ImportError:
        from run_circrna_discovery import circ_id
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    settings = dict({"filter_count": 5, "filter_sample": 3, "percentage": 0.1}, **(settings or {}))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    base = base_fields(preflight, ids, "dcc")
    n_a, n_b = (len(preflight["replicates"][side]) for side in ("a", "b"))
    if min(n_a, n_b) < 2:
        write_results(out, [])
        write_status(out / "status.json", "dcc", "unsupported", ["CircTest needs at least 2 replicates per side"],
                     **status_fields(base))
        return
    script = Path(__file__).resolve().parent / "circrna_circtest.R"
    with open(out / "circtest.log", "w") as log:
        subprocess.run([rscript, str(script), str(Path(dcc_dir) / "CircRNACount.tsv.gz"), str(linear), str(n_a),
                        str(n_b), str(settings["filter_count"]), str(settings["filter_sample"]),
                        str(settings["percentage"]), str(circtest_lib), str(out / "circtest.tsv")],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
    table = read_table(out / "circtest.tsv")
    rows = []
    for row in table:
        if row["status"] != "tested":
            continue
        start, end = int(row["Start"]), int(row["End"])
        a, b = number(row["ratio_a"]), number(row["ratio_b"])
        rows.append(dict(base, tool_version=SOFTWARE["dcc"], feature_id=circ_id(row["Chr"], start, end), gene_id="",
                         gene_name="", chrom=row["Chr"], start=start, end=end, coordinates=circ_id(row["Chr"], start, end),
                         strand=row["Strand"], event_class="back_splice_junction",
                         effect="" if a is None or b is None else "{:.6f}".format(a - b),
                         effect_definition="mean circ / (circ + linear) ratio per replicate, A minus B (CircTest)",
                         n_a=n_a, n_b=n_b, coverage=row["circ_reads"], native_statistic=row["p_value"],
                         statistic_type="CircTest beta-binomial likelihood-ratio p-value (aod::betabin), BH q",
                         p_value=row["p_value"], q_value=row["q_value"], status="tested"))
    write_results(out, rows)
    filtered = sum(row["status"] == "filtered" for row in table)
    failed = sum(row["status"] == "fit_failed" for row in table)
    reasons = ["CircTest Circ.test per circRNA on DCC CircRNACount and DCC host-gene LinearCount",
               "LinearCount from the umbrella HISAT2 BAM (DCC's own counting; the STAR caches keep no BAM)"]
    if n_a != n_b:
        reasons.append("unequal groups: Circ.filter's rule applied per group in circrna_circtest.R")
    write_status(out / "status.json", "dcc", "ok" if rows else "no_tests", reasons, circrnas=len(table),
                 tested=len(rows), filtered=filtered, fit_failed=failed, settings=settings, **status_fields(base))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    s = commands.add_parser("software")
    s.add_argument("--root", required=True)
    s.add_argument("--record", required=True)
    c = commands.add_parser("ciriquant")
    c.add_argument("--preflight", required=True)
    c.add_argument("--source", required=True, help="the comparison's CIRIquant analysis folder")
    c.add_argument("--gene-shard", action="append", default=[])
    c.add_argument("--ciri-de", required=True)
    c.add_argument("--out-dir", required=True)
    c.add_argument("--ids", default="{}")
    lin = commands.add_parser("dcc-linear")
    lin.add_argument("--preflight", required=True)
    lin.add_argument("--dcc-dir", required=True)
    lin.add_argument("--bam", action="append", default=[], help="unit=bam")
    lin.add_argument("--genome", required=True)
    lin.add_argument("--out", required=True)
    t = commands.add_parser("dcc-circtest")
    t.add_argument("--preflight", required=True)
    t.add_argument("--dcc-dir", required=True)
    t.add_argument("--linear", required=True)
    t.add_argument("--circtest-lib", required=True)
    t.add_argument("--out-dir", required=True)
    t.add_argument("--settings", default="{}")
    t.add_argument("--ids", default="{}")
    args = parser.parse_args(argv)
    if args.command == "software":
        software(args.root, args.record)
    elif args.command == "ciriquant":
        ciriquant(args.preflight, args.source, args.gene_shard, args.ciri_de, args.out_dir, json.loads(args.ids))
    elif args.command == "dcc-linear":
        dcc_linear(args.preflight, args.dcc_dir, dict(item.split("=", 1) for item in args.bam), args.genome, args.out)
    else:
        dcc_circtest(args.preflight, args.dcc_dir, args.linear, args.circtest_lib, args.out_dir,
                     json.loads(args.settings), json.loads(args.ids))


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    sys.exit(main())

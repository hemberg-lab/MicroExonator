"""MAJIQ v3 module: per-run SJ caches and per-comparison HET quantification.

MAJIQ v3 (package rna_majiq, https://bitbucket.org/biociphers/majiq_academic)
is licensed software: it runs from an installation Guillermo makes under its
licence (`majiq_bin_folder`), or one built from the source he downloads
(`majiq_source`, rule umbrella_majiq_software). Nothing here downloads or
redistributes it. Command syntax checked against the rna_majiq source
(majiq_academic main, read 2026-10-01).

sj        Per run, from the temporary HISAT2 BAM:
            majiq-build sj BAM SPLICEGRAPH OUT.sj --prefix RUN --strandness S --nthreads N
          The SJ file (a zarr folder) is MAJIQ's reusable per-run input; it is
          kept. --prefix names the experiment by run: every BAM is called
          aligned.bam, and MAJIQ would otherwise name them all "aligned".
compare   Per comparison, from the kept SJ files only:
            majiq-build update  splicegraph from the selected runs only, one
                                build group per side (--groups-tsv group/sj),
                                frozen under the analysis ID
            majiq psi-coverage  one PsiCoverage per biological replicate; the
                                technical runs of a replicate are summed with
                                rna_majiq's PsiCoverage.sum (the CLI keeps one
                                experiment per SJ file)
            majiq heterogen     replicate-aware HET, groups a (-psi1) and b (-psi2)
          Each step's extra arguments can be replaced from the config
          (majiq_update_args, majiq_psicov_args, majiq_heterogen_args); the
          `--help` of every command used is written to the log.
          Outputs: native.tsv.gz (HET table as written, with splicegraph
          annotation), results.tsv.gz (one row per LSV connection: effect =
          median PSI(A) - median PSI(B) of the per-replicate posterior means;
          p = the chosen HET test's raw p-value; q = Benjamini-Hochberg over
          the tested connections), status.json.
"""

import argparse
import csv
import gzip
import json
import os
import shlex
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

if __package__:
    from src.run_dapars2 import benjamini_hochberg
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from run_dapars2 import benjamini_hochberg
    from umbrella_modules import NORMALIZED, write_status

STRANDNESS = {"firststrand": "REVERSE", "secondstrand": "FORWARD",
              "unstranded": "NONE", "auto": "AUTO"}
DEFAULTS = {
    "update": "--min-experiments 0.5",
    "psicov": "",
    "heterogen": "--stats ttest mannwhitneyu tnom",
}
STATISTICS = ("ttest", "mannwhitneyu", "tnom", "infoscore")
HET_EFFECT = ("median over replicates of PSI(A) - median of PSI(B) for this connection of its LSV "
              "(MAJIQ HET raw_psi_quantile_0.500; PSI from 0 to 1)")
# sums a replicate's technical runs into one PsiCoverage, inside MAJIQ's own Python;
# min_experiments_f=1: a connection passes when at least one run passes
POOL = ("import sys, rna_majiq as nm; "
        "nm.PsiCoverage.from_zarr(sys.argv[3:]).sum(sys.argv[2], min_experiments_f=1).to_zarr(sys.argv[1])")


def _tool(bin_dir, name):
    return str(Path(bin_dir) / name) if bin_dir else name


def _run(command, log, cwd=None):
    log.write("$ {}\n".format(" ".join(shlex.quote(str(part)) for part in command)))
    log.flush()
    subprocess.run([str(part) for part in command], stdout=log, stderr=subprocess.STDOUT,
                   check=True, cwd=cwd)


def _help(command, log):
    log.write("--- {} --help\n".format(" ".join(command)))
    log.flush()
    subprocess.run([str(part) for part in command] + ["--help"], stdout=log,
                   stderr=subprocess.STDOUT, check=False)


def sj(bam, splicegraph, out_dir, run_id, strandedness, threads=1, bin_dir="", extra=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    strand = STRANDNESS.get(strandedness, strandedness)
    temporary = out_dir / "run.sj.tmp"
    with open(out_dir / "sj.log", "w") as log:
        _run([_tool(bin_dir, "majiq-build"), "sj", bam, splicegraph, temporary,
              "--prefix", run_id, "--strandness", strand, "--nthreads", threads], log)
    final = out_dir / "run.sj"
    if final.is_dir():          # left by an interrupted earlier attempt (no marker)
        shutil.rmtree(final)
    if temporary.is_dir():
        shutil.move(str(temporary), str(final))
    else:
        os.replace(temporary, final)
    (out_dir / "cache.json").write_text(json.dumps(dict(
        extra or {}, run_id=run_id, tool="majiq", strandness=strand), sort_keys=True, indent=2) + "\n")


def compare(preflight_path, caches, splicegraph, out_dir, threads=1, bin_dir="", args=None, ids=None,
            statistic="ttest"):
    if statistic not in STATISTICS:
        raise ValueError("majiq_het_statistic must be one of {}".format(", ".join(STATISTICS)))
    preflight = json.loads(Path(preflight_path).read_text())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    args = dict(DEFAULTS, **(args or {}))
    ids = ids or {}
    base = dict(project_id=preflight["project_id"], comparison_id=preflight["comparison_id"], **ids)
    if not preflight["inference_supported"]:
        for name in ("native.tsv.gz", "results.tsv.gz"):
            with gzip.open(out_dir / name, "wt") as stream:
                stream.write("# inference not supported\n")
        write_status(out_dir / "status.json", "majiq", "unsupported", preflight["reasons"], **base)
        return
    members = defaultdict(list)
    for run, replicate in preflight["collapse"].items():
        members[replicate].append(run)
    majiq, build = _tool(bin_dir, "majiq"), _tool(bin_dir, "majiq-build")
    # absolute: every path is used from inside the private folder's commands
    sj_of = {run: (Path(caches[run]) / "run.sj").resolve() for run in caches}
    with tempfile.TemporaryDirectory(prefix="majiq-", dir=str(out_dir)) as directory, \
            open(out_dir / "majiq.log", "w") as log:
        work = Path(directory).resolve()
        for command in ([build, "update"], [majiq, "psi-coverage"], [majiq, "heterogen"]):
            _help(command, log)
        with open(work / "groups.tsv", "w") as stream:
            stream.write("group\tsj\n")
            for side in ("a", "b"):
                for replicate in preflight["replicates"][side]:
                    for run in sorted(members[replicate]):
                        stream.write("{}\t{}\n".format(side, sj_of[run]))
        graph = work / "splicegraph.zarr"
        _run([build, "update", Path(splicegraph).resolve(), graph, "--groups-tsv", work / "groups.tsv",
              "--nthreads", threads] + shlex.split(args["update"]), log)
        coverage = {}
        for side in ("a", "b"):
            coverage[side] = []
            for replicate in preflight["replicates"][side]:
                runs = sorted(members[replicate])
                path = work / "{}.psicov".format(replicate)
                if len(runs) == 1:
                    _run([majiq, "psi-coverage", graph, path, sj_of[runs[0]], "--prefixes", replicate,
                          "--nthreads", threads] + shlex.split(args["psicov"]), log)
                else:
                    runs_path = work / "{}.runs.psicov".format(replicate)
                    _run([majiq, "psi-coverage", graph, runs_path] + [sj_of[run] for run in runs]
                         + ["--prefixes"] + runs + ["--nthreads", threads] + shlex.split(args["psicov"]), log)
                    _run([_tool(bin_dir, "python"), "-c", POOL, path, replicate, runs_path], log)
                coverage[side].append(path)
        native = work / "heterogen.tsv"
        _run([majiq, "heterogen", "-psi1"] + coverage["a"] + ["-psi2"] + coverage["b"]
             + ["-n", "a", "b", "--splicegraph", graph, "--output-tsv", native, "--nthreads", threads]
             + shlex.split(args["heterogen"]), log)
        with open(native, "rb") as source, gzip.open(out_dir / "native.tsv.gz", "wb") as target:
            shutil.copyfileobj(source, target)
        # the frozen splicegraph is small next to the SJ inputs; keep it for export and VOILA
        shutil.make_archive(str(out_dir / "splicegraph"), "gztar", root_dir=str(work), base_dir=graph.name)
    normalize(out_dir / "native.tsv.gz", out_dir, preflight, base, statistic)


def _number(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return None if value != value else value      # NaN is missing


def normalize(native, out_dir, preflight, base, statistic="ttest"):
    with gzip.open(native, "rt") as stream:
        lines = [line for line in stream if not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    header = reader.fieldnames or []
    median = {side: "{}-raw_psi_quantile_0.500".format(side) for side in ("a", "b")}
    # one test: MAJIQ writes raw_pvalue without the test's name in front
    pvalue = "{}-raw_pvalue".format(statistic)
    if pvalue not in header and "raw_pvalue" in header:
        pvalue = "raw_pvalue"
    needed = ["gene_id", "seqid", "start", "end", median["a"], median["b"], pvalue]
    missing = [name for name in needed if name not in header]
    records, pvalues = [], []
    if not missing:
        for row in reader:
            a, b = _number(row[median["a"]]), _number(row[median["b"]])
            p = _number(row[pvalue])
            intron = row.get("is_intron", "").lower() in ("true", "1")
            kind = "intron" if intron else "junction"
            records.append(dict(
                base, tool="majiq", tool_version=base.get("tool_version", ""),
                feature_id="{}:{}:{}-{}:{}:{}-{}".format(
                    row["gene_id"], row.get("event_type", ""), row.get("ref_exon_start", ""),
                    row.get("ref_exon_end", ""), kind, row["start"], row["end"]),
                gene_id=row["gene_id"], gene_name=row.get("gene_name", ""), chrom=row["seqid"],
                start=row["start"], end=row["end"], strand=row.get("strand", ""),
                coordinates="MAJIQ {} {}-{}; reference exon {}-{}; other exon {}-{}".format(
                    kind, row["start"], row["end"], row.get("ref_exon_start", ""),
                    row.get("ref_exon_end", ""), row.get("other_exon_start", ""),
                    row.get("other_exon_end", "")),
                event_class="lsv_" + kind,
                effect="" if a is None or b is None else "{:.6g}".format(a - b),
                effect_definition=HET_EFFECT,
                n_a=len(preflight["replicates"]["a"]), n_b=len(preflight["replicates"]["b"]),
                coverage="replicates passing: a {}, b {}".format(
                    row.get("a-num_passed", ""), row.get("b-num_passed", "")),
                native_statistic=row.get("tnom_score", ""),
                statistic_type="MAJIQ HET {} p-value on per-replicate PSI (native_statistic: TNOM score)".format(
                    statistic),
                p_value="" if p is None else "{:.6g}".format(p),
                status="tested" if p is not None else "untested"))
            pvalues.append(p)
    for record, q in zip(records, benjamini_hochberg(pvalues)):
        record["q_value"] = "" if q is None else "{:.6g}".format(q)
    with gzip.open(Path(out_dir) / "results.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    tested = sum(p is not None for p in pvalues)
    if missing:
        status, reasons = "native_only", ["HET columns not found: {} (have: {})".format(
            ", ".join(missing), ", ".join(header))]
    elif not records:
        status, reasons = "empty", ["MAJIQ HET reported no quantifiable LSV connections"]
    elif not tested:
        status, reasons = "no_tests", ["no connection has a {} p-value".format(statistic)]
    else:
        status, reasons = "ok", []
    write_status(Path(out_dir) / "status.json", "majiq", status, reasons, connections=len(records),
                 tested=tested, statistic=statistic, **base)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    s = commands.add_parser("sj")
    s.add_argument("--bam", required=True)
    s.add_argument("--splicegraph", required=True)
    s.add_argument("--out-dir", required=True)
    s.add_argument("--run-id", required=True)
    s.add_argument("--strandedness", required=True)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--bin-dir", default="")
    s.add_argument("--identity", default="{}")
    c = commands.add_parser("compare")
    c.add_argument("--preflight", required=True)
    c.add_argument("--cache", action="append", default=[], help="RUN=DIR")
    c.add_argument("--splicegraph", required=True)
    c.add_argument("--out-dir", required=True)
    c.add_argument("--threads", type=int, default=1)
    c.add_argument("--bin-dir", default="")
    c.add_argument("--args", default="{}", help="JSON overrides for update/psicov/heterogen")
    c.add_argument("--statistic", default="ttest", choices=STATISTICS,
                   help="HET test whose p-value fills p_value/q_value")
    c.add_argument("--ids", default="{}")
    args = parser.parse_args(argv)
    if args.command == "sj":
        sj(args.bam, args.splicegraph, args.out_dir, args.run_id, args.strandedness,
           args.threads, args.bin_dir, json.loads(args.identity))
    else:
        compare(args.preflight, dict(item.split("=", 1) for item in args.cache), args.splicegraph,
                args.out_dir, args.threads, args.bin_dir, json.loads(args.args), json.loads(args.ids),
                args.statistic)


if __name__ == "__main__":
    main()

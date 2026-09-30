"""MAJIQ v3 module: per-run SJ caches and per-comparison HET quantification.

MAJIQ v3 is licensed software: it runs from an installation built from the
source Guillermo downloads under its licence (rule umbrella_majiq_software)
or from `majiq_bin_folder`. Nothing here downloads or redistributes it.

sj        Per run, from the temporary HISAT2 BAM:
            majiq-build sj BAM SPLICEGRAPH OUT.sj --strandness S --nthreads N
          The .sj file is MAJIQ's reusable per-run input; it is kept.
compare   Per comparison, from the kept .sj files only:
            majiq-build update  (splicegraph from the selected runs only, frozen
                                 under the analysis ID)
            majiq psi-coverage  (one PsiCoverage per biological replicate; its
                                 technical runs' SJ are pooled at coverage level)
            majiq heterogen     (replicate-aware HET, A = psi1, B = psi2)
          Each step's arguments can be replaced from the config
          (majiq_update_args, majiq_psicov_args, majiq_heterogen_args), and the
          `--help` of every command used is written to the log, because the
          v3 option names were not verifiable without the licensed binary.
          Outputs: native.tsv.gz (HET table as written), results.tsv.gz
          (normalized; statistics kept as native TNOM/t-test/Wilcoxon values,
          never relabelled as adjusted p-values), status.json.
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
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from umbrella_modules import NORMALIZED, write_status

STRANDNESS = {"firststrand": "REVERSE", "secondstrand": "FORWARD",
              "unstranded": "NONE", "auto": "AUTO"}
DEFAULTS = {
    "update": "--min-experiments 0.5",
    "psicov": "",
    "heterogen": "--stats ttest wilcoxon tnom",
}
HET_EFFECT = "median PSI(A) - median PSI(B) of the junction within its LSV (MAJIQ HET)"


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
              "--strandness", strand, "--nthreads", threads], log)
    if temporary.is_dir():
        shutil.move(str(temporary), str(out_dir / "run.sj"))
    else:
        os.replace(temporary, out_dir / "run.sj")
    (out_dir / "cache.json").write_text(json.dumps(dict(
        extra or {}, run_id=run_id, tool="majiq", strandness=strand), sort_keys=True, indent=2) + "\n")


def compare(preflight_path, caches, splicegraph, out_dir, threads=1, bin_dir="", args=None, ids=None):
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
    with tempfile.TemporaryDirectory(prefix="majiq-", dir=str(out_dir)) as directory, \
            open(out_dir / "majiq.log", "w") as log:
        work = Path(directory)
        for command in ([build, "update"], [majiq, "psi-coverage"], [majiq, "heterogen"]):
            _help(command, log)
        with open(work / "experiments.tsv", "w") as stream:
            stream.write("path\tgroup\n")
            for side in ("a", "b"):
                for replicate in preflight["replicates"][side]:
                    for run in sorted(members[replicate]):
                        stream.write("{}\t{}\n".format(Path(caches[run]) / "run.sj", side))
        graph = work / "splicegraph.zarr"
        _run([build, "update", splicegraph, graph, "--experiments", work / "experiments.tsv",
              "--nthreads", threads] + shlex.split(args["update"]), log)
        coverage = {}
        for side in ("a", "b"):
            coverage[side] = []
            for replicate in preflight["replicates"][side]:
                path = work / "{}.psicov".format(replicate)
                _run([majiq, "psi-coverage", graph, path]
                     + [Path(caches[run]) / "run.sj" for run in sorted(members[replicate])]
                     + ["--nthreads", threads] + shlex.split(args["psicov"]), log)
                coverage[side].append(path)
        native = work / "heterogen.tsv"
        _run([majiq, "heterogen", graph, "--psi1"] + coverage["a"] + ["--psi2"] + coverage["b"]
             + ["--output-tsv", native, "--nthreads", threads] + shlex.split(args["heterogen"]), log)
        with open(native, "rb") as source, gzip.open(out_dir / "native.tsv.gz", "wb") as target:
            shutil.copyfileobj(source, target)
        # the frozen splicegraph is small next to the SJ inputs; keep it for export
        shutil.make_archive(str(out_dir / "splicegraph"), "gztar", root_dir=str(work), base_dir=graph.name)
    normalize(out_dir / "native.tsv.gz", out_dir, preflight, base)


def _pick(header, *candidates):
    lowered = {name.lower(): name for name in header}
    for candidate in candidates:
        for name in lowered:
            if candidate in name:
                return lowered[name]
    return None


def normalize(native, out_dir, preflight, base):
    with gzip.open(native, "rt") as stream:
        lines = [line for line in stream if not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    header = reader.fieldnames or []
    lsv = _pick(header, "lsv_id")
    gene = _pick(header, "gene_id")
    name = _pick(header, "gene_name")
    effect = _pick(header, "median_dpsi", "dpsi")
    tnom = _pick(header, "tnom")
    coordinates = _pick(header, "junction_coord", "junctions_coords", "coord")
    status = "ok" if lsv and effect else "native_only"
    records = []
    if status == "ok":
        for row in reader:
            records.append(dict(base, tool="majiq", tool_version=base.get("tool_version", ""),
                                feature_id=row.get(lsv, ""), gene_id=row.get(gene, "") if gene else "",
                                gene_name=row.get(name, "") if name else "",
                                coordinates=row.get(coordinates, "") if coordinates else "",
                                event_class="lsv", effect=row.get(effect, ""),
                                effect_definition=HET_EFFECT + " [column {}]".format(effect),
                                n_a=len(preflight["replicates"]["a"]),
                                n_b=len(preflight["replicates"]["b"]),
                                native_statistic=row.get(tnom, "") if tnom else "",
                                statistic_type="TNOM score (MAJIQ HET; heuristic, not an FDR)",
                                p_value="", q_value="", status="quantified"))
    with gzip.open(Path(out_dir) / "results.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    reasons = [] if status == "ok" else ["HET columns not recognised: {}".format(", ".join(header))]
    write_status(Path(out_dir) / "status.json", "majiq", status, reasons, lsv_rows=len(records), **base)


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
    c.add_argument("--ids", default="{}")
    args = parser.parse_args(argv)
    if args.command == "sj":
        sj(args.bam, args.splicegraph, args.out_dir, args.run_id, args.strandedness,
           args.threads, args.bin_dir, json.loads(args.identity))
    else:
        compare(args.preflight, dict(item.split("=", 1) for item in args.cache), args.splicegraph,
                args.out_dir, args.threads, args.bin_dir, json.loads(args.args), json.loads(args.ids))


if __name__ == "__main__":
    main()

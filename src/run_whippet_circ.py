"""Whippet with --circ: a parallel Whippet quantification that also reports back-splicing.

quant    Per run. `whippet-quant --circ` on the staged reads and the umbrella Whippet
         index, with the core quant's whippet_flags, so its linear events can be set
         against the core Whippet results. Keeps psi, jnc, gene and isoform TPM (not the
         read map) and cache.json with the Whippet version. Never touches the core
         Whippet outputs in umbrella/work/.../whippet/.
compare  Per comparison, from the kept caches only:
         native.diff.gz   whippet-delta itself on the --circ psi files, with its defaults
                          (as the core umbrella_whippet_delta). It skips BS lines
                          (src/diff.jl), so this is the linear events only.
         bs_counts.tsv.gz every BS line per run: Psi, Total_Reads, back-splice reads.
         results.tsv.gz   BS events tested with whippet-delta's model (src/me_delta.py, a
                          port of src/diff.jl): per sample Beta(Psi*N + 1, (1 - Psi)*N + 1),
                          replicate draws pooled and refitted, Probability =
                          max P(A - B > 0), P(B - A > 0). whippet-delta cannot do this
                          itself: it walks the sample files line by line and BS lines
                          differ between samples, so here they are matched by gene, node
                          pair, coordinates and strand. A sample without a BS line for an
                          event is not used for it, as whippet-delta treats a missing Psi.

Whippet's BS Psi (src/io.jl output_circular) is back-splice reads over back-splice reads
plus the linear reads leaving the same node; Total_Reads is that sum.
"""

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import zlib
from collections import defaultdict
from pathlib import Path

KEPT = ("psi.gz", "jnc.gz", "gene.tpm.gz", "isoform.tpm.gz")
BS_COLUMNS = ("gene", "node", "coord", "strand")
DIFF_HEADER = ("Gene\tNode\tCoord\tStrand\tType\tPsi_A\tPsi_B\tDeltaPsi\tProbability\tComplexity\tEntropy\n")


def whippet_version(whippet_bin):
    """Version line of the Whippet.jl checkout that holds bin/ (Project.toml)."""
    project = Path(whippet_bin).resolve().parent / "Project.toml"
    if project.is_file():
        for line in project.read_text().splitlines():
            if line.strip().startswith("version"):
                return line.split("=", 1)[1].strip().strip('"')
    return "unknown"


def count_lines(path):
    events = circular = 0
    with gzip.open(path, "rt") as stream:
        next(stream, None)
        for line in stream:
            events += 1
            circular += line.split("\t", 5)[4] == "BS"
    return events, circular


def quant(reads, index, out_dir, run_id, julia="julia", whippet_bin="", flags="", identity=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "tmp_quant"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir()
    command = ([julia, str(Path(whippet_bin) / "whippet-quant.jl")] + list(reads) +
               ["-x", str(index), "-o", str(work / "quant")] + flags.split() + ["--circ"])
    with (out / "quant.log").open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    for suffix in KEPT:
        source = work / ("quant." + suffix)
        if not source.is_file():
            raise SystemExit("whippet-quant --circ wrote no quant." + suffix)
        shutil.move(str(source), str(out / ("quant." + suffix)))
    shutil.rmtree(work)
    events, circular = count_lines(out / "quant.psi.gz")
    record = dict(identity or {}, run_id=run_id, tool="whippet_circ", whippet_version=whippet_version(whippet_bin),
                  command=" ".join(command[1:]).replace(str(work / "quant"), "<out>/quant"),
                  events=events, bs_events=circular)
    (out / "cache.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    return record


def read_bs(path):
    """{(gene, node, coord, strand): (psi, total_reads)} for the BS lines of one psi.gz."""
    events = {}
    with gzip.open(path, "rt") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        column = {name: i for i, name in enumerate(header)}
        for line in stream:
            row = line.rstrip("\n").split("\t")
            if row[column["Type"]] != "BS":
                continue
            psi, total = row[column["Psi"]], row[column["Total_Reads"]]
            if psi in ("NA", "") or total in ("NA", ""):
                continue
            key = (row[column["Gene"]], row[column["Node"]], row[column["Coord"]], row[column["Strand"]])
            events[key] = (float(psi), float(total))
    return events


def parse_coord(coord):
    chrom, span = coord.rsplit(":", 1)
    first, second = (int(value) for value in span.split("-"))
    return chrom, min(first, second), max(first, second)


def event_rng(seed, key):
    import numpy
    return numpy.random.default_rng([seed, zlib.crc32("|".join(key).encode())])


def bs_delta(samples, runs, settings):
    """[(key, result or None, values)] for every BS event seen in any run."""
    try:
        from src.me_delta import compare as delta
    except ImportError:
        from me_delta import compare as delta
    keys = sorted(set().union(*(samples[run] for side in ("a", "b") for run in runs[side])))
    tested = []
    for key in keys:
        values = {side: [samples[run][key] for run in runs[side] if key in samples[run]] for side in ("a", "b")}
        result = delta(values["a"], values["b"], event_rng(settings["seed"], key), settings["min_reads"],
                       settings["min_samples"], settings["size"])
        tested.append((key, result, values))
    return tested


def native_delta(julia, whippet_bin, files, prefix):
    # whippet-delta reads -a/-b without a comma as a glob pattern; the trailing comma keeps a
    # one-run side a file list
    subprocess.run([julia, str(Path(whippet_bin) / "whippet-delta.jl"), "-a", ",".join(files["a"]) + ",",
                    "-b", ",".join(files["b"]) + ",", "-o", str(prefix)], check=True)


def compare(preflight, caches, out_dir, julia="julia", whippet_bin="", ids=None, settings=None, selection=None,
            native=native_delta):
    from src.umbrella_modules import NORMALIZED, write_status
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    settings = dict({"min_reads": 5, "min_samples": 1, "size": 1000, "seed": 123456}, **(settings or {}))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    base = dict(ids or {}, project_id=preflight.get("project_id", ""),
                comparison_id=preflight.get("comparison_id", ""), selection=selection or {})
    runs = {side: list(preflight["runs"][side]) for side in ("a", "b")}
    if not all(runs.values()):
        _empty(out)
        write_status(out / "status.json", "whippet_circ", "unsupported", ["a side has no run"], **base)
        return
    records = {}
    for run in runs["a"] + runs["b"]:
        record = json.loads((Path(caches[run]) / "cache.json").read_text())
        if record.get("run_id") != run:
            raise SystemExit("whippet_circ cache {} belongs to {}, not {}".format(caches[run], record.get("run_id"), run))
        records[run] = record
    if len({record.get("module_reference_id") for record in records.values()}) > 1:
        raise SystemExit("whippet_circ caches come from different settings")
    versions = sorted({record.get("whippet_version", "unknown") for record in records.values()})
    version = "Whippet " + "/".join(versions) + " --circ"
    reasons = ["each run is one whippet-delta sample, as in the core umbrella_whippet_delta"]

    files = {side: [str((Path(caches[run]) / "quant.psi.gz").resolve()) for run in runs[side]] for side in ("a", "b")}
    native(julia, whippet_bin, files, out / "native")
    if not (out / "native.diff.gz").is_file():
        raise SystemExit("whippet-delta wrote no native.diff.gz")

    samples = {run: read_bs(Path(caches[run]) / "quant.psi.gz") for run in runs["a"] + runs["b"]}
    order = runs["a"] + runs["b"]
    with gzip.open(out / "bs_counts.tsv.gz", "wt") as stream:
        stream.write("\t".join(BS_COLUMNS + tuple("{}_{}".format(run, field) for run in order
                                                  for field in ("psi", "total_reads", "bs_reads"))) + "\n")
        for key in sorted(set().union(*samples.values())):
            cells = []
            for run in order:
                if key in samples[run]:
                    psi, total = samples[run][key]
                    cells += ["{:.6g}".format(psi), "{:.6g}".format(total), "{:.6g}".format(psi * total)]
                else:
                    cells += ["NA", "NA", "NA"]
            stream.write("\t".join(key + tuple(cells)) + "\n")

    rows = []
    for key, result, values in bs_delta(samples, runs, settings):
        if result is None:
            continue
        gene, node, coord, strand = key
        chrom, start, end = parse_coord(coord)
        rows.append(dict(base, tool="whippet_circ", tool_version=version, feature_id="{}:{}:{}".format(gene, node, coord),
                         gene_id=gene, gene_name="", chrom=chrom, start=start, end=end, coordinates=coord,
                         strand=strand, event_class="back_splice_junction",
                         effect="{:.6f}".format(result["DeltaPsi"]),
                         effect_definition="Whippet BS Psi (back-splice / back-splice + linear reads from the "
                                           "node), mean of the fitted Beta, A minus B",
                         n_a=result["Samples_A"], n_b=result["Samples_B"],
                         coverage="{:.6g}".format(result["Reads_A"] + result["Reads_B"]),
                         native_statistic="{:.3f}".format(result["Probability"]),
                         statistic_type="whippet-delta Beta-posterior Probability (src/me_delta.py port) on BS events",
                         p_value="", q_value="", status="tested"))
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("\t".join(NORMALIZED) + "\n")
        for row in rows:
            stream.write("\t".join(str(row.get(column, "")) for column in NORMALIZED) + "\n")
    detected = len(set().union(*samples.values()))
    if not rows:
        write_status(out / "status.json", "whippet_circ", "empty" if not detected else "no_tests",
                     reasons + ["no BS event had min_reads in min_samples runs of both sides"],
                     bs_events=detected, settings=settings, **base)
        return
    write_status(out / "status.json", "whippet_circ", "ok", reasons, bs_events=detected, tested=len(rows),
                 settings=settings, whippet_versions=versions, **base)


def _empty(out):
    from src.umbrella_modules import NORMALIZED
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("\t".join(NORMALIZED) + "\n")
    with gzip.open(out / "bs_counts.tsv.gz", "wt") as stream:
        stream.write("\t".join(BS_COLUMNS) + "\n")
    with gzip.open(out / "native.diff.gz", "wt") as stream:
        stream.write(DIFF_HEADER)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    q = commands.add_parser("quant")
    q.add_argument("--reads", nargs="+", required=True)
    q.add_argument("--index", required=True)
    q.add_argument("--out-dir", required=True)
    q.add_argument("--run-id", required=True)
    q.add_argument("--julia", default="julia")
    q.add_argument("--whippet-bin", default="")
    q.add_argument("--flags", default="")
    q.add_argument("--identity", default="{}")
    c = commands.add_parser("compare")
    c.add_argument("--preflight", required=True)
    c.add_argument("--cache", action="append", default=[], help="run=cache_dir")
    c.add_argument("--out-dir", required=True)
    c.add_argument("--julia", default="julia")
    c.add_argument("--whippet-bin", default="")
    c.add_argument("--ids", default="{}")
    c.add_argument("--settings", default="{}")
    c.add_argument("--selection", default="{}")
    args = parser.parse_args(argv)
    if args.command == "quant":
        quant(args.reads, args.index, args.out_dir, args.run_id, args.julia, args.whippet_bin, args.flags,
              json.loads(args.identity))
    else:
        compare(args.preflight, dict(item.split("=", 1) for item in args.cache), args.out_dir, args.julia,
                args.whippet_bin, json.loads(args.ids), json.loads(args.settings), json.loads(args.selection))


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    sys.exit(main())

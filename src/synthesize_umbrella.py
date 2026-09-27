"""One report per comparison: what each tool found, and how far to trust it.

For every tool the report gives its status (ran, unsupported, missing), how
many features it tested and called, the direction of the calls (A - B), the
calling rule, which upstream data it shares with other tools, and the size
of its output. It does not combine tools into a vote: tools that share an
upstream (one HISAT2 alignment, one Whippet index, one Salmon run) are not
independent, so the shared dependency is reported instead.

It also reports, per run, the share of junction reads captured by the Whippet
index, and QC outliers computed across all runs of each group in the
comparison (median and MAD; never within a single batch shard).

Tool formats and calling rules
  deseq2      padj < 0.05                          direction: log2FoldChange
  delta       Probability >= 0.9, |DeltaPsi| >= 0.1 direction: DeltaPsi
  rmats       FDR < 0.05, |IncLevelDifference| >= 0.1, over all *.MATS.JC.txt
  leafcutter  p.adjust < 0.05 (clusters; no direction)
  suppa       p-value < 0.05, |dPSI| >= 0.1, over all *.dpsi files (direction
              as SUPPA2 reports it; see the column names in the files)
"""

import argparse
import csv
import gzip
import json
import math
import statistics
from pathlib import Path


QC_METRICS = ("hisat2_overall_alignment", "featurecounts_assigned_fraction", "junction_fragments")
OUTLIER_Z = 3.5


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _marker(path):
    """A tool that did not run leaves one comment line saying why."""
    with _open(path) as stream:
        first = stream.readline()
    if first.startswith("# inference not supported"):
        return "unsupported"
    if first.startswith("# "):
        return "skipped: " + first[2:].strip()
    return None


def _count(rows, significant, direction):
    tested = called = up = down = 0
    for row in rows:
        tested += 1
        if significant(row):
            called += 1
            value = direction(row) if direction else None
            if value is not None:
                up += value > 0
                down += value < 0
    return tested, called, (up if direction else None), (down if direction else None)


def _table(path):
    with _open(path) as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def evaluate(tool):
    """tool: {name, format, paths: [...]} -> counts and status."""
    paths = [Path(path) for path in tool["paths"]]
    if not paths or not all(path.exists() for path in paths):
        return {"status": "missing"}
    for path in paths:
        status = _marker(path) if path.is_file() else None
        if status:
            return {"status": status}
    fmt = tool["format"]
    rows = []
    if fmt == "deseq2":
        rows = _table(paths[0])
        result = _count(rows, lambda r: (_float(r["padj"]) or 1) < 0.05,
                        lambda r: _float(r["log2FoldChange"]))
        rule = "padj < 0.05"
    elif fmt == "delta":
        rows = _table(paths[0])
        result = _count(rows, lambda r: (_float(r["Probability"]) or 0) >= 0.9
                        and abs(_float(r["DeltaPsi"]) or 0) >= 0.1,
                        lambda r: _float(r["DeltaPsi"]))
        rule = "Probability >= 0.9 and |DeltaPsi| >= 0.1"
    elif fmt == "rmats":
        for path in paths:
            rows += _table(path)
        result = _count(rows, lambda r: (_float(r["FDR"]) or 1) < 0.05
                        and abs(_float(r["IncLevelDifference"]) or 0) >= 0.1,
                        lambda r: _float(r["IncLevelDifference"]))
        rule = "FDR < 0.05 and |IncLevelDifference| >= 0.1"
    elif fmt == "leafcutter":
        rows = _table(paths[0])
        result = _count(rows, lambda r: (_float(r.get("p.adjust")) or 1) < 0.05, None)
        rule = "p.adjust < 0.05"
    elif fmt == "suppa":
        for path in paths:
            with _open(path) as stream:
                header = stream.readline().rstrip("\n").split("\t")
                dpsi = [i for i, name in enumerate(header) if name.endswith("_dPSI")][0] + 1
                pval = [i for i, name in enumerate(header) if name.endswith("_p-val")][0] + 1
                rows += [line.rstrip("\n").split("\t") for line in stream]
        result = _count(rows, lambda r: (_float(r[pval]) or 1) < 0.05 and abs(_float(r[dpsi]) or 0) >= 0.1,
                        lambda r: _float(r[dpsi]))
        rule = "p-value < 0.05 and |dPSI| >= 0.1"
    else:
        raise ValueError("unknown tool format: {}".format(fmt))
    tested, called, up, down = result
    return {"status": "ran", "tested": tested, "called": called, "up": up, "down": down,
            "rule": rule, "bytes": sum(path.stat().st_size for path in paths if path.is_file())}


def qc_outliers(qc_files, groups):
    """Robust z per metric within each group, across all of its runs."""
    runs = []
    for path in qc_files:
        runs += _table(path)
    flagged = []
    by_group = {}
    for run in runs:
        by_group.setdefault(groups.get(run["run_id"], "?"), []).append(run)
    for group, members in sorted(by_group.items()):
        for metric in QC_METRICS:
            values = {run["run_id"]: _float(run[metric]) for run in members}
            values = {k: (math.log10(v + 1) if metric == "junction_fragments" else v)
                      for k, v in values.items() if v is not None}
            if len(values) < 3:
                continue
            median = statistics.median(values.values())
            mad = statistics.median(abs(v - median) for v in values.values()) * 1.4826
            if not mad:
                continue
            for run_id, value in sorted(values.items()):
                z = (value - median) / mad
                if abs(z) > OUTLIER_Z:
                    flagged.append((group, run_id, metric, round(z, 2)))
        for run in members:
            declared, inferred = run["declared_strandedness"], run["inferred_strandedness"]
            if declared not in ("auto", inferred) and inferred != "undetermined":
                flagged.append((group, run["run_id"], "strandedness declared {} inferred {}".format(
                    declared, inferred), "NA"))
    return flagged


def report(config):
    lines = ["# comparison {} ({} - {}); inference supported: {}\n".format(
        config["comparison_id"], config["group_a"], config["group_b"], config["inference_supported"])]
    for reason in config.get("reasons", []):
        lines.append("# reason: {}\n".format(reason))
    lines.append("\t".join(("tool", "readout", "family", "upstream", "status", "tested", "called",
                            "up", "down", "rule", "bytes", "output")) + "\n")
    for tool in config["tools"]:
        result = evaluate(tool)
        lines.append("\t".join(str(value) for value in (
            tool["name"], tool["readout"], tool["family"], tool["upstream"], result["status"],
            result.get("tested", "NA"), result.get("called", "NA"),
            "NA" if result.get("up") is None else result["up"],
            "NA" if result.get("down") is None else result["down"],
            result.get("rule", "NA"), result.get("bytes", "NA"), ",".join(tool["paths"]))) + "\n")
    lines.append("\n# junction reads captured by the Whippet index, per run\nrun_id\tcaptured\n")
    rates = {}
    for path in config.get("capture_rates", []):
        with open(path) as stream:
            rates.update(json.load(stream))
    for run in sorted(set(config["runs"]) & set(rates)):
        lines.append("{}\t{}\n".format(run, "NA" if rates[run] is None else rates[run]))
    lines.append("\n# QC outliers (robust z > {} within group, or strandedness mismatch)\n"
                 "group\trun_id\tmetric\tz\n".format(OUTLIER_Z))
    for row in qc_outliers(config.get("qc", []), config.get("groups", {})):
        lines.append("\t".join(map(str, row)) + "\n")
    return "".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, help="JSON written by the workflow")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    with open(args.config) as stream:
        config = json.load(stream)
    Path(args.output).write_text(report(config))


if __name__ == "__main__":
    main()

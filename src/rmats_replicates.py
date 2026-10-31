"""rMATS on biological replicates: pool technical runs, order sides by block.

rMATS models replicates: every BAM in --b1/--b2 is one replicate. Runs of
one biological replicate that were not merged at staging (library_id) would
count as separate replicates. Their BAMs are temporary and only the per-run
prep files are kept, so instead of merging BAMs the counts are pooled: post
runs with --statoff on the per-run prep files (runs listed replicate by
replicate), the comma-separated per-run counts of JC.raw.input and
JCEC.raw.input are summed within each replicate, and --task stat then fits
rMATS's model on one count per biological replicate. individualCounts files
(per run) are removed before stat, so they are never joined to pooled counts.

With a paired design (src/comparison_design.py) the replicates are listed in
block order on both sides, as --paired-stats (PAIRADISE) pairs them by
position.
"""

import json
from pathlib import Path

if __package__:
    from src.comparison_design import paired_replicates
else:
    from comparison_design import paired_replicates


COUNT_COLUMNS = ("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2")
POOLED_FILES = ("JC.raw.input.{}.txt", "JCEC.raw.input.{}.txt")
EVENT_TYPES = ("SE", "MXE", "A3SS", "A5SS", "RI")


def replicate_runs(preflight):
    """{'a': [(replicate, [runs])], 'b': [...]} in the order rMATS gets them.

    Paired designs: block order. Otherwise the preflight's replicate order."""
    members = {}
    for run, replicate in sorted(preflight["collapse"].items()):
        members.setdefault(replicate, []).append(run)
    pairs = paired_replicates(preflight)
    order = {"a": pairs[0], "b": pairs[1]} if pairs else preflight["replicates"]
    return {side: [(replicate, members[replicate]) for replicate in order[side]] for side in ("a", "b")}


def plan(preflight):
    """(runs a, runs b, sizes a, sizes b, pooled) for one comparison.

    Without technical runs and without a paired design the run order is the
    preflight's, exactly as before, and nothing is pooled."""
    grouped = replicate_runs(preflight)
    pooled = {replicate: runs for side in ("a", "b") for replicate, runs in grouped[side] if len(runs) > 1}
    if not pooled and paired_replicates(preflight) is None:
        runs = {side: list(preflight["runs"][side]) for side in ("a", "b")}
        return runs["a"], runs["b"], None, None, {}
    flat = {side: [run for _, runs in grouped[side] for run in runs] for side in ("a", "b")}
    sizes = {side: [len(runs) for _, runs in grouped[side]] for side in ("a", "b")}
    return flat["a"], flat["b"], sizes["a"], sizes["b"], pooled


def _pool(field, sizes):
    values = field.split(",")
    if len(values) != sum(sizes):
        raise ValueError("rMATS count field {!r} has {} values, expected {}".format(field, len(values), sum(sizes)))
    out, start = [], 0
    for size in sizes:
        chunk = values[start:start + size]
        start += size
        if any(value in ("NA", "") for value in chunk):
            out.append("NA" if all(value in ("NA", "") for value in chunk) else
                       str(sum(int(float(value)) for value in chunk if value not in ("NA", ""))))
        else:
            out.append(str(sum(int(float(value)) for value in chunk)))
    return ",".join(out)


def pool_counts(od, sizes_a, sizes_b):
    """Sum per-run counts into per-replicate counts in od's raw input files."""
    od = Path(od)
    sizes = {"IJC_SAMPLE_1": sizes_a, "SJC_SAMPLE_1": sizes_a, "IJC_SAMPLE_2": sizes_b, "SJC_SAMPLE_2": sizes_b}
    for event in EVENT_TYPES:
        for template in POOLED_FILES:
            path = od / template.format(event)
            if not path.exists():
                continue
            lines = path.read_text().splitlines()
            if not lines:
                continue
            header = lines[0].split("\t")
            index = {name: header.index(name) for name in COUNT_COLUMNS}
            rewritten = [lines[0]]
            for line in lines[1:]:
                values = line.split("\t")
                for name, column in index.items():
                    values[column] = _pool(values[column], sizes[name])
                rewritten.append("\t".join(values))
            path.write_text("\n".join(rewritten) + "\n")
        individual = od / "individualCounts.{}.txt".format(event)
        if individual.exists():
            individual.unlink()


def commands(base, preflight, paired):
    """The rMATS --task commands for one comparison: [(task, extra args)].

    base: the common command (b1/b2 written by the caller from plan())."""
    runs_a, runs_b, sizes_a, sizes_b, pooled = plan(preflight)
    stat = ["--paired-stats"] if paired else []
    if not pooled:
        return [base + ["--task", "inte"], base + stat + ["--task", "post"]]
    return [base + ["--task", "inte"], base + ["--statoff", "--task", "post"], ("pool", sizes_a, sizes_b),
            base + stat + ["--task", "stat"]]


def pooling_record(preflight):
    """{'pooled_replicates': {replicate: [runs]}} when runs were pooled, else {}."""
    pooled = plan(preflight)[4]
    return {"pooled_replicates": pooled} if pooled else {}


def write_status(path, preflight, model):
    """A small status.json next to the core rMATS tables (not a declared output)."""
    record = dict(tool="rmats", **pooling_record(preflight), **model)
    Path(path).write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")

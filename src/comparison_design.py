"""Paired and blocked comparison designs (the optional `design` of a comparison).

A comparison opts in with

  design: {type: paired, block: subject_id}
  design: {type: blocked, block: donor, unbalanced: drop}

`block` names a manifest metadata column (an extra column, never one of the
required ones) holding the pairing identity: donor, subject, cell line, lot.
It is not part of any run identity, so adding it reprocesses nothing.

  paired   every block has exactly one biological replicate on each side
  blocked  every block has at least one biological replicate on each side

`unbalanced` says what to do with a block present on one side only: `refuse`
(the default) stops at start-up, `drop` leaves its replicates out of the
comparison and lists them in the preflight.

The preflight of a comparison with a design gains `design` (type, block,
residual df, the block of every replicate, what was dropped) and `blocks`
(ordered by block ID, the replicates of each side). Comparisons without a
design never get these keys, so their preflights are unchanged.

Each tool records the model it used (`model`), only for designed
comparisons:

  blocked_glm          the block is a term of the tool's GLM (~ block + group)
  paired_native        the tool's own paired mode, sides given in block order
  paired_test          our own paired statistic on per-block differences
  unpaired_adaptation  the tool has no paired/blocked model; it ran its
                       ordinary two-group analysis
"""

if __package__:
    from src.umbrella_manifest import REQUIRED
else:
    from umbrella_manifest import REQUIRED


DESIGN_TYPES = ("paired", "blocked")
UNBALANCED = ("refuse", "drop")
DESIGN_KEYS = ("type", "block", "unbalanced")
# columns the manifest reader derives for library units
DERIVED = ("library_id", "library_members")

BLOCKED_GLM = {
    "deseq2_featurecounts": "DESeq2 ~ block + group",
    "deseq2_tximport": "DESeq2 ~ block + group",
    "circrna_abundance": "edgeR quasi-likelihood ~ block + group",
    "qapa": "DEXSeq full ~ sample + exon + block:exon + condition:exon, reduced without condition:exon",
}
PAIRED_TEST = {
    "dapars2": "two-sided paired t-test on per-block PDUI differences (block means, A - B)",
}
# paired modes that need design type paired (one replicate per side per block)
PAIRED_NATIVE = {
    "suppa2": "SUPPA2 diffSplice -m classical --paired (Wilcoxon signed-rank), sides in block order",
    "rmats": "rMATS --paired-stats (PAIRADISE), b1 and b2 in block order",
    "rmats_conservative": "rMATS --paired-stats (PAIRADISE), b1 and b2 in block order",
    "vasttools": "vast-tools diff --paired, sides in block order",
    "vasttools_turbo": "vast-tools diff --paired, sides in block order",
}
NO_MODEL = "tool has no paired/blocked model"


def normalise_design(comparison_id, spec):
    """The validated design of a comparison: {type, block, unbalanced}."""
    if not isinstance(spec, dict):
        raise ValueError("{}: design must be a mapping, e.g. {{type: paired, block: subject_id}}".format(
            comparison_id))
    unknown = set(spec) - set(DESIGN_KEYS)
    if unknown:
        raise ValueError("{}: unknown design key(s) {}".format(comparison_id, ", ".join(sorted(unknown))))
    kind = spec.get("type")
    if kind not in DESIGN_TYPES:
        raise ValueError("{}: design type must be paired or blocked, not {!r}".format(comparison_id, kind))
    block = spec.get("block")
    if not isinstance(block, str) or not block.strip():
        raise ValueError("{}: design needs block, the manifest column with the pairing identity".format(
            comparison_id))
    block = block.strip()
    if block in REQUIRED or block in DERIVED:
        raise ValueError("{}: design block {} is a workflow column; put the pairing identity "
                         "(donor, subject, line, lot) in an extra metadata column".format(comparison_id, block))
    unbalanced = spec.get("unbalanced", "refuse")
    if unbalanced not in UNBALANCED:
        raise ValueError("{}: design unbalanced must be refuse or drop, not {!r}".format(
            comparison_id, unbalanced))
    return {"type": kind, "block": block, "unbalanced": unbalanced}


def apply_design(comparison_id, design, runs):
    """Check a design against the selected runs ({'a': [...], 'b': [...]}).

    Returns (runs, record, blocks): the runs left after `unbalanced: drop`,
    the preflight `design` record and the ordered `blocks`. Every violation
    raises ValueError."""
    column = design["block"]
    replicate_block = {}
    for side in ("a", "b"):
        for run in runs[side]:
            value = run.metadata.get(column, "").strip()
            if not value:
                raise ValueError("{}: run {} has no {} (design block)".format(comparison_id, run.run_id, column))
            if ";" in value:
                raise ValueError("{}: library {} merges runs with different {} values ({})".format(
                    comparison_id, run.run_id, column, value))
            previous = replicate_block.setdefault(run.biological_replicate_id, value)
            if previous != value:
                raise ValueError("{}: biological replicate {} has several {} values ({}, {})".format(
                    comparison_id, run.biological_replicate_id, column, *sorted((previous, value))))
    members = {}
    for side in ("a", "b"):
        for run in runs[side]:
            members.setdefault(replicate_block[run.biological_replicate_id], {"a": set(), "b": set()})[
                side].add(run.biological_replicate_id)
    one_sided = sorted(block for block, sides in members.items() if not (sides["a"] and sides["b"]))
    dropped = sorted(replicate for block in one_sided for side in ("a", "b") for replicate in members[block][side])
    if one_sided and design["unbalanced"] == "refuse":
        raise ValueError("{}: {} block(s) have replicates on one side only ({}); fix the selection or set "
                         "design unbalanced: drop".format(comparison_id, column, ", ".join(one_sided)))
    complete = sorted(block for block in members if block not in one_sided)
    if not complete:
        raise ValueError("{}: no {} block spans both sides".format(comparison_id, column))
    if len(complete) < 2:
        raise ValueError("{}: {} complete {} block(s); a {} design needs at least 2".format(
            comparison_id, len(complete), column, design["type"]))
    if design["type"] == "paired":
        uneven = [block for block in complete if len(members[block]["a"]) != 1 or len(members[block]["b"]) != 1]
        if uneven:
            raise ValueError("{}: a paired design needs exactly one biological replicate per side in each {} "
                             "block ({} has more); use design type blocked".format(
                                 comparison_id, column, ", ".join(uneven)))
    kept = {side: [run for run in runs[side] if replicate_block[run.biological_replicate_id] in complete]
            for side in ("a", "b")}
    n_replicates = sum(len(members[block][side]) for block in complete for side in ("a", "b"))
    residual_df = n_replicates - len(complete) - 1
    if residual_df < 1:
        raise ValueError("{}: ~ {} + group leaves {} residual degrees of freedom; at least 1 is needed".format(
            comparison_id, column, residual_df))
    blocks = [{"block": block, "a": sorted(members[block]["a"]), "b": sorted(members[block]["b"])}
              for block in complete]
    record = dict(design, residual_df=residual_df, dropped_blocks=one_sided, dropped_replicates=dropped,
                  replicate_block={replicate: replicate_block[replicate]
                                   for block in blocks for side in ("a", "b") for replicate in block[side]})
    return kept, record, blocks


def has_design(preflight):
    return bool(preflight.get("design"))


def paired_replicates(preflight):
    """(A replicates, B replicates) in block order for a paired design, else None."""
    design = preflight.get("design")
    if not design or design["type"] != "paired":
        return None
    blocks = preflight["blocks"]
    return [block["a"][0] for block in blocks], [block["b"][0] for block in blocks]


def paired_runs(preflight):
    """(A runs, B runs) in block order when every replicate has one run, else None."""
    pairs = paired_replicates(preflight)
    if pairs is None:
        return None
    runs = {}
    for run, replicate in preflight["collapse"].items():
        runs.setdefault(replicate, []).append(run)
    if any(len(runs[replicate]) != 1 for side in pairs for replicate in side):
        return None
    return [runs[r][0] for r in pairs[0]], [runs[r][0] for r in pairs[1]]


def fields(model, reason):
    return {"model": model, "model_reason": reason}


def tool_model(preflight, tool):
    """{} for a comparison without a design, else {model, model_reason} for `tool`."""
    design = preflight.get("design")
    if not design:
        return {}
    if tool in BLOCKED_GLM:
        return fields("blocked_glm", BLOCKED_GLM[tool].replace("block", design["block"], 1)
                      if tool != "qapa" else BLOCKED_GLM[tool])
    if tool in PAIRED_TEST:
        return fields("paired_test", PAIRED_TEST[tool])
    if tool in PAIRED_NATIVE:
        if design["type"] != "paired":
            return fields("unpaired_adaptation", "its paired mode needs design type paired "
                          "(one replicate per side per block); ran the ordinary two-group analysis")
        return fields("paired_native", PAIRED_NATIVE[tool])
    return fields("unpaired_adaptation", NO_MODEL)


def block_differences(preflight, values):
    """Per-block differences mean(A) - mean(B) of {replicate: value} (None values skipped).

    Blocks without a value on both sides are left out."""
    differences = []
    for block in preflight["blocks"]:
        sides = [[values[r] for r in block[side] if values.get(r) is not None] for side in ("a", "b")]
        if sides[0] and sides[1]:
            differences.append(sum(sides[0]) / len(sides[0]) - sum(sides[1]) / len(sides[1]))
    return differences


def wilcoxon_signed_rank(differences, exact_limit=50):
    """Two-sided Wilcoxon signed-rank p-value of paired differences.

    Zero differences are dropped (Wilcoxon's convention) and tied magnitudes
    get mid-ranks. Exact (sign-flip distribution by dynamic programming over
    doubled ranks) up to `exact_limit` non-zero differences, else the normal
    approximation with tie and continuity corrections. None below 2 values."""
    import math
    values = [d for d in differences if d is not None and d != 0]
    n = len(values)
    if n < 2:
        return None
    magnitudes = sorted(abs(d) for d in values)
    rank = {}
    position = 0
    while position < n:
        end = position
        while end + 1 < n and magnitudes[end + 1] == magnitudes[position]:
            end += 1
        rank[magnitudes[position]] = (position + end) / 2.0 + 1
        position = end + 1
    doubled = [int(round(2 * rank[abs(d)])) for d in values]
    total = sum(doubled)
    observed = sum(r for r, d in zip(doubled, values) if d > 0)
    deviation = abs(observed - total / 2.0)
    if n <= exact_limit:
        counts = [0] * (total + 1)
        counts[0] = 1
        for r in doubled:
            for s in range(total, r - 1, -1):
                counts[s] += counts[s - r]
        extreme = sum(c for s, c in enumerate(counts) if c and abs(s - total / 2.0) >= deviation - 1e-9)
        return min(1.0, extreme / float(2 ** n))
    ties = {}
    for d in values:
        ties[abs(d)] = ties.get(abs(d), 0) + 1
    variance = n * (n + 1) * (2 * n + 1) / 24.0 - sum(t ** 3 - t for t in ties.values()) / 48.0
    if variance <= 0:
        return 1.0
    z = max(0.0, deviation / 2.0 - 0.5) / math.sqrt(variance)
    return min(1.0, math.erfc(z / math.sqrt(2)))


WEIGHTING_REASON = ("per-sample posteriors are pooled; replicates with more inputs weigh more "
                    "(each run is one input; the computation is unchanged)")
# delta tools that take one input per run and pool per-sample posteriors per side
POOLED_POSTERIOR_TOOLS = ("whippet_delta", "microexonator_delta", "microexonator_whippet_delta", "whippet_circ")


def input_weighting(preflight):
    """{} unless a replicate has several runs; else the inputs per replicate.

    weighting: {replicate: inputs} for replicates with more than one input;
    weights: {side: {replicate: inputs / inputs on that side}} for every
    replicate. For tools that pool per-run posteriors (Whippet delta, the
    MicroExonator delta port, Whippet --circ), this is each replicate's share
    of its side."""
    inputs = {}
    for run, replicate in preflight["collapse"].items():
        inputs[replicate] = inputs.get(replicate, 0) + 1
    several = {replicate: n for replicate, n in sorted(inputs.items()) if n > 1}
    if not several:
        return {}
    weights = {}
    for side in ("a", "b"):
        total = sum(inputs[r] for r in preflight["replicates"][side])
        weights[side] = {r: round(inputs[r] / float(total), 6) for r in preflight["replicates"][side]}
    return {"weighting": several, "weights": weights, "weighting_reason": WEIGHTING_REASON}

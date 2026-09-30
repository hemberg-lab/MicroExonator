"""Check each two-group comparison before any statistics are run.

The comparisons file is YAML (JSON also works) with a list under
`comparisons`. Each comparison has a stable `comparison_id` (used in output
paths), a `project_id` (comparisons stay within one project), two sides A and
B (effects are reported as A - B), and optionally `exclude`, a list of run_id
or biological_replicate_id to leave out.

A side is either a manifest group, the storage label its shards are filed
under:

  group_a: case
  group_b: ctrl

or a selection of runs, with a `label` used in reports:

  a: {label: asd_female, where: {phenotype: ASD, sex: F}}
  b: {label: ctrl_female, groups: [ctrl], where: {sex: F}}

A selection keeps the runs that match every key it gives:
  groups   manifest group(s)
  samples  sample_id(s)
  runs     run_id(s)
  where    {column: value or [values]} on any manifest column, including
           extra metadata columns
New conditions are therefore new metadata columns plus new comparisons: the
runs stay in the shards of their original group and nothing is reprocessed.

For each comparison this writes a JSON preflight with the runs, the
run-to-replicate collapse map, the group-by-batch shards to join, and
`inference_supported`. Inference is not supported, and only descriptive
outputs are made, when a group has fewer than two biological replicates, or
when group is confounded with batch, layout or strandedness (the two groups
share no value of that factor). No covariates are ever added silently.
rMATS additionally needs every run of the comparison to have one layout.
"""

import argparse
import json
from pathlib import Path

if __package__:
    from src.umbrella_manifest import SAFE_ID, load_umbrella_manifest
else:
    from umbrella_manifest import SAFE_ID, load_umbrella_manifest


FACTORS = ("batch_id", "layout", "strandedness")
SELECTORS = ("groups", "samples", "runs", "where")


def _as_list(value, what):
    values = value if isinstance(value, list) else [value]
    if not values:
        raise ValueError("{} is empty".format(what))
    return [_as_text(item) for item in values]


def _as_text(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list)) or value is None:
        raise ValueError("selection values must be plain values, not {!r}".format(value))
    return str(value).strip()


def _side(comparison, side):
    """The normalised selection of one side: {label, groups?, samples?, runs?, where?}."""
    shorthand = comparison.get("group_" + side)
    spec = comparison.get(side)
    if (shorthand is None) == (spec is None):
        raise ValueError("{}: give either group_{} or {} (a selection)".format(
            comparison.get("comparison_id"), side, side))
    if spec is None:
        if not isinstance(shorthand, str) or not SAFE_ID.fullmatch(shorthand):
            raise ValueError("comparison needs a safe group_{}: {!r}".format(side, shorthand))
        return {"label": shorthand, "groups": [shorthand]}
    if not isinstance(spec, dict):
        raise ValueError("{}: side {} must be a mapping".format(comparison.get("comparison_id"), side))
    unknown = set(spec) - set(SELECTORS) - {"label"}
    if unknown:
        raise ValueError("{}: unknown selection key(s) {}".format(
            comparison.get("comparison_id"), ", ".join(sorted(unknown))))
    label = spec.get("label")
    if not isinstance(label, str) or not SAFE_ID.fullmatch(label):
        raise ValueError("{}: side {} needs a safe label: {!r}".format(
            comparison.get("comparison_id"), side, label))
    normal = {"label": label}
    for key in ("groups", "samples", "runs"):
        if key in spec:
            normal[key] = sorted(set(_as_list(spec[key], key)))
    if "where" in spec:
        if not isinstance(spec["where"], dict) or not spec["where"]:
            raise ValueError("{}: where must be a non-empty mapping".format(comparison.get("comparison_id")))
        normal["where"] = {str(column): sorted(set(_as_list(value, column)))
                           for column, value in sorted(spec["where"].items())}
    if len(normal) == 1:
        raise ValueError("{}: side {} selects nothing; give groups, samples, runs or where".format(
            comparison.get("comparison_id"), side))
    return normal


def load_comparisons(path):
    text = Path(path).read_text()
    try:
        import yaml
        data = yaml.safe_load(text)
    except ImportError:
        data = json.loads(text)
    comparisons = data.get("comparisons") if isinstance(data, dict) else None
    if not isinstance(comparisons, list) or not comparisons:
        raise ValueError("comparisons file needs a non-empty 'comparisons' list")
    seen = set()
    for comparison in comparisons:
        for field in ("comparison_id", "project_id"):
            value = comparison.get(field)
            if not isinstance(value, str) or not SAFE_ID.fullmatch(value):
                raise ValueError("comparison needs a safe {}: {!r}".format(field, value))
        if comparison["comparison_id"] in seen:
            raise ValueError("duplicate comparison_id: {}".format(comparison["comparison_id"]))
        seen.add(comparison["comparison_id"])
        sides = {side: _side(comparison, side) for side in ("a", "b")}
        if sides["a"] == dict(sides["b"], label=sides["a"]["label"]):
            raise ValueError("{}: sides A and B select the same runs".format(comparison["comparison_id"]))
        if sides["a"]["label"] == sides["b"]["label"]:
            raise ValueError("{}: group_a and group_b are the same".format(comparison["comparison_id"]))
        comparison["selection"] = sides
        comparison["group_a"] = sides["a"]["label"]
        comparison["group_b"] = sides["b"]["label"]
        comparison.setdefault("exclude", [])
    return comparisons


def _matches(run, selection):
    if "groups" in selection and run.group not in selection["groups"]:
        return False
    if "samples" in selection and run.sample_id not in selection["samples"]:
        return False
    if "runs" in selection and run.run_id not in selection["runs"]:
        return False
    for column, values in selection.get("where", {}).items():
        if run.metadata.get(column, "").strip() not in values:
            return False
    return True


def preflight(manifest, comparison):
    comparison_id = comparison["comparison_id"]
    selection = comparison.get("selection") or {
        side: _side(comparison, side) for side in ("a", "b")}
    columns = set().union(*(run.metadata.keys() for run in manifest.runs)) if manifest.runs else set()
    for side in ("a", "b"):
        missing = set(selection[side].get("where", {})) - columns
        if missing:
            raise ValueError("{}: no manifest column {}".format(comparison_id, ", ".join(sorted(missing))))
    excluded = set(comparison.get("exclude", []))
    runs = {"a": [], "b": []}
    for run in manifest.included_runs():
        if run.project_id != comparison["project_id"]:
            continue
        if run.run_id in excluded or run.biological_replicate_id in excluded:
            continue
        for side in ("a", "b"):
            if _matches(run, selection[side]):
                runs[side].append(run)
    reasons = []
    for side in ("a", "b"):
        if not runs[side]:
            raise ValueError("{}: no included runs for group {}".format(
                comparison_id, selection[side]["label"]))
    both = {run.run_id for run in runs["a"]} & {run.run_id for run in runs["b"]}
    if both:
        raise ValueError("{}: runs selected on both sides: {}".format(comparison_id, ", ".join(sorted(both))))
    split = ({run.biological_replicate_id for run in runs["a"]}
             & {run.biological_replicate_id for run in runs["b"]})
    if split:
        raise ValueError("{}: biological replicates split across sides: {}".format(
            comparison_id, ", ".join(sorted(split))))
    everything = runs["a"] + runs["b"]
    references = {run.reference_id for run in everything}
    if len(references) != 1:
        raise ValueError("{}: runs span several reference bundles".format(comparison_id))
    replicates = {side: sorted({run.biological_replicate_id for run in runs[side]})
                  for side in ("a", "b")}
    for side in ("a", "b"):
        if len(replicates[side]) < 2:
            reasons.append("group {} has {} biological replicate(s); at least 2 are needed".format(
                selection[side]["label"], len(replicates[side])))
    for factor in FACTORS:
        values = {side: {getattr(run, factor) for run in runs[side]} for side in ("a", "b")}
        if not values["a"] & values["b"]:
            reasons.append("group is confounded with {} ({} vs {})".format(
                factor, ",".join(sorted(values["a"])), ",".join(sorted(values["b"]))))
    technical = sorted(replicate for replicate in set(
        run.biological_replicate_id for run in everything)
        if sum(run.biological_replicate_id == replicate for run in everything) > 1)
    layouts = {run.layout for run in everything}
    rmats_reasons = list(reasons)
    if len(layouts) != 1:
        rmats_reasons.append("mixed layouts ({}); rMATS needs one".format(",".join(sorted(layouts))))
    # the shards that store the selected runs (their original group and batch)
    shards = sorted({(run.reference_id, run.project_id, run.group, run.batch_id) for run in everything})
    return {
        "comparison_id": comparison_id,
        "project_id": comparison["project_id"],
        "reference_id": references.pop(),
        "group_a": selection["a"]["label"],
        "group_b": selection["b"]["label"],
        "selection": selection,
        "effect": "group_a - group_b",
        "excluded": sorted(excluded),
        "runs": {side: sorted(run.run_id for run in runs[side]) for side in ("a", "b")},
        "replicates": replicates,
        "collapse": {run.run_id: run.biological_replicate_id for run in sorted(everything, key=lambda r: r.run_id)},
        "layouts": sorted(layouts),
        # Count matrices are summed per replicate; per-run PSI inputs are not
        # combined, so delta tools see these replicates' runs separately.
        "technical_run_replicates": technical,
        "shards": ["/".join(shard) for shard in shards],
        "inference_supported": not reasons,
        "reasons": reasons,
        "tools": {"rmats": {"supported": not rmats_reasons, "reasons": rmats_reasons}},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--comparisons", required=True)
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    manifest = load_umbrella_manifest(args.manifest)
    comparisons = {c["comparison_id"]: c for c in load_comparisons(args.comparisons)}
    if args.comparison_id not in comparisons:
        parser.error("unknown comparison_id {}".format(args.comparison_id))
    result = preflight(manifest, comparisons[args.comparison_id])
    with open(args.output, "w") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()

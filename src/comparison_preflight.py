"""Check each two-group comparison before any statistics are run.

The comparisons file is YAML (JSON also works) with a list under
`comparisons`, each with:

  comparison_id  stable ID used in output paths
  project_id     comparisons stay within one project
  group_a        e.g. the case group
  group_b        e.g. the control group; effects are reported as A - B
  exclude        optional list of run_id or biological_replicate_id to leave out

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
        for field in ("comparison_id", "project_id", "group_a", "group_b"):
            value = comparison.get(field)
            if not isinstance(value, str) or not SAFE_ID.fullmatch(value):
                raise ValueError("comparison needs a safe {}: {!r}".format(field, value))
        if comparison["comparison_id"] in seen:
            raise ValueError("duplicate comparison_id: {}".format(comparison["comparison_id"]))
        seen.add(comparison["comparison_id"])
        if comparison["group_a"] == comparison["group_b"]:
            raise ValueError("{}: group_a and group_b are the same".format(comparison["comparison_id"]))
        comparison.setdefault("exclude", [])
    return comparisons


def preflight(manifest, comparison):
    excluded = set(comparison.get("exclude", []))
    runs = {"a": [], "b": []}
    for run in manifest.included_runs():
        if run.project_id != comparison["project_id"]:
            continue
        if run.run_id in excluded or run.biological_replicate_id in excluded:
            continue
        for side in ("a", "b"):
            if run.group == comparison["group_" + side]:
                runs[side].append(run)
    reasons = []
    for side in ("a", "b"):
        if not runs[side]:
            raise ValueError("{}: no included runs for group {}".format(
                comparison["comparison_id"], comparison["group_" + side]))
    everything = runs["a"] + runs["b"]
    references = {run.reference_id for run in everything}
    if len(references) != 1:
        raise ValueError("{}: runs span several reference bundles".format(comparison["comparison_id"]))
    replicates = {side: sorted({run.biological_replicate_id for run in runs[side]})
                  for side in ("a", "b")}
    for side in ("a", "b"):
        if len(replicates[side]) < 2:
            reasons.append("group {} has {} biological replicate(s); at least 2 are needed".format(
                comparison["group_" + side], len(replicates[side])))
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
    shards = sorted({(run.reference_id, run.project_id, run.group, run.batch_id) for run in everything})
    return {
        "comparison_id": comparison["comparison_id"],
        "project_id": comparison["project_id"],
        "reference_id": references.pop(),
        "group_a": comparison["group_a"],
        "group_b": comparison["group_b"],
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

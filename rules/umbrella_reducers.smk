"""Durable splicing shards, comparison preflight, and joined comparison inputs.

MicroExonator and Whippet per-run outputs are copied into group-by-batch
shards (src/reduce_splicing.py). The per-run files are left in place: they
stay the inputs of the existing delta rules until rebuilding from shards has
been checked on real data.

With `umbrella_comparisons` set (YAML, see src/comparison_preflight.py), each
comparison gets a preflight JSON and one joined matrix per kind under
comparisons/{reference_id}/{project_id}/{comparison_id}/. Count matrices are
summed per biological replicate; the others keep one column block per run.
"""

import hashlib
import json

from src.comparison_preflight import load_comparisons, preflight as comparison_preflight


rule umbrella_microexonator_shard:
    input:
        tables=lambda w: [downstream_PSI(record.run_id) for record in umbrella_shard_runs(w)]
    output:
        shard=protected("splicing/" + UMBRELLA_SHARD + ".microexonator.tsv.gz"),
        checksums=protected("splicing/" + UMBRELLA_SHARD + ".microexonator_checksums.json")
    params:
        inputs=lambda w: " ".join("--input {}={}".format(record.run_id, downstream_PSI(record.run_id))
                                  for record in umbrella_shard_runs(w)),
        prefix="splicing/" + UMBRELLA_SHARD,
        manifest_sha256=umbrella_shard_sha
    run:
        shell("python3 src/reduce_splicing.py microexonator --prefix {params.prefix} "
              "--reference-id {wildcards.reference_id} "
              "--manifest-sha256 {params.manifest_sha256} {params.inputs}")


rule umbrella_whippet_shard:
    input:
        psi=umbrella_run_files("/whippet/quant.psi.gz")
    output:
        shard=protected("splicing/" + UMBRELLA_SHARD + ".whippet.tsv.gz"),
        checksums=protected("splicing/" + UMBRELLA_SHARD + ".whippet_checksums.json")
    params:
        inputs=lambda w: " ".join("--input {}={}/whippet/quant.psi.gz".format(record.run_id, record.work_dir)
                                  for record in umbrella_shard_runs(w)),
        prefix="splicing/" + UMBRELLA_SHARD,
        manifest_sha256=umbrella_shard_sha
    run:
        shell("python3 src/reduce_splicing.py whippet --prefix {params.prefix} "
              "--reference-id {wildcards.reference_id} "
              "--manifest-sha256 {params.manifest_sha256} {params.inputs}")


UMBRELLA_SPLICING_SHARDS = {"microexonator": [], "whippet": []}
for reference_id, project_id, group, batch_id in sorted({
        (run.reference_id, run.project_id, run.group, run.batch_id)
        for run in UMBRELLA_MANIFEST.included_runs()}):
    shard = "{}/{}/{}/{}".format(reference_id, project_id, group, batch_id)
    for kind in UMBRELLA_SPLICING_SHARDS:
        UMBRELLA_SPLICING_SHARDS[kind].append("splicing/{}.{}.tsv.gz".format(shard, kind))


# ------------------------------------------------------------ comparisons

# kind -> (shard path pattern, joiner kind, collapse technical runs by summing)
UMBRELLA_JOIN_KINDS = {
    "junctions": ("junctions/{}.junctions.tsv.gz", "junctions", True),
    "featurecounts": ("genes/{}.featurecounts.tsv.gz", "featurecounts", True),
    "salmon_counts": ("genes/{}.salmon_counts.tsv.gz", "salmon", True),
    "salmon_tpm": ("genes/{}.salmon_tpm.tsv.gz", "salmon", False),
    "salmon_length": ("genes/{}.salmon_length.tsv.gz", "salmon", False),
    "salmon_tx_counts": ("genes/{}.salmon_tx_counts.tsv.gz", "salmon", True),
    "salmon_tx_tpm": ("genes/{}.salmon_tx_tpm.tsv.gz", "salmon", False),
    "microexonator": ("splicing/{}.microexonator.tsv.gz", "microexonator", False),
    "whippet": ("splicing/{}.whippet.tsv.gz", "whippet", False),
}

UMBRELLA_COMPARISONS = {}
if config.get("umbrella_comparisons"):
    try:
        for comparison in load_comparisons(config["umbrella_comparisons"]):
            UMBRELLA_COMPARISONS[comparison["comparison_id"]] = comparison_preflight(
                UMBRELLA_MANIFEST, comparison)
    except (OSError, ValueError) as error:
        raise WorkflowError("invalid umbrella_comparisons: {}".format(error))
UMBRELLA_COMPARISON_ROOT = "comparisons/{reference_id}/{project_id}/{comparison_id}"


def umbrella_comparison(wildcards):
    result = UMBRELLA_COMPARISONS.get(wildcards.comparison_id)
    if result is None or (result["reference_id"], result["project_id"]) != (
            wildcards.reference_id, wildcards.project_id):
        raise WorkflowError("unknown umbrella comparison {}".format(wildcards.comparison_id))
    return result


rule umbrella_comparison_preflight:
    output:
        UMBRELLA_COMPARISON_ROOT + "/preflight.json"
    params:
        membership_sha256=lambda w: hashlib.sha256(json.dumps(
            umbrella_comparison(w), sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    run:
        result = umbrella_comparison(wildcards)
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")


rule umbrella_comparison_join:
    input:
        shards=lambda w: [UMBRELLA_JOIN_KINDS[w.kind][0].format(shard)
                          for shard in umbrella_comparison(w)["shards"]],
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json"
    output:
        joined=UMBRELLA_COMPARISON_ROOT + "/joined/{kind}.tsv.gz",
        provenance=UMBRELLA_COMPARISON_ROOT + "/joined/{kind}.provenance.json"
    wildcard_constraints:
        kind="|".join(UMBRELLA_JOIN_KINDS)
    params:
        joiner=lambda w: UMBRELLA_JOIN_KINDS[w.kind][1],
        collapse=lambda w: ("--collapse " + UMBRELLA_COMPARISON_ROOT.format(**w) + "/preflight.json")
            if UMBRELLA_JOIN_KINDS[w.kind][2] else ""
    shell:
        "python3 src/join_umbrella_shards.py --kind {params.joiner} "
        "--reference-id {wildcards.reference_id} --output {output.joined} "
        "--provenance {output.provenance} --select {input.preflight} "
        "{params.collapse} {input.shards}"


UMBRELLA_COMPARISON_TARGETS = []
for comparison_id, result in sorted(UMBRELLA_COMPARISONS.items()):
    root = "comparisons/{}/{}/{}".format(result["reference_id"], result["project_id"], comparison_id)
    UMBRELLA_COMPARISON_TARGETS += [root + "/joined/{}.tsv.gz".format(kind) for kind in UMBRELLA_JOIN_KINDS]

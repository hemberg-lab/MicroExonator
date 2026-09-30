"""Per-comparison differential analyses from joined shards, and one synthesis.

Every tool writes fixed outputs. When a tool cannot run (inference not
supported by the preflight, rMATS with mixed layouts) its output is a single comment line saying why, and no
inferential command is run.

  deseq2_featurecounts   DESeq2 ~ group on featureCounts gene counts
  deseq2_tximport        DESeq2 ~ group on Salmon via DESeqDataSetFromTximport
  microexonator_delta    src/me_delta.py on per-run tables rebuilt from shards
  whippet_delta          whippet-delta.jl on per-run psi files rebuilt from shards
  rmats                  rMATS-turbo inte + post on the kept prep files
  leafcutter             conda-managed Python LeafCutter clustering and DS
  suppa2                 (pilot) SUPPA2 diffSplice on Salmon transcript TPM
"""

import json


UMBRELLA_RMATS_TYPES = ("SE", "A5SS", "A3SS", "MXE", "RI")
UMBRELLA_SUPPA_TYPES = ("SE", "A5", "A3", "MX", "RI", "AF", "AL")
UMBRELLA_SUPPA_EVENTS = UMBRELLA_REFERENCE_ROOT + "/suppa/events"
UMBRELLA_UNSUPPORTED = "# inference not supported"


def umbrella_runs_of(wildcards, side):
    return umbrella_comparison(wildcards)["runs"][side]


def umbrella_write_marker(paths, text):
    import gzip as _gzip
    for path in paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        opener = _gzip.open if str(path).endswith(".gz") else open
        with opener(path, "wt") as stream:
            stream.write(text + "\n")


rule umbrella_deseq2:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        matrices=lambda w: [UMBRELLA_COMPARISON_ROOT.format(**w) + "/joined/{}.tsv.gz".format(kind)
                            for kind in (("featurecounts",) if w.route == "featurecounts"
                                         else ("salmon_counts", "salmon_tpm", "salmon_length"))]
    output:
        results=UMBRELLA_COMPARISON_ROOT + "/deseq2_{route}.results.tsv.gz",
        descriptive=UMBRELLA_COMPARISON_ROOT + "/deseq2_{route}.descriptive.tsv.gz",
        status=UMBRELLA_COMPARISON_ROOT + "/deseq2_{route}.status.txt"
    wildcard_constraints:
        route="featurecounts|tximport"
    params:
        prefix=UMBRELLA_COMPARISON_ROOT + "/deseq2_{route}",
        inputs=lambda w, input: ("--counts {}".format(input.matrices[0]) if w.route == "featurecounts"
                                 else "--counts {} --tpm {} --length {}".format(*input.matrices))
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/deseq2_{route}.log"
    conda:
        "../envs/umbrella-comparisons.yaml"
    shell:
        "Rscript src/expression_deseq2.R --route {wildcards.route} --preflight {input.preflight} "
        "{params.inputs} --out {params.prefix} > {log} 2>&1"


rule umbrella_microexonator_delta:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        joined=UMBRELLA_COMPARISON_ROOT + "/joined/microexonator.tsv.gz",
        microexons=UMBRELLA_FIXED_MICROEXONS
    output:
        UMBRELLA_COMPARISON_ROOT + "/microexonator_delta.tsv"
    params:
        work=UMBRELLA_COMPARISON_ROOT + "/delta_inputs/microexonator",
        options=config.get("umbrella_me_delta_options", "")
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/microexonator_delta.log"
    conda:
        "../envs/umbrella-delta.yaml"
    shell:
        # --options=: an empty value must stay attached to its flag (Snakemake
        # writes nothing for an empty {params.options:q})
        "python3 src/run_comparison_tools.py microexonator "
        "--preflight {input.preflight:q} --joined {input.joined:q} "
        "--microexons {input.microexons:q} --work {params.work:q} "
        "--options={params.options:q} --outputs {output:q} --log {log:q} 2> {log:q}"


# Legacy MegaSearch route, optional (umbrella_optional: microexonator_whippet_delta):
# Whippet delta on Whippet PSI with MicroExonator's corrected PSI swapped into
# the microexon nodes. Runs next to the Whippet-free microexonator_delta, to
# validate it.
rule umbrella_microexonator_whippet_delta:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        joined=UMBRELLA_COMPARISON_ROOT + "/joined/microexonator.tsv.gz",
        joined_whippet=UMBRELLA_COMPARISON_ROOT + "/joined/whippet.tsv.gz",
        exons=UMBRELLA_WHIPPET_MEMBERS[1],
        microexons=UMBRELLA_FIXED_MICROEXONS
    output:
        full=UMBRELLA_COMPARISON_ROOT + "/microexonator_whippet_delta.diff.gz",
        microexons=UMBRELLA_COMPARISON_ROOT + "/microexonator_whippet_delta.microexons.tsv"
    params:
        work=UMBRELLA_COMPARISON_ROOT + "/delta_inputs/microexonator_whippet",
        julia=config.get("julia", "julia"),
        whippet_bin=config.get("whippet_bin_folder", "")
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/microexonator_whippet_delta.log"
    shell:
        "python3 src/run_comparison_tools.py microexonator_whippet "
        "--preflight {input.preflight:q} --joined {input.joined:q} "
        "--joined-whippet {input.joined_whippet:q} --whippet-exons {input.exons:q} "
        "--microexons {input.microexons:q} --work {params.work:q} "
        "--julia {params.julia:q} --whippet-bin={params.whippet_bin:q} "
        "--outputs {output.full:q} {output.microexons:q} --log {log:q} 2>> {log:q}"


rule umbrella_whippet_delta:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        joined=UMBRELLA_COMPARISON_ROOT + "/joined/whippet.tsv.gz"
    output:
        UMBRELLA_COMPARISON_ROOT + "/whippet_delta.diff.gz"
    params:
        work=UMBRELLA_COMPARISON_ROOT + "/delta_inputs/whippet",
        prefix=UMBRELLA_COMPARISON_ROOT + "/whippet_delta",
        julia=config.get("julia", "julia"),
        whippet_bin=config.get("whippet_bin_folder", "")
    run:
        result = umbrella_comparison(wildcards)
        if not result["inference_supported"]:
            umbrella_write_marker(output, UMBRELLA_UNSUPPORTED)
        else:
            paths = {side: ["{}/{}.psi.gz".format(params.work, run) for run in result["runs"][side]]
                     for side in ("a", "b")}
            pairs = " ".join("--run {}={}".format(run, path) for side in ("a", "b")
                             for run, path in zip(result["runs"][side], paths[side]))
            shell("python3 src/umbrella_tool_inputs.py delta --kind whippet --joined {input.joined} " + pairs)
            shell("{params.julia} {params.whippet_bin}/whippet-delta.jl -a " + ",".join(paths["a"])
                  + " -b " + ",".join(paths["b"]) + " -o {params.prefix}")
            shell("rm -rf {params.work}")


def umbrella_rmats_records(wildcards):
    result = umbrella_comparison(wildcards)
    records = []
    for side in ("a", "b"):
        for run_id in result["runs"][side]:
            record = UMBRELLA_MANIFEST.by_run[run_id]
            records.append("rmats/{}/{}/{}/{}/{}.rmats".format(
                record.reference_id, record.project_id, record.group, record.batch_id, run_id))
    return records


rule umbrella_rmats_post:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        prep=umbrella_rmats_records,
        gtf=UMBRELLA_RMATS_GTF
    output:
        [UMBRELLA_COMPARISON_ROOT + "/rmats/{}.MATS.JC.txt".format(kind) for kind in UMBRELLA_RMATS_TYPES]
    params:
        work=UMBRELLA_COMPARISON_ROOT + "/rmats_work"
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/rmats.log"
    threads: 4
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "python3 src/run_comparison_tools.py rmats --preflight {input.preflight} --work {params.work} "
        "--log {log} --gtf {input.gtf} --threads {threads} --prep {input.prep} --outputs {output}"


rule umbrella_leafcutter:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        joined=UMBRELLA_COMPARISON_ROOT + "/joined/junctions.tsv.gz"
    output:
        UMBRELLA_COMPARISON_ROOT + "/leafcutter_cluster_significance.txt"
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/leafcutter.log"
    threads: 4
    conda:
        "../envs/umbrella-leafcutter.yaml"
    shell:
        "python3 src/run_leafcutter.py --preflight {input.preflight} "
        "--joined {input.joined} --output {output} --log {log} --threads {threads}"


rule umbrella_suppa_events:
    input:
        UMBRELLA_SALMON_GTF
    output:
        [UMBRELLA_SUPPA_EVENTS + "_{}_strict.ioe".format(kind) for kind in UMBRELLA_SUPPA_TYPES]
    params:
        prefix=UMBRELLA_SUPPA_EVENTS
    conda:
        "../envs/umbrella-suppa.yaml"
    shell:
        "suppa.py generateEvents -i {input} -o {params.prefix} -f ioe -e SE SS MX RI FL"


rule umbrella_suppa:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        joined=UMBRELLA_COMPARISON_ROOT + "/joined/salmon_tx_tpm.tsv.gz",
        events=[UMBRELLA_SUPPA_EVENTS + "_{}_strict.ioe".format(kind) for kind in UMBRELLA_SUPPA_TYPES]
    output:
        [UMBRELLA_COMPARISON_ROOT + "/suppa2/{}.dpsi".format(kind) for kind in UMBRELLA_SUPPA_TYPES]
    params:
        work=UMBRELLA_COMPARISON_ROOT + "/suppa_work"
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/suppa2.log"
    conda:
        "../envs/umbrella-suppa.yaml"
    shell:
        "python3 src/run_comparison_tools.py suppa --preflight {input.preflight} --work {params.work} "
        "--log {log} --joined {input.joined} --events {input.events} --outputs {output}"


def umbrella_tool_table(wildcards):
    root = UMBRELLA_COMPARISON_ROOT.format(**wildcards)
    tools = [
        ("deseq2_featurecounts", "expression", "annotation_based", "hisat2", "deseq2",
         [root + "/deseq2_featurecounts.results.tsv.gz"]),
        ("deseq2_tximport", "expression", "annotation_based", "salmon", "deseq2",
         [root + "/deseq2_tximport.results.tsv.gz"]),
        ("microexonator_delta", "splicing", "annotation_based", "microexonator", "delta",
         [root + "/microexonator_delta.tsv"]),
        ("whippet_delta", "splicing", "annotation_based", "whippet", "delta",
         [root + "/whippet_delta.diff.gz"]),
        ("rmats", "splicing", "annotation_based", "hisat2", "rmats",
         [root + "/rmats/{}.MATS.JC.txt".format(kind) for kind in UMBRELLA_RMATS_TYPES]),
    ]
    if UMBRELLA_OPTIONAL.get("microexonator_whippet_delta", False):
        tools.append(("microexonator_whippet_delta", "splicing", "annotation_based",
                      "microexonator+whippet", "delta",
                      [root + "/microexonator_whippet_delta.microexons.tsv"]))
    if UMBRELLA_OPTIONAL.get("leafcutter", True):
        tools.append(("leafcutter", "splicing", "annotation_free", "hisat2", "leafcutter",
                      [root + "/leafcutter_cluster_significance.txt"]))
    if UMBRELLA_OPTIONAL.get("suppa2", True):
        tools.append(("suppa2", "splicing", "annotation_based", "salmon", "suppa",
                      [root + "/suppa2/{}.dpsi".format(kind) for kind in UMBRELLA_SUPPA_TYPES]))
    if not UMBRELLA_OPTIONAL.get("rmats"):
        tools = [tool for tool in tools if tool[0] != "rmats"]
    return tools


rule umbrella_synthesis:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        outputs=lambda w: [path for tool in umbrella_tool_table(w) for path in tool[5]],
        capture=lambda w: ["junctions/{}.capture_rates.json".format(shard)
                           for shard in umbrella_comparison(w)["shards"]],
        qc=lambda w: ["qc/{}.qc.tsv.gz".format(shard) for shard in umbrella_comparison(w)["shards"]]
    output:
        report=UMBRELLA_COMPARISON_ROOT + "/synthesis.tsv",
        config=UMBRELLA_COMPARISON_ROOT + "/synthesis_config.json"
    run:
        result = umbrella_comparison(wildcards)
        config_out = dict(
            comparison_id=result["comparison_id"], group_a=result["group_a"],
            group_b=result["group_b"], inference_supported=result["inference_supported"],
            reasons=result["reasons"], runs=result["runs"]["a"] + result["runs"]["b"],
            groups={run: UMBRELLA_MANIFEST.by_run[run].group for run in result["collapse"]},
            capture_rates=list(input.capture), qc=list(input.qc),
            tools=[dict(name=name, readout=readout, family=family, upstream=upstream,
                        format=fmt, paths=paths)
                   for name, readout, family, upstream, fmt, paths in umbrella_tool_table(wildcards)])
        Path(output.config).write_text(json.dumps(config_out, sort_keys=True, indent=2) + "\n")
        shell("python3 src/synthesize_umbrella.py --config {output.config} --output {output.report}")


UMBRELLA_SYNTHESIS_TARGETS = [
    "comparisons/{}/{}/{}/synthesis.tsv".format(result["reference_id"], result["project_id"], comparison_id)
    for comparison_id, result in sorted(UMBRELLA_COMPARISONS.items())]

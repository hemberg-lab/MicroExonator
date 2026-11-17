"""Whippet BS analysis from retained core --circ outputs by default.

Select it with umbrella_modules: [whippet_circ]. Its scope is "selected" by default, like
the circRNA modules. Per run, import the checked core quantification into a cache.
whippet_circ_source: parallel explicitly retains the former independent quant.
Per comparison: whippet-delta on those psi files (native.diff.gz, linear
events only) and whippet-delta's model on the BS events (results.tsv.gz). See
src/run_whippet_circ.py.

Config: whippet_circ_min_reads (5), whippet_circ_min_samples (1), whippet_circ_draws (1000),
whippet_circ_seed (123456); these change only the comparison outputs.
"""

WHIPPET_CIRC = "whippet_circ"
WHIPPET_CIRC_JULIA = config.get("julia", "julia")
WHIPPET_CIRC_BIN = config.get("whippet_bin_folder", "")
MODULE_SETTINGS[WHIPPET_CIRC] = {"min_reads": int(config.get("whippet_circ_min_reads", 5)),
                                 "min_samples": int(config.get("whippet_circ_min_samples", 1)),
                                 "size": int(config.get("whippet_circ_draws", 1000)),
                                 "seed": int(config.get("whippet_circ_seed", 123456)),
                                 "native": "whippet-delta defaults", "report": "bs-delta-v1"}
MODULE_FINAL[WHIPPET_CIRC] = "results.tsv.gz"


def whippet_circ_check():
    if not WHIPPET_CIRC_BIN:
        raise WorkflowError("whippet_circ runs the lab Whippet: set whippet_bin_folder (and julia)")


MODULE_EXTRA_CHECKS[WHIPPET_CIRC] = whippet_circ_check


WHIPPET_CIRC_SOURCE = config.get("whippet_circ_source", "core")


def whippet_core_identity(wildcards):
    import shlex
    identity = _json.loads(module_identity(wildcards, WHIPPET_CIRC))
    identity["flags"] = shlex.split(str(config.get("whippet_flags", "") or ""))
    identity["processing_sha256"] = UMBRELLA_MANIFEST.run_sha256(umbrella_quant_run(wildcards).run_id)
    return _json.dumps(identity, sort_keys=True)


def whippet_core_inputs(wildcards):
    record = module_run(wildcards, WHIPPET_CIRC)
    return [record.work_dir + "/whippet/quant." + suffix for suffix in
            ("psi.gz", "jnc.gz", "gene.tpm.gz", "isoform.tpm.gz", "provenance.json")]


WHIPPET_CORE_WILDCARDS = dict(MODULE_WILDCARDS)
WHIPPET_CORE_IMPORT_RUNS = []
if not MODULE_PRODUCERS and WHIPPET_CIRC_SOURCE == "core":
    import re
    WHIPPET_CORE_IMPORT_RUNS = [record.run_id for record in module_runs(WHIPPET_CIRC)
                               if not _os.path.exists(module_cache_dir(record, WHIPPET_CIRC) + "/cache.json")]
    WHIPPET_CORE_WILDCARDS["run_id"] = "(?:" + "|".join(re.escape(run) for run in WHIPPET_CORE_IMPORT_RUNS) + ")"

if WHIPPET_CIRC_SOURCE == "core" and (MODULE_PRODUCERS or WHIPPET_CORE_IMPORT_RUNS):
    rule umbrella_whippet_circ_from_core:
        input:
            whippet_core_inputs
        output:
            [protected(MODULE_CACHE.replace("{tool}", WHIPPET_CIRC) + "/" + name)
             for name in CACHE_FILES[WHIPPET_CIRC]]
        wildcard_constraints:
            **WHIPPET_CORE_WILDCARDS
        params:
            source=lambda w: umbrella_quant_run(w).work_dir + "/whippet",
            out=MODULE_CACHE.replace("{tool}", WHIPPET_CIRC),
            identity=whippet_core_identity
        shell:
            "python3 src/whippet_core.py import --source {params.source:q} --out-dir {params.out:q} "
            "--run-id {wildcards.run_id:q} --identity {params.identity:q}"


# Explicit compatibility route for an independently quantified legacy cache.
if MODULE_PRODUCERS and WHIPPET_CIRC_SOURCE == "parallel":
    rule umbrella_whippet_circ_quant:
        input:
            reads=module_ingest(WHIPPET_CIRC, lambda w: umbrella_quant_reads(w)),
            valid=module_ingest(WHIPPET_CIRC, lambda w: umbrella_quant_run(w).work_dir + "/reads.valid"),
            index=module_ingest(WHIPPET_CIRC, lambda w: _os.path.abspath(UMBRELLA_WHIPPET_MEMBERS[0]))
        output:
            [protected(MODULE_CACHE.replace("{tool}", WHIPPET_CIRC) + "/" + name)
             for name in CACHE_FILES[WHIPPET_CIRC]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", WHIPPET_CIRC),
            reads=lambda w: umbrella_read_arguments(w, "whippet"),
            julia=WHIPPET_CIRC_JULIA,
            whippet_bin=WHIPPET_CIRC_BIN,
            # --flags=VALUE: an empty whippet_flags must still be one argument
            flags=str(config.get("whippet_flags", "") or ""),
            identity=lambda w: module_identity(w, WHIPPET_CIRC)
        threads: 1
        log:
            "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.whippet_circ.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.whippet_circ.{cache_id}.tsv"
        shell:
            "python3 src/run_whippet_circ.py quant --reads {params.reads} --index {input.index:q} "
            "--out-dir {params.out:q} --run-id {wildcards.run_id:q} --julia {params.julia:q} "
            "--whippet-bin {params.whippet_bin:q} --flags={params.flags:q} --identity {params.identity:q} "
            "> {log} 2>&1"


rule umbrella_whippet_circ_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches(WHIPPET_CIRC),
        psi=module_comparison_caches(WHIPPET_CIRC, "quant.psi.gz")
    output:
        [MODULE_ANALYSIS.replace("{tool}", WHIPPET_CIRC) + "/" + name
         for name in ("results.tsv.gz", "bs_counts.tsv.gz", "native.diff.gz", "status.json")]
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", WHIPPET_CIRC),
        caches=lambda w: module_cache_arguments(w, WHIPPET_CIRC),
        julia=WHIPPET_CIRC_JULIA,
        whippet_bin=WHIPPET_CIRC_BIN,
        settings=_json.dumps(MODULE_SETTINGS[WHIPPET_CIRC], sort_keys=True),
        selection=lambda w: module_selection_row(w, WHIPPET_CIRC),
        ids=lambda w: module_ids(w, WHIPPET_CIRC, SOFTWARE[WHIPPET_CIRC])
    log:
        MODULE_ANALYSIS.replace("{tool}", WHIPPET_CIRC) + "/compare.log"
    benchmark:
        MODULE_ANALYSIS.replace("{tool}", WHIPPET_CIRC) + "/benchmark.tsv"
    conda:
        "../envs/umbrella-circrna-diff.yaml"
    shell:
        "python3 src/run_whippet_circ.py compare --preflight {input.preflight:q} {params.caches} "
        "--out-dir {params.out:q} --julia {params.julia:q} --whippet-bin {params.whippet_bin:q} "
        "--settings {params.settings:q} --selection {params.selection:q} --ids {params.ids:q} > {log} 2>&1"


rule prepare_whippet_circ:
    input: module_target(WHIPPET_CIRC, module_prepare_targets)


rule quant_whippet_circ:
    input: module_target(WHIPPET_CIRC, module_quant_targets)

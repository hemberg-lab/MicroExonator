"""Optional Whippet --circ module: a parallel Whippet quant that also reports back-splicing.

Select it with umbrella_modules: [whippet_circ]. Its scope is "selected" by default, like
the circRNA modules. Per run, `whippet-quant --circ` on the staged reads and the umbrella
Whippet index (whippet_flags as in the core quant); the core Whippet outputs and shards are
never touched. Per comparison: whippet-delta on those psi files (native.diff.gz, linear
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


# the index path is absolute: the existing index is evidence, never a reason to rebuild the reference
if MODULE_PRODUCERS:
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
            "--whippet-bin {params.whippet_bin:q} --flags {params.flags:q} --identity {params.identity:q} "
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

"""Optional vast-tools-turbo module: `align` split into independent chunk jobs.

Select it with `umbrella_modules: [vasttools_turbo]`, beside or instead of the
native `vasttools` module (they share the pinned software, VASTDB and staged
reads, and keep separate caches and results). Per run:

  split  strandedness probe as vast-tools does it, then N read chunks
  chunk  `vast-tools align` per chunk, strandedness forced (N jobs)
  merge  `vast-tools merge` into to_combine.tar + cache.json, the native layout

`compare` is the native one (src/run_vasttools.py) with the label vasttools_turbo.
Config: vasttools_turbo_chunks (default 8; changes the identity),
vasttools_turbo_threads (default 4 per chunk; does not).
"""

TURBO = "vasttools_turbo"
TURBO_CHUNKS = int(config.get("vasttools_turbo_chunks", 8))
TURBO_THREADS = int(config.get("vasttools_turbo_threads", 4))
if TURBO_CHUNKS < 1 or TURBO_THREADS < 1:
    raise WorkflowError("vasttools_turbo_chunks and vasttools_turbo_threads must be at least 1")
if TURBO in UMBRELLA_MODULES and MODULE_SETTINGS[TURBO]["ir_version"] != 2:
    raise WorkflowError("vasttools_turbo supports vasttools_ir_version 2 only")

TURBO_WORK = ("umbrella/modules/work/vasttools_turbo/{reference_id}/{project_id}/{batch_id}/"
              "{run_id}/{cache_id}")
TURBO_ID = "{reference_id}/{project_id}/{batch_id}/{run_id}"


def turbo_work(wildcards):
    return TURBO_WORK.format(reference_id=wildcards.reference_id, project_id=wildcards.project_id,
                             batch_id=wildcards.batch_id, run_id=wildcards.run_id,
                             cache_id=wildcards.cache_id)


def turbo_chunk_files(wildcards, suffix):
    return ["{}/chunk_{}.{}".format(turbo_work(wildcards), k, suffix) for k in range(TURBO_CHUNKS)]


def turbo_mate_files(wildcards, directory, chunk):
    record = umbrella_quant_run(wildcards)
    mates = 2 if record.layout == "PE" else 1
    return " ".join("{}/c{}_R{}.fastq.gz".format(directory, chunk, m + 1) for m in range(mates))


if MODULE_PRODUCERS:
    rule umbrella_vasttools_turbo_split:
        input:
            reads=module_ingest(TURBO, lambda w: umbrella_quant_reads(w)),
            valid=module_ingest(TURBO, lambda w: umbrella_quant_run(w).work_dir + "/reads.valid"),
            vastdb=module_ingest(TURBO, VASTDB)
        output:
            split=temp(TURBO_WORK + "/split.json"),
            chunks=temp(directory(TURBO_WORK + "/chunks"))
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=TURBO_WORK,
            dbdir=VASTDB_DIR,
            species=VASTTOOLS_SPECIES,
            n=TURBO_CHUNKS
        threads: 4
        log:
            "umbrella/logs/" + TURBO_ID + ".vasttools_turbo_split.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/" + TURBO_ID + ".vasttools_turbo_split.{cache_id}.tsv"
        conda:
            "../envs/umbrella-vasttools.yaml"
        shell:
            "python3 src/run_vasttools_turbo.py split --dbdir {params.dbdir} --species {params.species} "
            "--run-id {wildcards.run_id} --reads {input.reads} --chunks {params.n} "
            "--out-dir {params.out} --threads {threads} > {log} 2>&1"


    rule umbrella_vasttools_turbo_chunk:
        input:
            split=module_ingest(TURBO, lambda w: turbo_work(w) + "/split.json"),
            chunks=module_ingest(TURBO, lambda w: turbo_work(w) + "/chunks"),
            software=module_ingest(TURBO, VASTTOOLS),
            vastdb=module_ingest(TURBO, VASTDB)
        output:
            tar=temp(TURBO_WORK + "/chunk_{chunk}.tar"),
            meta=temp(TURBO_WORK + "/chunk_{chunk}.json")
        wildcard_constraints:
            chunk="[0-9]+",
            **MODULE_WILDCARDS
        params:
            software=VASTTOOLS_SOFTWARE,
            dbdir=VASTDB_DIR,
            species=VASTTOOLS_SPECIES,
            expr="--expr" if config.get("vasttools_expr", False) else "",
            reads=lambda w, input: turbo_mate_files(w, input.chunks, w.chunk)
        threads: TURBO_THREADS
        log:
            "umbrella/logs/" + TURBO_ID + ".vasttools_turbo_chunk{chunk}.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/" + TURBO_ID + ".vasttools_turbo_chunk{chunk}.{cache_id}.tsv"
        conda:
            "../envs/umbrella-vasttools.yaml"
        shell:
            "python3 src/run_vasttools_turbo.py chunk --vast-tools {params.software} --dbdir {params.dbdir} "
            "--species {params.species} --run-id {wildcards.run_id} --split-json {input.split} "
            "--index {wildcards.chunk} --reads {params.reads} --out-tar {output.tar} "
            "--out-meta {output.meta} --threads {threads} {params.expr} > {log} 2>&1"


    rule umbrella_vasttools_turbo_merge:
        input:
            split=module_ingest(TURBO, lambda w: turbo_work(w) + "/split.json"),
            tars=module_ingest(TURBO, lambda w: turbo_chunk_files(w, "tar")),
            metas=module_ingest(TURBO, lambda w: turbo_chunk_files(w, "json")),
            software=module_ingest(TURBO, VASTTOOLS),
            vastdb=module_ingest(TURBO, VASTDB)
        output:
            [protected(MODULE_CACHE.replace("{tool}", TURBO) + "/" + name) for name in CACHE_FILES[TURBO]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", TURBO),
            software=VASTTOOLS_SOFTWARE,
            dbdir=VASTDB_DIR,
            species=VASTTOOLS_SPECIES,
            expr="--expr" if config.get("vasttools_expr", False) else "",
            identity=lambda w: module_identity(w, TURBO)
        threads: 1
        log:
            "umbrella/logs/" + TURBO_ID + ".vasttools_turbo_merge.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/" + TURBO_ID + ".vasttools_turbo_merge.{cache_id}.tsv"
        conda:
            "../envs/umbrella-vasttools.yaml"
        shell:
            "python3 src/run_vasttools_turbo.py merge --vast-tools {params.software} --dbdir {params.dbdir} "
            "--species {params.species} --run-id {wildcards.run_id} --split-json {input.split} "
            "--tars {input.tars} --metas {input.metas} --out-dir {params.out} {params.expr} "
            "--identity {params.identity:q} > {log} 2>&1"


rule umbrella_vasttools_turbo_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches(TURBO),
        software=VASTTOOLS,
        vastdb=VASTDB,
        me_db=UMBRELLA_REFERENCE["me_db"]
    output:
        [MODULE_ANALYSIS.replace("{tool}", TURBO) + "/" + name
         for name in ("inclusion.tsv.gz", "tidy.tsv.gz", "diff.tsv.gz", "results.tsv.gz",
                      "microexons.tsv.gz", "status.json")]
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", TURBO),
        caches=lambda w: module_cache_arguments(w, TURBO),
        software=VASTTOOLS_SOFTWARE,
        dbdir=VASTDB_DIR,
        species=VASTTOOLS_SPECIES,
        settings=_json.dumps(MODULE_SETTINGS[TURBO], sort_keys=True),
        ids=lambda w: module_ids(w, TURBO, SOFTWARE[TURBO])
    threads: 4
    benchmark:
        MODULE_ANALYSIS.replace("{tool}", TURBO) + "/benchmark.tsv"
    conda:
        "../envs/umbrella-vasttools.yaml"
    shell:
        "python3 src/run_vasttools.py compare --preflight {input.preflight} {params.caches} "
        "--vast-tools {params.software} --dbdir {params.dbdir} --species {params.species} "
        "--me-db {input.me_db} --out-dir {params.out} --threads {threads} "
        "--settings {params.settings:q} --ids {params.ids:q} --tool-name vasttools_turbo"


rule prepare_vasttools_turbo:
    input: module_target(TURBO, module_prepare_targets)


rule quant_vasttools_turbo:
    input: module_target(TURBO, module_quant_targets)

"""Optional circRNA modules: CIRIquant (catalogue quantification), CIRCexplorer2 and DCC
(discovery on one shared STAR chimeric alignment per run).

Select them like any module (umbrella_modules: [ciriquant, circexplorer2, dcc]). Their
scope is "selected" by default: they run only on the comparisons listed for them in
umbrella_module_selection (a TSV with comparison_id, tool, role, library_evidence,
decided_by, approved_by, reason), and keep per-run caches only for those comparisons'
runs, in the same pass as every other module.

Catalogue (CIRIquant; also used to flag catalogue-supported discovery calls):
  circrna_catalogue_sources:
    - {path: Homo-sapiens.txt, format: circpedia, tissue: Central nervous system, name: CIRCpedia_v3_CNS}
    - {path: human_bed_v3.0.txt, format: circatlas, role: annotate, name: circAtlas_v3}
    - {path: hsa_hg19_circRNA.txt, format: circbase, assembly: hg19, role: annotate, name: circBase}
  (format: bed6 | circatlas | circbase | circpedia; role: quantify (default) | annotate:
   annotate sources only add their names to BSJs a quantify source has)
  circrna_liftover_chain: hg19ToHg38.over.chain.gz   (needed for sources not on circrna_assembly)
  circrna_assembly: hg38 (default), circrna_max_span (0 = no cap), circrna_max_pseudo_bases
Sources are BED6 (0-based start). Sources on another assembly are lifted boundary by
boundary; failed rows are kept beside the lifted file. The frozen catalogue lives in
umbrella/modules/reference/<reference_id>/circrna_catalogue/<catalogue_id>/.

Per run (ingest only): ciriquant reads the staged FASTQs and the shared temporary HISAT2
BAM (paired-end only); circrna_star runs STAR chimeric detection from the staged reads.
circexplorer2 per-run calls are computed from the kept STAR cache, so they can be
(re)made in cached_only. DCC runs per comparison on the kept STAR caches.
Comparison outputs are descriptive or detection-only: no p-values are claimed.
"""

import shlex
import hashlib

# Comparison-only additions: neither statistical option changes producer cache IDs.
CIRC_DIFF_ENABLED = flag(config.get("circrna_differential", False))
CIRC_HISAT_ENABLED = flag(config.get("circrna_hisat_delta", False))
CIRC_DIFF_SETTINGS = dict(seed=int(config.get("circrna_delta_seed", 123456)),
                          size=int(config.get("circrna_delta_draws", 10000)),
                          min_reads=int(config.get("circrna_delta_min_reads", 5)),
                          min_samples=int(config.get("circrna_delta_min_samples", 2)),
                          threshold=float(config.get("circrna_delta_threshold", 0.1)))
# edgeR needs enough events after filterByExpr for normalization and dispersion; fewer is a status
CIRC_ABUNDANCE_MIN_EVENTS = int(config.get("circrna_abundance_min_events", 10))
CIRC_DIFF_ID = hashlib.sha256(_json.dumps(dict(version="boundary-delta-edger-v2",
    settings=CIRC_DIFF_SETTINGS, hisat=CIRC_HISAT_ENABLED, abundance_min_events=CIRC_ABUNDANCE_MIN_EVENTS),
    sort_keys=True).encode()).hexdigest()[:16]

from src.umbrella_modules import circrna_catalogue_id, circrna_catalogue_settings

CIRC_TOOLS = ("ciriquant", "circexplorer2", "dcc")
try:
    CIRC_CATALOGUE_SETTINGS = circrna_catalogue_settings(config)
except ValueError as error:
    raise WorkflowError(str(error))
CIRC_ASSEMBLY = CIRC_CATALOGUE_SETTINGS["assembly"]
CIRC_SOURCES = [entry if isinstance(entry, dict) else {"path": entry}
                for entry in config.get("circrna_catalogue_sources", []) or []]
CIRC_CHAIN = config.get("circrna_liftover_chain", "")
CIRC_CATALOGUE = "umbrella/modules/reference/{}/circrna_catalogue/{}".format(
    UMBRELLA_REFERENCE_ID, circrna_catalogue_id(UMBRELLA_REFERENCE_ID, config))
CIRC_CATALOGUE_BED = CIRC_CATALOGUE + "/catalogue.bed"
CIRC_GENOME = "umbrella/modules/reference/{}/circrna_genome/genome.fa".format(UMBRELLA_REFERENCE_ID)
CIRC_STAR_INDEX = MODULE_REFERENCE["circrna_star"] + "/star_index"
CIRC_ANNOTATION = MODULE_REFERENCE["circrna_star"] + "/annotation.gtf"
CIRC_REFFLAT = MODULE_REFERENCE["circexplorer2"] + "/refflat.txt"


def circ_catalogue_check():
    if not CIRC_SOURCES:
        raise WorkflowError("ciriquant quantifies a frozen circRNA catalogue: set circrna_catalogue_sources "
                            "(BED6 files; e.g. circAtlas hg38, circBase lifted from hg19)")
    missing = [entry["path"] for entry in CIRC_SOURCES if not _os.path.isfile(entry["path"])]
    if missing:
        raise WorkflowError("circrna_catalogue_sources not found: {}".format(", ".join(missing)))
    if CIRC_CATALOGUE_SETTINGS["liftover_chain"] == "missing" or (
            CIRC_CATALOGUE_SETTINGS["liftover_chain"].startswith("missing:")):
        raise WorkflowError("circrna_catalogue_sources has sources not on {} (e.g. circBase hg19): set "
                            "circrna_liftover_chain to the UCSC chain, e.g. hg19ToHg38.over.chain.gz".format(
                                CIRC_ASSEMBLY))


def ciriquant_check():
    circ_catalogue_check()
    single = [comparison_id for comparison_id, result in module_comparisons("ciriquant")
              if set(result.get("layouts", [])) != {"PE"}]
    if single:
        raise WorkflowError("ciriquant needs paired-end runs (its realignment uses -1/-2); single-end or mixed "
                            "comparisons: {}. Select circexplorer2 and dcc for them instead".format(", ".join(single)))


MODULE_EXTRA_CHECKS["ciriquant"] = ciriquant_check


def circ_catalogue_flag(wildcards):
    """--catalogue for the discovery support tables, when a catalogue is configured."""
    return "--catalogue " + CIRC_CATALOGUE_BED if CIRC_SOURCES else ""


def circ_catalogue_inputs(wildcards):
    return [CIRC_CATALOGUE_BED] if CIRC_SOURCES else []


# ------------------------------------------------------------ references

localrules: umbrella_circrna_genome


rule umbrella_circrna_genome:
    input:
        UMBRELLA_REFERENCE["genome_fasta"]
    output:
        # not protected(): when this is a link, protecting it would chmod the shared genome itself
        genome=CIRC_GENOME,
        fai=protected(CIRC_GENOME + ".fai")
    conda:
        "../envs/umbrella-ciriquant.yaml"
    shell:
        # CIRIquant and CIRCexplorer2 read the FASTA by byte offset: a plain copy or link with its .fai
        "if gzip -t {input} 2>/dev/null; then gzip -dc {input} > {output.genome}; "
        "else ln -sf $(cd $(dirname {input}) && pwd)/$(basename {input}) {output.genome}; fi "
        "&& samtools faidx {output.genome} --fai-idx {output.fai}"


def circ_source(wildcards):
    entry = CIRC_SOURCES[int(wildcards.index)]
    return {"bed": entry["path"], "chain": [CIRC_CHAIN] if entry.get("assembly", CIRC_ASSEMBLY) != CIRC_ASSEMBLY else []}


def circ_source_params(wildcards):
    entry = CIRC_SOURCES[int(wildcards.index)]
    tissue = entry.get("tissue", "")
    return "--format {} {}".format(entry.get("format", "bed6"), "--tissue {}".format(shlex.quote(tissue)) if tissue else "")


CIRC_QUANTIFY = [index for index, entry in enumerate(CIRC_SOURCES) if entry.get("role", "quantify") == "quantify"]
CIRC_ANNOTATE = [index for index, entry in enumerate(CIRC_SOURCES) if entry.get("role", "quantify") == "annotate"]
CIRC_INDEX = MODULE_REFERENCE["ciriquant"] + "/pseudo_index"


rule umbrella_circrna_source:
    input:
        unpack(circ_source)
    output:
        CIRC_CATALOGUE + "/sources/{index}.bed"
    wildcard_constraints:
        index="[0-9]+"
    params:
        prefix=CIRC_CATALOGUE + "/sources/{index}",
        convert=circ_source_params
    shell:
        # each database's layout -> BED6; sources on another assembly are then lifted boundary by boundary
        "python3 src/circrna_catalogue.py convert --input {input.bed} {params.convert} --out {params.prefix}.bed6 "
        "--report {params.prefix}.convert.json && "
        "if [ -n '{input.chain}' ]; then python3 src/circrna_catalogue.py lift --input {params.prefix}.bed6 "
        "--chain {input.chain} --out {output} --failed {params.prefix}.failed.tsv --report {params.prefix}.lift.json "
        "&& rm {params.prefix}.bed6; else mv {params.prefix}.bed6 {output}; fi"


rule umbrella_circrna_catalogue:
    input:
        sources=[CIRC_CATALOGUE + "/sources/{}.bed".format(index) for index in CIRC_QUANTIFY],
        annotate=[CIRC_CATALOGUE + "/sources/{}.bed".format(index) for index in CIRC_ANNOTATE],
        genome=CIRC_GENOME,
        fai=CIRC_GENOME + ".fai"
    output:
        bed=protected(CIRC_CATALOGUE_BED),
        annotation=protected(CIRC_CATALOGUE + "/annotation.tsv.gz"),
        report=CIRC_CATALOGUE + "/report.json"
    params:
        out=CIRC_CATALOGUE,
        max_span=CIRC_CATALOGUE_SETTINGS["max_span"],
        max_pseudo=CIRC_CATALOGUE_SETTINGS["max_pseudo_bases"]
    shell:
        "python3 src/circrna_catalogue.py build --input {input.sources} --genome {input.genome} "
        "--out-dir {params.out} --max-span {params.max_span} --max-pseudo-bases {params.max_pseudo} "
        "--annotate {input.annotate}"


rule umbrella_ciriquant_index:
    input:
        catalogue=CIRC_CATALOGUE_BED,
        genome=CIRC_GENOME,
        fai=CIRC_GENOME + ".fai"
    output:
        protected(directory(CIRC_INDEX))
    threads: 8
    log:
        MODULE_REFERENCE["ciriquant"] + "/pseudo_index.log"
    conda:
        "../envs/umbrella-ciriquant.yaml"
    shell:
        "python3 src/run_ciriquant.py index --catalogue {input.catalogue} --genome {input.genome} "
        "--out-dir {output} --threads {threads} > {log} 2>&1"


rule umbrella_circrna_annotation:
    input:
        UMBRELLA_GTF["annotation"]
    output:
        protected(CIRC_ANNOTATION)
    shell:
        "gzip -cdf {input} > {output}"


rule umbrella_circrna_star_index:
    input:
        genome=CIRC_GENOME,
        gtf=CIRC_ANNOTATION
    output:
        protected(directory(CIRC_STAR_INDEX))
    params:
        overhang=int(config.get("circrna_star_sjdb_overhang", 100)),
        sa_bases=int(config.get("circrna_star_sa_index_nbases", 14))
    threads: 8
    log:
        MODULE_REFERENCE["circrna_star"] + "/star_index.log"
    conda:
        "../envs/umbrella-circrna-discovery.yaml"
    shell:
        # run_circrna_discovery.py picks a STAR build the node can run (see STAR_BUILDS there)
        "python3 src/run_circrna_discovery.py star-index --genome {input.genome} --gtf {input.gtf} "
        "--out-dir {output} --threads {threads} --overhang {params.overhang} --sa-bases {params.sa_bases} "
        "> {log} 2>&1"


rule umbrella_circrna_refflat:
    input:
        CIRC_ANNOTATION
    output:
        protected(CIRC_REFFLAT)
    shell:
        "python3 src/run_circrna_discovery.py refflat --gtf {input} --out {output}"


# ------------------------------------------------------------ per run

if MODULE_PRODUCERS:
    rule umbrella_ciriquant_prep:
        input:
            reads=module_ingest("ciriquant", lambda w: umbrella_quant_reads(w)),
            valid=module_ingest("ciriquant", lambda w: umbrella_quant_run(w).work_dir + "/reads.valid"),
            bam=module_ingest("ciriquant", lambda w: umbrella_quant_run(w).work_dir + "/align/aligned.bam"),
            bai=module_ingest("ciriquant", lambda w: umbrella_quant_run(w).work_dir + "/align/aligned.bam.bai"),
            summary=module_ingest("ciriquant", lambda w: umbrella_quant_run(w).work_dir +
                                  "/align/junctions.summary.json"),
            catalogue=module_ingest("ciriquant", CIRC_CATALOGUE_BED),
            index=module_ingest("ciriquant", CIRC_INDEX),
            genome=module_ingest("ciriquant", CIRC_GENOME),
            fai=module_ingest("ciriquant", CIRC_GENOME + ".fai")
        output:
            [protected(MODULE_CACHE.replace("{tool}", "ciriquant") + "/" + name) for name in CACHE_FILES["ciriquant"]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", "ciriquant"),
            strandedness=lambda w: umbrella_quant_run(w).strandedness,
            layout=lambda w: umbrella_quant_run(w).layout,
            anchor=int(config.get("ciriquant_anchor", 5)),
            identity=lambda w: module_identity(w, "ciriquant")
        threads: 8
        log:
            "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.ciriquant.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.ciriquant.{cache_id}.tsv"
        conda:
            "../envs/umbrella-ciriquant.yaml"
        shell:
            "python3 src/run_ciriquant.py prep --reads {input.reads} --bam {input.bam} --catalogue {input.catalogue} "
            "--genome {input.genome} --out-dir {params.out} --run-id {wildcards.run_id} --threads {threads} "
            "--anchor {params.anchor} --strandedness {params.strandedness} --summary {input.summary} "
            "--layout {params.layout} --index-dir {input.index} --identity {params.identity:q} > {log} 2>&1"

    rule umbrella_circrna_star:
        input:
            reads=module_ingest("circrna_star", lambda w: umbrella_quant_reads(w)),
            valid=module_ingest("circrna_star", lambda w: umbrella_quant_run(w).work_dir + "/reads.valid"),
            index=module_ingest("circrna_star", CIRC_STAR_INDEX)
        output:
            [protected(MODULE_CACHE.replace("{tool}", "circrna_star") + "/" + name)
             for name in CACHE_FILES["circrna_star"]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", "circrna_star"),
            identity=lambda w: module_identity(w, "circrna_star")
        threads: 8
        log:
            "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.circrna_star.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.circrna_star.{cache_id}.tsv"
        conda:
            "../envs/umbrella-circrna-discovery.yaml"
        shell:
            "python3 src/run_circrna_discovery.py star --reads {input.reads} --index {input.index} "
            "--out-dir {params.out} --run-id {wildcards.run_id} --threads {threads} "
            "--identity {params.identity:q} > {log} 2>&1"


def circ_star_cache(wildcards, name="cache.json"):
    return module_cache_dir(umbrella_quant_run(wildcards), "circrna_star") + "/" + name


# computed from the kept STAR cache only, so it exists in both modes: no reads, no alignment
rule umbrella_circexplorer2_calls:
    input:
        star=lambda w: [circ_star_cache(w, name) for name in CACHE_FILES["circrna_star"]],
        refflat=CIRC_REFFLAT,
        genome=CIRC_GENOME,
        fai=CIRC_GENOME + ".fai"
    output:
        [protected(MODULE_CACHE.replace("{tool}", "circexplorer2") + "/" + name)
         for name in CACHE_FILES["circexplorer2"]]
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        star=lambda w: module_cache_dir(module_run(w, "circexplorer2"), "circrna_star"),
        out=MODULE_CACHE.replace("{tool}", "circexplorer2"),
        identity=lambda w: module_identity(w, "circexplorer2")
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.circexplorer2.{cache_id}.log"
    conda:
        "../envs/umbrella-circrna-discovery.yaml"
    shell:
        "python3 src/run_circrna_discovery.py circexplorer2 --star-cache {params.star} --refflat {input.refflat} "
        "--genome {input.genome} --out-dir {params.out} --run-id {wildcards.run_id} "
        "--identity {params.identity:q} > {log} 2>&1"


# ------------------------------------------------------------ per comparison

rule umbrella_ciriquant_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches("ciriquant")
    output:
        [MODULE_ANALYSIS.replace("{tool}", "ciriquant") + "/" + name
         for name in ("results.tsv.gz", "bsj_counts.tsv.gz", "fsj_counts.tsv.gz", "library.tsv", "status.json")]
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "ciriquant"),
        caches=lambda w: module_cache_arguments(w, "ciriquant"),
        settings=_json.dumps(MODULE_SETTINGS["ciriquant"], sort_keys=True),
        selection=lambda w: module_selection_row(w, "ciriquant"),
        ids=lambda w: module_ids(w, "ciriquant", SOFTWARE["ciriquant"])
    conda:
        "../envs/umbrella-ciriquant.yaml"
    shell:
        "python3 src/run_ciriquant.py compare --preflight {input.preflight} {params.caches} --out-dir {params.out} "
        "--settings {params.settings:q} --selection {params.selection:q} --ids {params.ids:q}"


rule umbrella_circexplorer2_support:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches("circexplorer2"),
        catalogue=circ_catalogue_inputs
    output:
        [MODULE_ANALYSIS.replace("{tool}", "circexplorer2") + "/" + name
         for name in ("results.tsv.gz", "support.tsv.gz", "status.json")]
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "circexplorer2"),
        caches=lambda w: module_cache_arguments(w, "circexplorer2"),
        catalogue=circ_catalogue_flag,
        selection=lambda w: module_selection_row(w, "circexplorer2"),
        ids=lambda w: module_ids(w, "circexplorer2", SOFTWARE["circexplorer2"])
    conda:
        "../envs/umbrella-circrna-discovery.yaml"
    shell:
        "python3 src/run_circrna_discovery.py support --preflight {input.preflight} {params.caches} "
        "--out-dir {params.out} {params.catalogue} --selection {params.selection:q} --ids {params.ids:q}"


def dcc_star_caches(wildcards):
    result = module_comparison(wildcards, "dcc")
    return [module_cache_dir(UMBRELLA_MANIFEST.by_run[unit], "circrna_star") + "/" + name
            for unit in sorted(result["collapse"]) for name in CACHE_FILES["circrna_star"]]


def dcc_cache_arguments(wildcards):
    result = module_comparison(wildcards, "dcc")
    return " ".join("--cache {}={}".format(unit, module_cache_dir(UMBRELLA_MANIFEST.by_run[unit], "circrna_star"))
                    for unit in sorted(result["collapse"]))


rule umbrella_dcc_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=dcc_star_caches,
        gtf=CIRC_ANNOTATION,
        catalogue=circ_catalogue_inputs
    output:
        [MODULE_ANALYSIS.replace("{tool}", "dcc") + "/" + name
         for name in ("results.tsv.gz", "support.tsv.gz", "CircRNACount.tsv.gz", "CircCoordinates.tsv.gz",
                      "status.json")]
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "dcc"),
        caches=dcc_cache_arguments,
        catalogue=circ_catalogue_flag,
        settings=_json.dumps(MODULE_SETTINGS["dcc"], sort_keys=True),
        selection=lambda w: module_selection_row(w, "dcc"),
        ids=lambda w: module_ids(w, "dcc", SOFTWARE["dcc"])
    threads: 4
    conda:
        "../envs/umbrella-circrna-discovery.yaml"
    shell:
        "python3 src/run_circrna_discovery.py dcc --preflight {input.preflight} {params.caches} --gtf {input.gtf} "
        "--out-dir {params.out} {params.catalogue} --settings {params.settings:q} --threads {threads} "
        "--selection {params.selection:q} --ids {params.ids:q}"


rule prepare_ciriquant:
    input: module_target("ciriquant", module_prepare_targets)

rule quant_ciriquant:
    input: module_target("ciriquant", module_quant_targets)

rule prepare_circexplorer2:
    input: module_target("circexplorer2", module_prepare_targets)

rule quant_circexplorer2:
    input: module_target("circexplorer2", module_quant_targets)

rule prepare_dcc:
    input: module_target("dcc", module_prepare_targets)

rule quant_dcc:
    input: module_target("dcc", module_quant_targets)


def circ_diff_source(wildcards):
    module_comparison(wildcards, wildcards.tool)
    if wildcards.diff_id != CIRC_DIFF_ID:
        raise WorkflowError("circRNA differential settings ID mismatch")
    root = MODULE_ANALYSIS.format(**wildcards)
    names = (["bsj_counts.tsv.gz", "results.tsv.gz"] if wildcards.tool == "ciriquant"
             else ["support.tsv.gz"])
    return [root + "/" + name for name in names]


def circ_diff_shards(wildcards):
    if not CIRC_HISAT_ENABLED:
        return []
    result = module_comparison(wildcards, wildcards.tool)
    records = [UMBRELLA_MANIFEST.by_run[unit] for unit in result["collapse"]]
    paths = sorted(set("junctions/{}/{}/{}/{}.junctions.tsv.gz".format(
        r.reference_id, r.project_id, r.group, r.batch_id) for r in records))
    missing = [path for path in paths if not _os.path.isfile(path)]
    if missing:
        raise WorkflowError("Bayesian circRNA usage requires retained HISAT2 junction shards; "
                            "no read restaging is permitted for this step. Missing: " + ", ".join(missing))
    # Treat existing immutable shards as external evidence. Relative paths would
    # match umbrella_junction_shard and traverse back into missing BAM/FASTQs.
    return [_os.path.abspath(path) for path in paths]


rule umbrella_circrna_differential:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        tables=circ_diff_source,
        shards=circ_diff_shards
    output:
        status=MODULE_ANALYSIS + "/differential/{diff_id}/status.tsv",
        abundance=MODULE_ANALYSIS + "/differential/{diff_id}/abundance.tsv.gz",
        bayesian=MODULE_ANALYSIS + "/differential/{diff_id}/hisat_delta.tsv.gz",
        evidence=MODULE_ANALYSIS + "/differential/{diff_id}/hisat_evidence.tsv.gz",
        metadata=MODULE_ANALYSIS + "/differential/{diff_id}/preparation.json"
    wildcard_constraints:
        tool="ciriquant|circexplorer2|dcc",
        diff_id="[0-9a-f]{16}"
    params:
        source=MODULE_ANALYSIS,
        out=MODULE_ANALYSIS + "/differential/{diff_id}",
        settings=_json.dumps(CIRC_DIFF_SETTINGS, sort_keys=True),
        shards=lambda w, input: " ".join("--shard " + shlex.quote(p) for p in input.shards),
        min_events=CIRC_ABUNDANCE_MIN_EVENTS
    threads: 1
    conda:
        "../envs/umbrella-circrna-diff.yaml"
    shell:
        "python3 src/circrna_hisat_delta.py --tool {wildcards.tool} --preflight {input.preflight:q} "
        "--source {params.source:q} --out {params.out:q} --settings {params.settings:q} {params.shards} "
        "&& Rscript src/circrna_abundance.R {params.out:q} {params.min_events}"

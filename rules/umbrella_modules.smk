"""Optional umbrella modules: MAJIQ v3, DaPars2 and QAPA (src/umbrella_modules.py).

Selection (config):
  umbrella_modules: [majiq, dapars2, qapa]   any subset; default none
  umbrella_modules_mode: cached_only | ingest (default cached_only)

Targets:
  prepare_majiq / prepare_dapars2 / prepare_qapa   per-run caches of that tool
  quant_majiq / quant_dapars2 / quant_qapa         caches plus every comparison
  quant_umbrella_modules                           the selected modules, plus an
                                                   inventory of their outputs
  quant_umbrella                                   also includes the selected ones

A named target implies its tool even when it is not selected, and pulls no
other tool: MAJIQ and DaPars2 read the shared temporary HISAT2 BAM; QAPA reads
the staged reads with its own Salmon 3' UTR index. Comparison rules read only
the per-run caches and the comparison preflight, never reads or BAMs.
HISAT2's inputs are deliberately unchanged (Snakemake would re-align every
existing run if its input set changed), so on a brand-new reference the
alignment still waits for the full reference manifest; on an existing one
nothing extra is built.

Layout (IDs are computed, see src/umbrella_modules.py):
  umbrella/modules/reference/<reference_id>/<tool>/<module_reference_id>/
  umbrella/modules/cache/<reference_id>/<project>/<batch>/<run>/<tool>/<cache_id>/
  comparisons/<reference_id>/<project>/<comparison>/modules/<tool>/<analysis_id>/
Per-run caches are protected; comparison outputs are recomputed freely.
"""

import json as _json
import os as _os

from src.umbrella_modules import (CACHE_FILES, DAPARS2_COMMIT, MODULES, SOFTWARE, analysis_id,
                                  cache_id, missing_caches, module_reference_id, parse_selection)

try:
    UMBRELLA_MODULES, UMBRELLA_MODULES_MODE = parse_selection(config)
except ValueError as error:
    raise WorkflowError(str(error))

MODULE_REFERENCE_ID = {tool: module_reference_id(UMBRELLA_REFERENCE_ID, tool, config) for tool in MODULES}
# DaPars2 and QAPA locate 3' UTRs from coding ends, so they may need a different
# GTF than the splicing tools (e.g. GENCODE comprehensive next to an exon-only
# Whippet build). It changes only their module_reference_id.
APA_GTF = config.get("apa_annotation_gtf") or UMBRELLA_REFERENCE["annotation_gtf"]
APA_TRANSCRIPTS = str(config.get("apa_transcripts", "basic"))
APA_MIN_UTRS = int(config.get("apa_min_utrs", 5000))
if APA_TRANSCRIPTS not in ("basic", "all"):
    raise WorkflowError("apa_transcripts must be basic or all, not {}".format(APA_TRANSCRIPTS))
if bool(config.get("qapa_gencode_polya")) != bool(config.get("qapa_polyasite")):
    raise WorkflowError("qapa_gencode_polya and qapa_polyasite go together (QAPA's -g and -p); "
                        "use qapa_polya_sites alone for a custom BED")
MODULE_REFERENCE = {tool: "umbrella/modules/reference/{}/{}/{}".format(
    UMBRELLA_REFERENCE_ID, tool, MODULE_REFERENCE_ID[tool]) for tool in MODULES}
MODULE_CACHE = ("umbrella/modules/cache/{reference_id}/{project_id}/{batch_id}/{run_id}/"
                "{tool}/{cache_id}")
MODULE_ANALYSIS = UMBRELLA_COMPARISON_ROOT + "/modules/{tool}/{analysis_id}"
MODULE_SETTINGS = {
    "majiq": dict({key: config.get("majiq_" + key + "_args") for key in ("update", "psicov", "heterogen")},
                  statistic=str(config.get("majiq_het_statistic", "ttest"))),
    "dapars2": {"coverage_threshold": int(config.get("dapars2_coverage_threshold", 10))},
    "qapa": {"test": "DEXSeq ~ sample + exon + condition:exon, genes with >= 2 sites"},
}
# Per-run cache producers are defined only in ingest mode. In cached_only the
# kept caches are plain inputs with no rule that could remake them, so neither
# Snakemake's records (an earlier ingest's input set) nor a rebuilt reference
# can ever schedule a download, alignment or re-quantification.
MODULE_PRODUCERS = UMBRELLA_MODULES_MODE == "ingest"
MODULE_WILDCARDS = dict(reference_id="[^/]+", project_id="[^/]+", batch_id="[^/]+",
                        run_id="[^/]+", cache_id="[0-9a-f]{16}", analysis_id="[0-9a-f]{16}")


def module_cache_dir(record, tool):
    return MODULE_CACHE.format(reference_id=record.reference_id, project_id=record.project_id,
                               batch_id=record.batch_id, run_id=record.run_id, tool=tool,
                               cache_id=cache_id(record, MODULE_REFERENCE_ID[tool]))


def module_cache_outputs(record, tool):
    return [module_cache_dir(record, tool) + "/" + name for name in CACHE_FILES[tool]]


def module_run(wildcards, tool):
    """The manifest run of a cache path; refuses stale IDs and, in cached_only, making caches."""
    record = umbrella_quant_run(wildcards)
    expected = cache_id(record, MODULE_REFERENCE_ID[tool])
    if wildcards.cache_id != expected:
        raise WorkflowError("{} cache {} of run {} is not the current one ({})".format(
            tool, wildcards.cache_id, record.run_id, expected))
    return record


def module_ingest(tool, inputs):
    """A per-run producer's inputs, after checking its cache ID is the current one.

    Producers exist only in ingest mode (see MODULE_PRODUCERS), so this never
    changes a producer's input set between modes: Snakemake would otherwise
    rerun it with "Set of input files has changed"."""
    def function(wildcards):
        module_run(wildcards, tool)
        return inputs(wildcards) if callable(inputs) else inputs
    return function


def module_identity(wildcards, tool):
    record = umbrella_quant_run(wildcards)
    return _json.dumps({"cache_id": wildcards.cache_id, "reference_id": record.reference_id,
                        "module_reference_id": MODULE_REFERENCE_ID[tool], "software": SOFTWARE[tool]},
                       sort_keys=True)


def module_analysis_id(result, tool):
    return analysis_id(result, MODULE_REFERENCE_ID[tool], MODULE_SETTINGS[tool])


def module_analysis_dir(result, tool):
    return MODULE_ANALYSIS.format(reference_id=result["reference_id"], project_id=result["project_id"],
                                  comparison_id=result["comparison_id"], tool=tool,
                                  analysis_id=module_analysis_id(result, tool))


def module_comparison(wildcards, tool):
    result = umbrella_comparison(wildcards)
    if wildcards.analysis_id != module_analysis_id(result, tool):
        raise WorkflowError("{} analysis {} of {} is not the current one".format(
            tool, wildcards.analysis_id, result["comparison_id"]))
    return result


def module_comparison_caches(tool, name="cache.json"):
    def paths(wildcards):
        result = module_comparison(wildcards, tool)
        return [module_cache_dir(UMBRELLA_MANIFEST.by_run[run], tool) + "/" + name
                for run in sorted(result["collapse"])]
    return paths


def module_cache_arguments(wildcards, tool):
    result = module_comparison(wildcards, tool)
    return " ".join("--cache {}={}".format(run, module_cache_dir(UMBRELLA_MANIFEST.by_run[run], tool))
                    for run in sorted(result["collapse"]))


def module_ids(wildcards, tool, version):
    result = module_comparison(wildcards, tool)
    return _json.dumps({"reference_id": result["reference_id"], "module_reference_id":
                        MODULE_REFERENCE_ID[tool], "analysis_id": wildcards.analysis_id,
                        "tool_version": version}, sort_keys=True)


# ------------------------------------------------------------ DaPars2

DAPARS2_SOFTWARE = "umbrella/modules/software/dapars2/" + DAPARS2_COMMIT
DAPARS2_UTR = MODULE_REFERENCE["dapars2"] + "/utr.bed"
DAPARS2_WINDOWS = MODULE_REFERENCE["dapars2"] + "/utr_windows.bed"

localrules: umbrella_dapars2_software


rule umbrella_dapars2_software:
    output:
        DAPARS2_SOFTWARE + "/src/DaPars2_Multi_Sample_Multi_Chr.py"
    params:
        root=DAPARS2_SOFTWARE,
        url="https://github.com/3UTR/DaPars2/archive/{}.tar.gz".format(DAPARS2_COMMIT)
    shell:
        # pinned upstream commit (GPL-2.0), unpacked once
        "mkdir -p {params.root} && curl -fsSL {params.url} | tar xz --strip-components 1 -C {params.root}"


rule umbrella_dapars2_reference:
    input:
        APA_GTF
    output:
        bed=protected(DAPARS2_UTR),
        windows=protected(DAPARS2_WINDOWS),
        report=MODULE_REFERENCE["dapars2"] + "/utr_report.json"
    params:
        transcripts=APA_TRANSCRIPTS,
        min_utrs=APA_MIN_UTRS
    shell:
        "python3 src/apa_annotation.py dapars-utr --gtf {input} --bed {output.bed} "
        "--windows {output.windows} --report {output.report} --transcripts {params.transcripts} "
        "--min-utrs {params.min_utrs}"


if MODULE_PRODUCERS:
    rule umbrella_dapars2_coverage:
        input:
            bam=module_ingest("dapars2", lambda w: UMBRELLA_WORK.format(**w) + "/align/aligned.bam"),
            bai=module_ingest("dapars2", lambda w: UMBRELLA_WORK.format(**w) + "/align/aligned.bam.bai"),
            windows=module_ingest("dapars2", DAPARS2_WINDOWS)
        output:
            [protected(MODULE_CACHE.replace("{tool}", "dapars2") + "/" + name)
             for name in CACHE_FILES["dapars2"]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", "dapars2"),
            identity=lambda w: module_identity(w, "dapars2")
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.dapars2.{cache_id}.tsv"
        conda:
            "../envs/umbrella-dapars2.yaml"
        shell:
            "python3 src/run_dapars2.py coverage --bam {input.bam} --windows {input.windows} "
            "--out-dir {params.out} --run-id {wildcards.run_id} --identity {params.identity:q}"


rule umbrella_dapars2_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches("dapars2"),
        utr=DAPARS2_UTR,
        software=rules.umbrella_dapars2_software.output
    output:
        native=MODULE_ANALYSIS.replace("{tool}", "dapars2") + "/native.tsv.gz",
        results=MODULE_ANALYSIS.replace("{tool}", "dapars2") + "/results.tsv.gz",
        status=MODULE_ANALYSIS.replace("{tool}", "dapars2") + "/status.json"
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "dapars2"),
        caches=lambda w: module_cache_arguments(w, "dapars2"),
        software=DAPARS2_SOFTWARE + "/src",
        threshold=MODULE_SETTINGS["dapars2"]["coverage_threshold"],
        ids=lambda w: module_ids(w, "dapars2", SOFTWARE["dapars2"])
    threads: 4
    benchmark:
        MODULE_ANALYSIS.replace("{tool}", "dapars2") + "/benchmark.tsv"
    conda:
        "../envs/umbrella-dapars2.yaml"
    shell:
        "python3 src/run_dapars2.py compare --preflight {input.preflight} {params.caches} "
        "--utr-bed {input.utr} --dapars2-dir {params.software} --out-dir {params.out} "
        "--threads {threads} --coverage-threshold {params.threshold} --ids {params.ids:q}"


# ------------------------------------------------------------ QAPA

QAPA_DB = MODULE_REFERENCE["qapa"] + "/ensembl_identifiers.txt"
QAPA_UTRS = MODULE_REFERENCE["qapa"] + "/qapa_3utrs.bed"
QAPA_FASTA = MODULE_REFERENCE["qapa"] + "/qapa_3utrs.fa"
QAPA_INDEX = MODULE_REFERENCE["qapa"] + "/salmon_index"
QAPA_DECOYS = bool(config.get("qapa_decoys", False))
QAPA_GENCODE_POLYA = MODULE_REFERENCE["qapa"] + "/gencode_polya_sites.bed"


def qapa_site_inputs(wildcards):
    if config.get("qapa_polya_sites"):
        return {"custom": config["qapa_polya_sites"]}
    if config.get("qapa_gencode_polya"):
        return {"gencode": QAPA_GENCODE_POLYA, "polyasite": config["qapa_polyasite"]}
    return {}


def qapa_site_arguments(wildcards, input):
    if hasattr(input, "custom"):
        return "-o {}".format(input.custom)
    if hasattr(input, "gencode"):
        return "-g {} -p {}".format(input.gencode, input.polyasite)
    return "-N"        # annotation-only: no poly(A) database


rule umbrella_qapa_polya:
    input:
        lambda w: config["qapa_gencode_polya"]
    output:
        protected(QAPA_GENCODE_POLYA)
    shell:
        "python3 src/apa_annotation.py polya-bed --gtf {input} --output {output}"


rule umbrella_qapa_db:
    input:
        APA_GTF
    output:
        protected(QAPA_DB)
    shell:
        "python3 src/apa_annotation.py qapa-db --gtf {input} --output {output}"


rule umbrella_qapa_build:
    input:
        unpack(lambda w: dict(qapa_site_inputs(w), gtf=APA_GTF, db=QAPA_DB))
    output:
        protected(QAPA_UTRS)
    params:
        selected=MODULE_REFERENCE["qapa"] + "/selected_transcripts.gtf",
        genepred=MODULE_REFERENCE["qapa"] + "/genes.genePred",
        transcripts=APA_TRANSCRIPTS,
        min_utrs=APA_MIN_UTRS,
        sites=qapa_site_arguments
    log:
        MODULE_REFERENCE["qapa"] + "/qapa_build.log"
    conda:
        "../envs/umbrella-qapa.yaml"
    shell:
        "python3 src/apa_annotation.py qapa-gtf --gtf {input.gtf} --output {params.selected} "
        "--transcripts {params.transcripts} "
        "&& gtfToGenePred -genePredExt {params.selected} {params.genepred} "
        "&& qapa build {params.sites} --db {input.db} {params.genepred} > {output} 2> {log} "
        "&& rm -f {params.selected} {params.genepred} "
        "&& python3 src/apa_annotation.py check-bed --bed {output} --min {params.min_utrs} "
        "--label 'QAPA 3 prime UTR library'"


rule umbrella_qapa_fasta:
    input:
        genome=UMBRELLA_REFERENCE["genome_fasta"],
        utrs=QAPA_UTRS
    output:
        protected(QAPA_FASTA)
    params:
        decoys="--decoys -d {}/decoys.txt".format(MODULE_REFERENCE["qapa"]) if QAPA_DECOYS else ""
    conda:
        "../envs/umbrella-qapa.yaml"
    shell:
        "qapa fasta -f {input.genome} {params.decoys} {input.utrs} {output}"


rule umbrella_qapa_index:
    input:
        QAPA_FASTA
    output:
        protected(directory(QAPA_INDEX))
    params:
        decoys="-d {}/decoys.txt".format(MODULE_REFERENCE["qapa"]) if QAPA_DECOYS else ""
    threads: 8
    log:
        MODULE_REFERENCE["qapa"] + "/salmon_index.log"
    conda:
        "../envs/umbrella-quant.yaml"
    shell:
        "salmon index -t {input} -i {output} {params.decoys} -p {threads} > {log} 2>&1"


if MODULE_PRODUCERS:
    rule umbrella_qapa_quant:
        input:
            reads=module_ingest("qapa", lambda w: umbrella_quant_reads(w)),
            valid=module_ingest("qapa", lambda w: umbrella_quant_run(w).work_dir + "/reads.valid"),
            index=module_ingest("qapa", QAPA_INDEX)
        output:
            [protected(MODULE_CACHE.replace("{tool}", "qapa") + "/" + name) for name in CACHE_FILES["qapa"]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", "qapa"),
            reads=lambda w: umbrella_read_arguments(w, "salmon"),
            library=config.get("qapa_library_type", "A"),
            identity=lambda w: module_identity(w, "qapa")
        threads: 4
        log:
            "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.qapa.{cache_id}.salmon.log"
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.qapa.{cache_id}.tsv"
        conda:
            "../envs/umbrella-quant.yaml"
        shell:
            "rm -rf {params.out}/salmon && salmon quant -i {input.index} -l {params.library} {params.reads} "
            "--validateMappings -p {threads} -o {params.out}/salmon 2> {log} "
            "&& python3 src/run_qapa.py keep --salmon {params.out}/salmon --out-dir {params.out} "
            "--run-id {wildcards.run_id} --identity {params.identity:q} && rm -rf {params.out}/salmon"


rule umbrella_qapa_pau:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches("qapa"),
        db=QAPA_DB
    output:
        pau=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/pau.tsv",
        counts=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/site_counts.tsv.gz",
        samples=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/samples.tsv"
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "qapa"),
        caches=lambda w: module_cache_arguments(w, "qapa")
    conda:
        "../envs/umbrella-qapa.yaml"
    shell:
        "python3 src/run_qapa.py pau --preflight {input.preflight} {params.caches} "
        "--db {input.db} --out-dir {params.out}"


rule umbrella_qapa_usage:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        pau=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/pau.tsv",
        counts=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/site_counts.tsv.gz",
        samples=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/samples.tsv"
    output:
        dexseq=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/dexseq.tsv.gz",
        results=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/results.tsv.gz",
        status=MODULE_ANALYSIS.replace("{tool}", "qapa") + "/status.json"
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "qapa"),
        ids=lambda w: module_ids(w, "qapa", SOFTWARE["qapa"])
    threads: 2
    log:
        MODULE_ANALYSIS.replace("{tool}", "qapa") + "/dexseq.log"
    benchmark:
        MODULE_ANALYSIS.replace("{tool}", "qapa") + "/benchmark.tsv"
    conda:
        "../envs/umbrella-apa-stats.yaml"
    shell:
        "Rscript src/qapa_dexseq.R {input.counts} {input.samples} {output.dexseq} {threads} > {log} 2>&1 "
        "&& python3 src/run_qapa.py normalize --preflight {input.preflight} --pau {input.pau} "
        "--dexseq {output.dexseq} --out-dir {params.out} --ids {params.ids:q}"


# ------------------------------------------------------------ MAJIQ v3

MAJIQ_BIN_FOLDER = config.get("majiq_bin_folder", "")
MAJIQ_SOFTWARE_ROOT = "umbrella/modules/software/majiq/" + MODULE_REFERENCE_ID["majiq"]
MAJIQ_INSTALLED = [] if MAJIQ_BIN_FOLDER else [MAJIQ_SOFTWARE_ROOT + "/installed.json"]
MAJIQ_BIN = MAJIQ_BIN_FOLDER or MAJIQ_SOFTWARE_ROOT + "/env/bin"
MAJIQ_GFF3 = MODULE_REFERENCE["majiq"] + "/annotation.gff3"
MAJIQ_SPLICEGRAPH = MODULE_REFERENCE["majiq"] + "/splicegraph.zarr"


def majiq_source(wildcards):
    if not config.get("majiq_source"):
        raise WorkflowError("MAJIQ v3 is licensed: download its source after registering with "
                            "BioCiphers and set majiq_source (or majiq_bin_folder for an existing "
                            "installation)")
    return config["majiq_source"]


rule umbrella_majiq_software:
    input:
        source=majiq_source,
        env="envs/umbrella-majiq.yaml"
    output:
        protected(MAJIQ_SOFTWARE_ROOT + "/installed.json")
    params:
        prefix=MAJIQ_SOFTWARE_ROOT + "/env",
        conda=config.get("conda_executable", "conda")
    log:
        MAJIQ_SOFTWARE_ROOT + "/install.log"
    shell:
        "python3 src/install_majiq.py --source {input.source} --env-file {input.env} "
        "--prefix {params.prefix} --record {output} --conda {params.conda} > {log} 2>&1"


rule umbrella_majiq_reference:
    input:
        gtf=UMBRELLA_GTF["whippet"],
        software=MAJIQ_INSTALLED
    output:
        gff3=protected(MAJIQ_GFF3),
        splicegraph=protected(directory(MAJIQ_SPLICEGRAPH))
    params:
        majiq_build=MAJIQ_BIN + "/majiq-build"
    log:
        MODULE_REFERENCE["majiq"] + "/splicegraph.log"
    conda:
        "../envs/umbrella-quant.yaml"
    shell:
        # the microexon-inserted Whippet GTF, so microexons are in the graph
        "gzip -dcf {input.gtf} | gffread - --keep-genes -o {output.gff3} "
        "&& {params.majiq_build} gff3 {output.gff3} {output.splicegraph} > {log} 2>&1"


if MODULE_PRODUCERS:
    rule umbrella_majiq_sj:
        input:
            bam=module_ingest("majiq", lambda w: UMBRELLA_WORK.format(**w) + "/align/aligned.bam"),
            bai=module_ingest("majiq", lambda w: UMBRELLA_WORK.format(**w) + "/align/aligned.bam.bai"),
            splicegraph=module_ingest("majiq", MAJIQ_SPLICEGRAPH),
            software=module_ingest("majiq", MAJIQ_INSTALLED)
        output:
            # run.sj next to it is MAJIQ's own format (file or directory)
            protected(MODULE_CACHE.replace("{tool}", "majiq") + "/cache.json")
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", "majiq"),
            strandedness=lambda w: config.get("majiq_strandness") or umbrella_quant_run(w).strandedness,
            bin=MAJIQ_BIN,
            identity=lambda w: module_identity(w, "majiq")
        threads: 4
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.majiq.{cache_id}.tsv"
        shell:
            "python3 src/run_majiq.py sj --bam {input.bam} --splicegraph {input.splicegraph} "
            "--out-dir {params.out} --run-id {wildcards.run_id} --strandedness {params.strandedness} "
            "--threads {threads} --bin-dir {params.bin} --identity {params.identity:q}"


rule umbrella_majiq_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches("majiq"),
        splicegraph=MAJIQ_SPLICEGRAPH,
        software=MAJIQ_INSTALLED
    output:
        native=MODULE_ANALYSIS.replace("{tool}", "majiq") + "/native.tsv.gz",
        results=MODULE_ANALYSIS.replace("{tool}", "majiq") + "/results.tsv.gz",
        status=MODULE_ANALYSIS.replace("{tool}", "majiq") + "/status.json"
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "majiq"),
        caches=lambda w: module_cache_arguments(w, "majiq"),
        bin=MAJIQ_BIN,
        args=_json.dumps({key: value for key, value in MODULE_SETTINGS["majiq"].items()
                          if value and key in ("update", "psicov", "heterogen")}),
        statistic=MODULE_SETTINGS["majiq"]["statistic"],
        ids=lambda w: module_ids(w, "majiq", SOFTWARE["majiq"])
    threads: 4
    benchmark:
        MODULE_ANALYSIS.replace("{tool}", "majiq") + "/benchmark.tsv"
    shell:
        "python3 src/run_majiq.py compare --preflight {input.preflight} {params.caches} "
        "--splicegraph {input.splicegraph} --out-dir {params.out} --threads {threads} "
        "--bin-dir {params.bin} --args {params.args:q} --statistic {params.statistic} --ids {params.ids:q}"


# ------------------------------------------------------------ targets

MODULE_FINAL = {"majiq": "results.tsv.gz", "dapars2": "results.tsv.gz", "qapa": "results.tsv.gz"}


def module_prepare_targets(tool):
    # MAJIQ's run.sj (file or directory) is checked by module_check, not declared
    return [path for record in UMBRELLA_MANIFEST.included_runs()
            for path in (module_cache_outputs(record, tool) if tool != "majiq"
                         else [module_cache_dir(record, tool) + "/cache.json"])]


def module_quant_targets(tool):
    return module_prepare_targets(tool) + [
        module_analysis_dir(result, tool) + "/" + MODULE_FINAL[tool]
        for _, result in sorted(UMBRELLA_COMPARISONS.items())]


def module_check(tool):
    """Target-scoped cache checks, run only when a target asks for this tool.

    A cache whose marker (cache.json) exists must have every other file too;
    in cached_only, every run must have a cache. Other tools' caches are never
    looked at, so a missing MAJIQ cache cannot block a QAPA target."""
    runs = UMBRELLA_MANIFEST.included_runs()
    incomplete = [module_cache_dir(record, tool) for record in runs
                  if _os.path.exists(module_cache_dir(record, tool) + "/cache.json")
                  and not all(_os.path.exists(module_cache_dir(record, tool) + "/" + name)
                              for name in CACHE_FILES[tool])]
    if incomplete:
        raise WorkflowError("{} {} cache(s) have their marker but not all their files (first: {}); "
                            "delete those folders and ingest them again".format(
                                len(incomplete), tool, incomplete[0]))
    if UMBRELLA_MODULES_MODE == "cached_only":
        missing = missing_caches(runs, module_cache_dir, [tool])
        if missing:
            raise WorkflowError(
                "umbrella_modules_mode is cached_only, but {} run/tool caches are missing (first: run {}, "
                "{}). Ingest them once with umbrella_modules_mode: ingest, plus umbrella_allow_restage: "
                "true if their reads are gone.".format(len(missing), missing[0][0], missing[0][1]))


def module_target(tool, targets):
    def function(wildcards):
        module_check(tool)
        return targets(tool)
    return function


def umbrella_module_targets(wildcards=None):
    """The selected modules' targets (for quant_umbrella and quant_umbrella_modules)."""
    paths = []
    for tool in UMBRELLA_MODULES:
        module_check(tool)
        paths += module_quant_targets(tool)
    return paths


rule prepare_majiq:
    input: module_target("majiq", module_prepare_targets)

rule prepare_dapars2:
    input: module_target("dapars2", module_prepare_targets)

rule prepare_qapa:
    input: module_target("qapa", module_prepare_targets)

rule quant_majiq:
    input: module_target("majiq", module_quant_targets)

rule quant_dapars2:
    input: module_target("dapars2", module_quant_targets)

rule quant_qapa:
    input: module_target("qapa", module_quant_targets)


UMBRELLA_MODULE_INVENTORY = "umbrella/modules/inventory/{}.{}.json".format(
    UMBRELLA_REFERENCE_ID, "_".join(UMBRELLA_MODULES) or "none")


def module_inventory_inputs(wildcards):
    if not UMBRELLA_MODULES:
        raise WorkflowError("quant_umbrella_modules needs umbrella_modules, e.g. "
                            "umbrella_modules: [majiq, dapars2, qapa]; or use quant_majiq, "
                            "quant_dapars2 or quant_qapa")
    return umbrella_module_targets(wildcards)


localrules: umbrella_module_inventory


rule umbrella_module_inventory:
    input:
        module_inventory_inputs
    output:
        UMBRELLA_MODULE_INVENTORY
    run:
        inventory = {"reference_id": UMBRELLA_REFERENCE_ID, "modules": UMBRELLA_MODULES,
                     "mode": UMBRELLA_MODULES_MODE,
                     "module_reference_ids": {tool: MODULE_REFERENCE_ID[tool] for tool in UMBRELLA_MODULES},
                     "software": {tool: SOFTWARE[tool] for tool in UMBRELLA_MODULES},
                     "comparisons": {}}
        for comparison_id, result in sorted(UMBRELLA_COMPARISONS.items()):
            for tool in UMBRELLA_MODULES:
                folder = module_analysis_dir(result, tool)
                status = _json.loads(Path(folder + "/status.json").read_text())
                inventory["comparisons"].setdefault(comparison_id, {})[tool] = {
                    "folder": folder, "status": status["status"], "reasons": status.get("reasons", []),
                    "bytes": sum(p.stat().st_size for p in Path(folder).rglob("*") if p.is_file())}
        inventory["cache_bytes"] = {tool: sum(
            Path(path).stat().st_size for record in UMBRELLA_MANIFEST.included_runs()
            for path in module_prepare_targets(tool) if record.run_id in path) for tool in UMBRELLA_MODULES}
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text(_json.dumps(inventory, sort_keys=True, indent=2) + "\n")


rule quant_umbrella_modules:
    input:
        UMBRELLA_MODULE_INVENTORY

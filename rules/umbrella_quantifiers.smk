"""Manifest-routed native quantifiers and guarded group Salmon reduction."""

from src.reduce_salmon import reduce_group
from src.whippet_flags import read_paths as whippet_read_paths
def umbrella_quant_run(wildcards):
    record = UMBRELLA_MANIFEST.by_run.get(wildcards.run_id)
    if record is None or not record.include:
        raise WorkflowError("umbrella quantification run is absent or excluded")
    if (record.reference_id, record.project_id, record.batch_id) != (
            wildcards.reference_id, wildcards.project_id, wildcards.batch_id):
        raise WorkflowError("umbrella quantification path conflicts with manifest")
    return record


def umbrella_quant_reads(wildcards):
    return UMBRELLA_MANIFEST.native_reads(umbrella_quant_run(wildcards).run_id)


def umbrella_read_arguments(wildcards, tool):
    record = umbrella_quant_run(wildcards)
    reads = UMBRELLA_MANIFEST.native_reads(record.run_id)
    if tool == "salmon":
        return "-1 {} -2 {}".format(*reads) if record.layout == "PE" else "-r {}".format(reads[0])
    return " ".join(reads)


# Whippet's per-run outputs are kept and protected: junction counts and gene
# and isoform TPM exist nowhere else, and remaking any of them needs the reads.
# Explicit old-library backfill stages a separate copy of the reads. It never
# changes the shared validation marker or input paths of ME/Salmon producers.
WHIPPET_BACKFILL = "umbrella/whippet_backfill/{reference_id}/{project_id}/{batch_id}/{run_id}"
import re as _whippet_re
WHIPPET_BACKFILL_PATTERNS = {
    layout: "(?:" + "|".join(_whippet_re.escape(r.run_id) for r in UMBRELLA_MANIFEST.included_runs()
                                if r.layout == layout and r.metadata.get("whippet_read_source") == "backfill") + ")"
    for layout in ("PE", "SE")}


def whippet_backfill_run(wildcards):
    record = umbrella_quant_run(wildcards)
    if record.metadata.get("whippet_read_source") != "backfill":
        raise WorkflowError("Whippet backfill staging requires whippet_read_source=backfill")
    return record


rule umbrella_whippet_backfill_pe:
    input:
        sources=umbrella_sources
    output:
        r1=temp(WHIPPET_BACKFILL + "/R1.fastq.gz"),
        r2=temp(WHIPPET_BACKFILL + "/R2.fastq.gz")
    wildcard_constraints:
        run_id=WHIPPET_BACKFILL_PATTERNS["PE"]
    params:
        identity=lambda w: UMBRELLA_MANIFEST.run_sha256(whippet_backfill_run(w).run_id),
        tmpdir=UMBRELLA_STAGE_TMPDIR
    threads: 4
    resources:
        get_data=1
    conda:
        "../envs/umbrella-inputs.yaml"
    shell:
        "python3 src/umbrella_stage_reads.py --manifest {UMBRELLA_MANIFEST.path:q} --run-id {wildcards.run_id:q} "
        "--r1 {output.r1:q} --r2 {output.r2:q} --threads {threads} {params.tmpdir} --allow-restage"


rule umbrella_whippet_backfill_se:
    input:
        sources=umbrella_sources
    output:
        r1=temp(WHIPPET_BACKFILL + "/R1.fastq.gz")
    wildcard_constraints:
        run_id=WHIPPET_BACKFILL_PATTERNS["SE"]
    params:
        identity=lambda w: UMBRELLA_MANIFEST.run_sha256(whippet_backfill_run(w).run_id),
        tmpdir=UMBRELLA_STAGE_TMPDIR
    threads: 4
    resources:
        get_data=1
    conda:
        "../envs/umbrella-inputs.yaml"
    shell:
        "python3 src/umbrella_stage_reads.py --manifest {UMBRELLA_MANIFEST.path:q} --run-id {wildcards.run_id:q} "
        "--r1 {output.r1:q} --threads {threads} {params.tmpdir} --allow-restage"


rule umbrella_whippet_backfill_validate:
    input:
        reads=lambda w: whippet_read_paths(whippet_backfill_run(w))
    output:
        protected(WHIPPET_BACKFILL + "/reads.valid")
    shell:
        "python3 src/validate_fastq_pairs.py {input.reads:q} --marker {output:q} --run-id {wildcards.run_id:q}"


ruleorder: umbrella_whippet_backfill_pe > umbrella_whippet_backfill_se


rule umbrella_whippet_quant:
    input:
        reads=lambda w: whippet_read_paths(umbrella_quant_run(w)),
        valid=lambda w: whippet_read_paths(umbrella_quant_run(w))[0].rsplit("/", 1)[0] + "/reads.valid",
        index=UMBRELLA_WHIPPET_MEMBERS[0],
        reference=UMBRELLA_REFERENCE_MANIFEST
    output:
        gene=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.gene.tpm.gz"),
        isoform=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.isoform.tpm.gz"),
        jnc=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.jnc.gz"),
        mapping=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.map.gz"),
        psi=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.psi.gz"),
        provenance=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.provenance.json")
    params:
        julia=config.get("julia", "julia"),
        whippet_bin=config.get("whippet_bin_folder", ""),
        reads=lambda w: " ".join(whippet_read_paths(umbrella_quant_run(w))),
        prefix="umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant",
        processing_sha256=lambda w: UMBRELLA_MANIFEST.run_sha256(umbrella_quant_run(w).run_id),
        flags=config.get("whippet_flags", "")
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.whippet.log"
    shell:
        "python3 src/whippet_core.py quant --reads {params.reads} "
        "--index {input.index:q} --prefix {params.prefix:q} --run-id {wildcards.run_id:q} "
        "--reference-id {wildcards.reference_id:q} --processing-sha256 {params.processing_sha256:q} --julia {params.julia:q} "
        "--whippet-bin {params.whippet_bin:q} --flags={params.flags:q} > {log} 2>&1"


rule umbrella_salmon_quant:
    input:
        reads=umbrella_quant_reads,
        valid=lambda w: umbrella_quant_run(w).work_dir + "/reads.valid",
        index=UMBRELLA_SALMON_INDEX,
        reference=UMBRELLA_REFERENCE_MANIFEST
    output:
        temp(directory("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/salmon"))
    params:
        reads=lambda w: umbrella_read_arguments(w, "salmon"),
        library=UMBRELLA_REFERENCE.get("salmon_library_type", "A")
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.salmon.log"
    threads: 4
    conda:
        "../envs/umbrella-quant.yaml"
    shell:
        "salmon quant -i {input.index} -l {params.library} {params.reads} "
        "--validateMappings -p {threads} -o {output} 2> {log}"


def umbrella_group_runs(wildcards):
    if wildcards.reference_id != UMBRELLA_REFERENCE_ID:
        raise WorkflowError("Salmon shard reference_id conflicts with frozen bundle")
    records = UMBRELLA_MANIFEST.runs_for(
        wildcards.project_id, wildcards.group, wildcards.batch_id)
    if not records:
        raise WorkflowError("Salmon shard has no included runs")
    return records


def umbrella_group_quants(wildcards):
    return [record.work_dir + "/salmon" for record in umbrella_group_runs(wildcards)]


rule umbrella_salmon_shard:
    input:
        quants=umbrella_group_quants,
        tx2gene=UMBRELLA_TX2GENE,
        reference=UMBRELLA_REFERENCE_MANIFEST
    output:
        counts=protected("genes/{reference_id}/{project_id}/{group}/{batch_id}.salmon_counts.tsv.gz"),
        tpm=protected("genes/{reference_id}/{project_id}/{group}/{batch_id}.salmon_tpm.tsv.gz"),
        length=protected("genes/{reference_id}/{project_id}/{group}/{batch_id}.salmon_length.tsv.gz"),
        tx_counts=protected("genes/{reference_id}/{project_id}/{group}/{batch_id}.salmon_tx_counts.tsv.gz"),
        tx_tpm=protected("genes/{reference_id}/{project_id}/{group}/{batch_id}.salmon_tx_tpm.tsv.gz"),
        checksums=protected("genes/{reference_id}/{project_id}/{group}/{batch_id}.salmon_checksums.json")
    params:
        manifest_sha256=umbrella_shard_sha
    run:
        records = umbrella_group_runs(wildcards)
        quant_by_run = {record.run_id: record.work_dir + "/salmon/quant.sf"
                        for record in records}
        prefix = "genes/{}/{}/{}/{}".format(
            wildcards.reference_id, wildcards.project_id,
            wildcards.group, wildcards.batch_id)
        reduce_group(quant_by_run, input.tx2gene, prefix,
                     wildcards.reference_id, params.manifest_sha256)


UMBRELLA_WHIPPET_PSI = [record.work_dir + "/whippet/quant.psi.gz"
                       for record in UMBRELLA_MANIFEST.included_runs()]
UMBRELLA_SALMON_SHARDS = []
for reference_id, project_id, group, batch_id in sorted({
        (run.reference_id, run.project_id, run.group, run.batch_id)
        for run in UMBRELLA_MANIFEST.included_runs()}):
    prefix = "genes/{}/{}/{}/{}".format(reference_id, project_id, group, batch_id)
    UMBRELLA_SALMON_SHARDS.extend(prefix + ".salmon_" + suffix for suffix in
                                  ("counts.tsv.gz", "tpm.tsv.gz", "length.tsv.gz",
                                   "tx_counts.tsv.gz", "tx_tpm.tsv.gz",
                                   "checksums.json"))


# Public targets. quant and get_whippet_psi keep their legacy meaning.
UMBRELLA_ME_QUANT_INPUTS = [path for path in rules.quant.input
                           if str(path) != str(FILTERED_ME_OUTPUT)]


rule quant_microexonator:
    input:
        UMBRELLA_ME_QUANT_INPUTS,
        UMBRELLA_SPLICING_SHARDS["microexonator"]


rule quant_whippet:
    input:
        UMBRELLA_ME_QUANT_INPUTS,
        UMBRELLA_SPLICING_SHARDS["microexonator"],
        UMBRELLA_WHIPPET_PSI,
        UMBRELLA_SPLICING_SHARDS["whippet"]


def umbrella_selected(tool, paths):
    return paths if UMBRELLA_OPTIONAL[tool] else []


# The full umbrella; umbrella_optional switches tools off (rules/umbrella_alignment.smk).
rule quant_umbrella:
    input:
        umbrella_selected("microexonator", UMBRELLA_ME_QUANT_INPUTS),
        umbrella_selected("microexonator", ["Report/out.robustly_detected.txt"]),
        umbrella_selected("whippet", UMBRELLA_WHIPPET_PSI),
        umbrella_selected("salmon", UMBRELLA_SALMON_SHARDS),
        UMBRELLA_ALIGNMENT_TARGETS,
        umbrella_selected("microexonator", UMBRELLA_SPLICING_SHARDS["microexonator"]),
        umbrella_selected("whippet", UMBRELLA_SPLICING_SHARDS["whippet"]),
        UMBRELLA_COMPARISON_TARGETS,
        UMBRELLA_SYNTHESIS_TARGETS,
        UMBRELLA_PCA_TARGETS,
        UMBRELLA_QC_TARGETS,
        # the selected modules, through their inventory (which needs all their outputs)
        [UMBRELLA_MODULE_INVENTORY] if UMBRELLA_MODULES else [],
        [UMBRELLA_REFERENCE_MANIFEST] if any(UMBRELLA_OPTIONAL[tool] for tool in ("whippet", "salmon", "hisat2")) else []


# differential_inclusion: MicroExonator per run, its group shards and its own
# Whippet-free delta (src/me_delta.py) for every comparison in
# umbrella_comparisons. Needs no HISAT2, Salmon or Whippet index.
# delta_method: whippet gives Whippet's delta instead.
UMBRELLA_DELTA_FILE = {"microexonator": "microexonator_delta.tsv", "whippet": "whippet_delta.diff.gz"}
UMBRELLA_DELTA_TARGETS = [
    "comparisons/{}/{}/{}/{}".format(result["reference_id"], result["project_id"], comparison_id,
                                     UMBRELLA_DELTA_FILE[DELTA_METHOD])
    for comparison_id, result in sorted(UMBRELLA_COMPARISONS.items())]


def umbrella_differential_inclusion(wildcards):
    if not UMBRELLA_COMPARISONS:
        raise WorkflowError("differential_inclusion needs umbrella_comparisons (a comparisons YAML)")
    if not UMBRELLA_OPTIONAL[DELTA_METHOD]:
        raise WorkflowError("differential_inclusion with delta_method: {0} needs {0}, which "
                            "umbrella_optional switches off".format(DELTA_METHOD))
    if DELTA_METHOD == "microexonator":
        return UMBRELLA_ME_QUANT_INPUTS + UMBRELLA_SPLICING_SHARDS["microexonator"] + UMBRELLA_DELTA_TARGETS
    return UMBRELLA_WHIPPET_PSI + UMBRELLA_SPLICING_SHARDS["whippet"] + UMBRELLA_DELTA_TARGETS


rule differential_inclusion:
    input:
        umbrella_differential_inclusion

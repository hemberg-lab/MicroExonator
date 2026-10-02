"""Manifest-routed native quantifiers and guarded group Salmon reduction."""

from src.reduce_salmon import reduce_group
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
rule umbrella_whippet_quant:
    input:
        reads=umbrella_quant_reads,
        valid=lambda w: umbrella_quant_run(w).work_dir + "/reads.valid",
        index=UMBRELLA_WHIPPET_MEMBERS[0],
        reference=UMBRELLA_REFERENCE_MANIFEST
    output:
        gene=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.gene.tpm.gz"),
        isoform=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.isoform.tpm.gz"),
        jnc=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.jnc.gz"),
        mapping=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.map.gz"),
        psi=protected("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant.psi.gz")
    params:
        julia=config.get("julia", "julia"),
        whippet_bin=config.get("whippet_bin_folder", ""),
        reads=lambda w: umbrella_read_arguments(w, "whippet"),
        prefix="umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/whippet/quant",
        flags=config.get("whippet_flags", "")
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.whippet.log"
    shell:
        "{params.julia} {params.whippet_bin}/whippet-quant.jl {params.reads} "
        "-x {input.index} -o {params.prefix} {params.flags} 2> {log}"


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
        UMBRELLA_QC_TARGETS,
        umbrella_module_targets,
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

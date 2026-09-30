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


rule quant_umbrella:
    input:
        UMBRELLA_ME_QUANT_INPUTS,
        "Report/out.robustly_detected.txt",
        UMBRELLA_WHIPPET_PSI,
        UMBRELLA_SALMON_SHARDS,
        UMBRELLA_ALIGNMENT_TARGETS,
        UMBRELLA_SPLICING_SHARDS["microexonator"],
        UMBRELLA_SPLICING_SHARDS["whippet"],
        UMBRELLA_COMPARISON_TARGETS,
        UMBRELLA_SYNTHESIS_TARGETS,
        UMBRELLA_QC_TARGETS,
        UMBRELLA_MODULE_TARGETS,
        UMBRELLA_REFERENCE_MANIFEST

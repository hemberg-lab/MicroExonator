"""Lightweight read QC and one MultiQC report per project (umbrella_optional: multiqc).

For each run, about 50 KB is kept next to the QC shards, under
qc/{reference_id}/{project_id}/{group}/{batch_id}/runs/{run_id}/:
FastQC's fastqc_data.txt for the first `umbrella_fastqc_reads` reads of each
file (default 2,000,000; 0 for all), and the HISAT2, featureCounts, Salmon and
Whippet summaries, copied before their temporary sources are deleted
(src/umbrella_fastqc.py). No FASTQ or BAM is kept.

New batches only: a run whose group-by-batch QC shard was written before these
summaries existed has no reads left, so it never triggers FastQC (that would
mean downloading it again). The report lists such runs in
runs_without_qc.txt. Since the inputs are kept, a report can be rebuilt at any
time.
"""

import os

UMBRELLA_QC_RUN = "qc/" + UMBRELLA_SHARD + "/runs/{run_id}"
UMBRELLA_FASTQC_READS = int(config.get("umbrella_fastqc_reads", 2000000))
UMBRELLA_MULTIQC = "multiqc/{reference_id}/{project_id}"
UMBRELLA_QC_KEPT = ("{run_id}.hisat2.summary.txt", "{run_id}.featurecounts.txt.summary",
                    "salmon/{run_id}/aux_info/meta_info.json",
                    "salmon/{run_id}/lib_format_counts.json", "{run_id}.whippet.map.gz")


def umbrella_qc_root(record):
    return "qc/{}/{}/{}/{}/runs/{}".format(record.reference_id, record.project_id,
                                            record.group, record.batch_id, record.run_id)


def umbrella_qc_is_new(record):
    """True while the run's QC shard is not written yet, i.e. its reads still exist."""
    return not os.path.exists("qc/{}/{}/{}/{}.qc.tsv.gz".format(
        record.reference_id, record.project_id, record.group, record.batch_id))


def umbrella_qc_outputs(record):
    root = umbrella_qc_root(record)
    return [root + "/fastqc"] + [root + "/" + name.format(run_id=record.run_id)
                                 for name in UMBRELLA_QC_KEPT]


rule umbrella_fastqc:
    input:
        reads=lambda w: UMBRELLA_MANIFEST.native_reads(umbrella_align_run(w).run_id),
        valid=lambda w: umbrella_align_run(w).work_dir + "/reads.valid"
    output:
        directory(UMBRELLA_QC_RUN + "/fastqc")
    wildcard_constraints:
        reference_id="[^/]+", project_id="[^/]+", group="[^/]+", batch_id="[^/]+", run_id="[^/]+"
    params:
        reads=UMBRELLA_FASTQC_READS,
        scratch=lambda w: umbrella_align_run(w).work_dir
    threads: 2
    conda:
        "../envs/umbrella-qc.yaml"
    shell:
        "python3 src/umbrella_fastqc.py fastqc --run-id {wildcards.run_id} --reads {input.reads} "
        "--out {output} --reads-per-file {params.reads} --threads {threads} --scratch {params.scratch}"


rule umbrella_qc_keep:
    input:
        hisat2=lambda w: umbrella_align_run(w).work_dir + "/align/hisat2.summary.txt",
        featurecounts=lambda w: umbrella_align_run(w).work_dir + "/align/featurecounts.txt.summary",
        salmon=lambda w: umbrella_align_run(w).work_dir + "/salmon",
        whippet_map=lambda w: umbrella_align_run(w).work_dir + "/whippet/quant.map.gz"
    output:
        [UMBRELLA_QC_RUN + "/" + name for name in UMBRELLA_QC_KEPT]
    wildcard_constraints:
        reference_id="[^/]+", project_id="[^/]+", group="[^/]+", batch_id="[^/]+", run_id="[^/]+"
    params:
        out=UMBRELLA_QC_RUN
    shell:
        "python3 src/umbrella_fastqc.py keep --run-id {wildcards.run_id} --hisat2 {input.hisat2} "
        "--featurecounts {input.featurecounts} --salmon {input.salmon} "
        "--whippet-map {input.whippet_map} --out {params.out}"


def umbrella_project_runs(wildcards):
    return [record for record in UMBRELLA_MANIFEST.included_runs()
            if (record.reference_id, record.project_id) == (wildcards.reference_id, wildcards.project_id)]


def umbrella_multiqc_inputs(wildcards):
    paths = []
    for record in umbrella_project_runs(wildcards):
        if umbrella_qc_is_new(record):
            paths += umbrella_qc_outputs(record)
        elif os.path.isdir(umbrella_qc_root(record)):
            # kept earlier: the folder as a plain input, so no rule is asked to remake it
            paths.append(umbrella_qc_root(record))
    return paths


def umbrella_multiqc_missing(wildcards):
    return " ".join(record.run_id for record in umbrella_project_runs(wildcards)
                    if not umbrella_qc_is_new(record) and not os.path.isdir(umbrella_qc_root(record)))


rule umbrella_multiqc:
    input:
        kept=umbrella_multiqc_inputs,
        config="src/umbrella_multiqc_config.yaml"
    output:
        report=UMBRELLA_MULTIQC + "/multiqc_report.html",
        data=directory(UMBRELLA_MULTIQC + "/multiqc_data"),
        missing=UMBRELLA_MULTIQC + "/runs_without_qc.txt"
    wildcard_constraints:
        reference_id="[^/]+", project_id="[^/]+"
    params:
        out=UMBRELLA_MULTIQC,
        missing=umbrella_multiqc_missing
    log:
        UMBRELLA_MULTIQC + "/multiqc.log"
    conda:
        "../envs/umbrella-qc.yaml"
    shell:
        "find {input.kept} -type f | sort > {params.out}/files.txt "
        "&& multiqc -f -o {params.out} -n multiqc_report.html -c {input.config} "
        "--data-format tsv -l {params.out}/files.txt > {log} 2>&1 "
        "&& printf '%s\\n' '# runs without kept QC summaries (processed before they were kept)' "
        "{params.missing} > {output.missing}"


UMBRELLA_QC_TARGETS = []
if UMBRELLA_OPTIONAL.get("multiqc"):
    _projects = set()
    for record in UMBRELLA_MANIFEST.included_runs():
        if umbrella_qc_is_new(record):
            UMBRELLA_QC_TARGETS += umbrella_qc_outputs(record)
            _projects.add((record.reference_id, record.project_id))
        elif os.path.isdir(umbrella_qc_root(record)):
            _projects.add((record.reference_id, record.project_id))
    UMBRELLA_QC_TARGETS += [UMBRELLA_MULTIQC.format(reference_id=reference_id, project_id=project_id)
                            + "/multiqc_report.html" for reference_id, project_id in sorted(_projects)]

"""Opt-in, mate-preserving intake with a legacy FASTQ compatibility path."""

import re
import shlex
import subprocess

from src.umbrella_manifest import link_legacy_fastq


def umbrella_run(wildcards):
    run = UMBRELLA_MANIFEST.by_run.get(wildcards.run_id)
    if run is None or not run.include:
        raise WorkflowError("umbrella run_id is absent or excluded: {}".format(wildcards.run_id))
    actual = (run.reference_id, run.project_id, run.batch_id)
    requested = (wildcards.reference_id, wildcards.project_id, wildcards.batch_id)
    if actual != requested:
        raise WorkflowError("umbrella run_id {} metadata does not match staged path".format(run.run_id))
    return run


def umbrella_sources(wildcards):
    run = umbrella_run(wildcards)
    if run.source_type == "sra":
        return []
    return [path for path in (run.source_1, run.source_2) if path]


# `umbrella_stage_tmpdir`: scratch for uncompressed reads while they are gzipped
# (node-local disk keeps them off the shared file system); default: the run's
# work directory.
UMBRELLA_STAGE_TMPDIR = ("--tmpdir " + shlex.quote(str(config["umbrella_stage_tmpdir"]))
                         if config.get("umbrella_stage_tmpdir") else "")


PE_RUN_PATTERN = "(?:{})".format("|".join(
    re.escape(run.run_id) for run in UMBRELLA_MANIFEST.included_runs() if run.layout == "PE"))
SE_RUN_PATTERN = "(?:{})".format("|".join(
    re.escape(run.run_id) for run in UMBRELLA_MANIFEST.included_runs() if run.layout == "SE"))


rule umbrella_stage_reads:
    input:
        sources=umbrella_sources
    output:
        # temporary: removed once every tool that reads them has run
        r1=temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/R1.fastq.gz"),
        r2=temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/R2.fastq.gz")
    wildcard_constraints:
        run_id=PE_RUN_PATTERN
    resources:
        # like the legacy download rules: `--resources get_data=N` caps parallel SRA downloads
        get_data=lambda wildcards: 1 if umbrella_run(wildcards).source_type == "sra" else 0
    threads: 6
    # SRA downloads fail now and then; a retry starts the run again from scratch
    retries: 2
    params:
        tmpdir=UMBRELLA_STAGE_TMPDIR,
        manifest=str(UMBRELLA_MANIFEST.path),
        run_sha256=lambda w: UMBRELLA_MANIFEST.run_sha256(w.run_id)
    conda:
        "../envs/umbrella-inputs.yaml"
    shell:
        "python3 src/umbrella_stage_reads.py --manifest {params.manifest:q} "
        "--run-id {wildcards.run_id:q} --r1 {output.r1:q} --r2 {output.r2:q} "
        "--threads {threads} {params.tmpdir}"


rule umbrella_stage_single:
    input:
        sources=umbrella_sources
    output:
        temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/R1.fastq.gz")
    wildcard_constraints:
        run_id=SE_RUN_PATTERN
    resources:
        # like the legacy download rules: `--resources get_data=N` caps parallel SRA downloads
        get_data=lambda wildcards: 1 if umbrella_run(wildcards).source_type == "sra" else 0
    threads: 6
    # SRA downloads fail now and then; a retry starts the run again from scratch
    retries: 2
    params:
        tmpdir=UMBRELLA_STAGE_TMPDIR,
        manifest=str(UMBRELLA_MANIFEST.path),
        run_sha256=lambda w: UMBRELLA_MANIFEST.run_sha256(w.run_id)
    conda:
        "../envs/umbrella-inputs.yaml"
    shell:
        "python3 src/umbrella_stage_reads.py --manifest {params.manifest:q} "
        "--run-id {wildcards.run_id:q} --r1 {output:q} --threads {threads} {params.tmpdir}"


def umbrella_validation_inputs(wildcards):
    record = umbrella_run(wildcards)
    return UMBRELLA_MANIFEST.native_reads(record.run_id)


rule umbrella_validate_reads:
    input:
        reads=umbrella_validation_inputs
    output:
        "umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/reads.valid"
    params:
        run_sha256=lambda w: UMBRELLA_MANIFEST.run_sha256(w.run_id)
    run:
        record = umbrella_run(wildcards)
        subprocess.run(["python3", "src/validate_fastq_pairs.py", *[str(path) for path in input.reads],
                        "--marker", str(output[0]), "--run-id", record.run_id], check=True)


rule umbrella_legacy_fastq:
    input:
        reads=umbrella_validation_inputs,
        valid="umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/reads.valid"
    output:
        temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/legacy.fastq.gz")
    threads: 4
    run:
        record = umbrella_run(wildcards)
        if record.layout != "PE":
            raise WorkflowError("single-end runs use native R1 as legacy input")
        subprocess.run(["bash", "src/concat_mates.sh", str(input.reads[0]),
                        str(input.reads[1]), str(output[0]), str(threads)], check=True)


def umbrella_legacy_bridge_input(wildcards):
    run = UMBRELLA_MANIFEST.by_run.get(wildcards.sample)
    if run is None or not run.include:
        raise WorkflowError("umbrella sample is absent or excluded: {}".format(wildcards.sample))
    return UMBRELLA_MANIFEST.legacy_fastq(run.run_id)


def umbrella_legacy_bridge_validation(wildcards):
    run = UMBRELLA_MANIFEST.by_run[wildcards.sample]
    return run.work_dir + "/reads.valid"


rule umbrella_legacy_bridge:
    input:
        fastq=umbrella_legacy_bridge_input,
        valid=umbrella_legacy_bridge_validation
    output:
        # like the legacy download rules: temporary unless Keep_fastq_gz
        "FASTQ/{sample}.fastq.gz" if str2bool(config.get("Keep_fastq_gz", False))
        else temp("FASTQ/{sample}.fastq.gz")
    run:
        # PE: the concatenation is temporary; SE: the staged R1 is temporary.
        # Either way MicroExonator gets a path that outlives the staged file.
        link_legacy_fastq(input.fastq, output[0])


rule umbrella_intake:
    input:
        [run.work_dir + "/reads.valid" for run in UMBRELLA_MANIFEST.included_runs()]

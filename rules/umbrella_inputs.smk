"""Opt-in, mate-preserving intake. Only paths below umbrella/work are outputs."""

import bz2
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


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


def stage_local_fastq(source, destination):
    """Link compressed inputs; make a workflow-owned gzip for other FASTQs."""
    source = Path(source)
    destination = Path(destination)
    if source.suffix == ".gz":
        destination.symlink_to(source)
        return
    opener = bz2.open if source.suffix == ".bz2" else open
    with opener(source, "rb") as incoming, gzip.open(destination, "wb") as outgoing:
        shutil.copyfileobj(incoming, outgoing)


def compress_fastq(source, destination):
    with open(source, "rb") as incoming, gzip.open(destination, "wb") as outgoing:
        shutil.copyfileobj(incoming, outgoing)


PE_RUN_PATTERN = "(?:{})".format("|".join(
    re.escape(run.run_id) for run in UMBRELLA_MANIFEST.included_runs() if run.layout == "PE"))
SE_RUN_PATTERN = "(?:{})".format("|".join(
    re.escape(run.run_id) for run in UMBRELLA_MANIFEST.included_runs() if run.layout == "SE"))


rule umbrella_stage_reads:
    input:
        umbrella_sources
    output:
        r1=temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/R1.fastq.gz"),
        r2=temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/R2.fastq.gz")
    wildcard_constraints:
        run_id=PE_RUN_PATTERN
    run:
        record = umbrella_run(wildcards)
        if record.layout != "PE":
            raise WorkflowError("umbrella_stage_reads requires PE; use umbrella_stage_single")
        directory = Path(output.r1).parent
        directory.mkdir(parents=True, exist_ok=True)
        if record.source_type == "fastq":
            stage_local_fastq(record.source_1, output.r1)
            stage_local_fastq(record.source_2, output.r2)
        else:
            with tempfile.TemporaryDirectory(dir=str(directory)) as work:
                work = Path(work)
                if record.source_type == "sra":
                    subprocess.run(["fasterq-dump", "--split-files", "-O", str(work),
                                    record.source_1], check=True)
                    mates = [work / (record.source_1 + suffix) for suffix in
                             ("_1.fastq", "_2.fastq")]
                else:
                    mates = [work / "R1.fastq", work / "R2.fastq"]
                    sorted_bam = work / "name_sorted.bam"
                    subprocess.run(["samtools", "sort", "-n", "-o", str(sorted_bam),
                                    record.source_1], check=True)
                    subprocess.run(["samtools", "fastq", "-1", str(mates[0]), "-2",
                                    str(mates[1]), "-0", os.devnull, "-s", os.devnull,
                                    "-n", "-F", "0x900", str(sorted_bam)], check=True)
                for source, target in zip(mates, (output.r1, output.r2)):
                    if not source.is_file():
                        raise WorkflowError("missing mate from {}: {}".format(record.run_id, source))
                    compress_fastq(source, target)


rule umbrella_stage_single:
    input:
        umbrella_sources
    output:
        temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/R1.fastq.gz")
    wildcard_constraints:
        run_id=SE_RUN_PATTERN
    run:
        record = umbrella_run(wildcards)
        if record.layout != "SE":
            raise WorkflowError("umbrella_stage_single requires SE; use umbrella_stage_reads")
        target = Path(output[0])
        target.parent.mkdir(parents=True, exist_ok=True)
        if record.source_type == "fastq":
            stage_local_fastq(record.source_1, target)
        else:
            with tempfile.TemporaryDirectory(dir=str(target.parent)) as work:
                work = Path(work)
                if record.source_type == "sra":
                    subprocess.run(["fasterq-dump", "--split-files", "-O", str(work),
                                    record.source_1], check=True)
                    source = work / (record.source_1 + ".fastq")
                else:
                    source = work / "R1.fastq"
                    with source.open("wb") as stream:
                        subprocess.run(["samtools", "fastq", "-F", "0x900",
                                        record.source_1], stdout=stream, check=True)
                if not source.is_file():
                    raise WorkflowError("missing single-end FASTQ from {}".format(record.run_id))
                compress_fastq(source, target)


def umbrella_validation_inputs(wildcards):
    record = umbrella_run(wildcards)
    return UMBRELLA_MANIFEST.native_reads(record.run_id)


rule umbrella_validate_reads:
    input:
        umbrella_validation_inputs
    output:
        "umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/reads.valid"
    run:
        record = umbrella_run(wildcards)
        subprocess.run(["python3", "src/validate_fastq_pairs.py", *[str(path) for path in input],
                        "--marker", str(output[0]), "--run-id", record.run_id], check=True)


rule umbrella_legacy_fastq:
    input:
        reads=umbrella_validation_inputs,
        valid="umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/reads.valid"
    output:
        temp("umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}/legacy.fastq.gz")
    run:
        record = umbrella_run(wildcards)
        if record.layout != "PE":
            raise WorkflowError("single-end runs use native R1 as legacy input")
        subprocess.run(["bash", "src/concat_mates.sh", str(input.reads[0]),
                        str(input.reads[1]), str(output[0])], check=True)


rule umbrella_intake:
    input:
        [run.work_dir + "/reads.valid" for run in UMBRELLA_MANIFEST.included_runs()]

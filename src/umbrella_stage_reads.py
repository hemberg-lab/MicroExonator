"""Stage one umbrella run's native FASTQ reads under a rule-local Conda env."""

import argparse
import bz2
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from umbrella_manifest import load_umbrella_manifest


def compress_fastq(source, destination):
    with open(source, "rb") as incoming, gzip.open(destination, "wb") as outgoing:
        shutil.copyfileobj(incoming, outgoing)


def stage_local_fastq(source, destination):
    """Link gzip inputs; make a workflow-owned gzip for other FASTQs."""
    source = Path(source)
    destination = Path(destination)
    if source.suffix == ".gz":
        destination.symlink_to(source)
        return
    opener = bz2.open if source.suffix == ".bz2" else open
    with opener(source, "rb") as incoming, gzip.open(destination, "wb") as outgoing:
        shutil.copyfileobj(incoming, outgoing)


def stage_run(manifest_path, run_id, r1, r2=None):
    manifest = load_umbrella_manifest(manifest_path)
    record = manifest.by_run[run_id]
    if not record.include:
        raise ValueError("umbrella run_id is excluded: {}".format(run_id))
    if (record.layout == "PE") != (r2 is not None):
        raise ValueError("output mate count does not match {} layout".format(record.layout))

    targets = [Path(r1)] + ([Path(r2)] if r2 is not None else [])
    targets[0].parent.mkdir(parents=True, exist_ok=True)
    if record.source_type == "fastq":
        sources = [record.source_1] + ([record.source_2] if r2 is not None else [])
        for source, target in zip(sources, targets):
            stage_local_fastq(source, target)
        return

    with tempfile.TemporaryDirectory(dir=str(targets[0].parent)) as work_name:
        work = Path(work_name)
        if record.source_type == "sra":
            subprocess.run(["fasterq-dump", "--split-files", "-O", str(work),
                            record.source_1], check=True)
            suffixes = ("_1.fastq", "_2.fastq") if r2 is not None else (".fastq",)
            sources = [work / (record.source_1 + suffix) for suffix in suffixes]
        elif r2 is not None:
            sources = [work / "R1.fastq", work / "R2.fastq"]
            sorted_bam = work / "name_sorted.bam"
            subprocess.run(["samtools", "sort", "-n", "-o", str(sorted_bam),
                            record.source_1], check=True)
            subprocess.run(["samtools", "fastq", "-1", str(sources[0]), "-2",
                            str(sources[1]), "-0", os.devnull, "-s", os.devnull,
                            "-n", "-F", "0x900", str(sorted_bam)], check=True)
        else:
            sources = [work / "R1.fastq"]
            with sources[0].open("wb") as stream:
                subprocess.run(["samtools", "fastq", "-F", "0x900",
                                record.source_1], stdout=stream, check=True)
        for source, target in zip(sources, targets):
            if not source.is_file():
                raise FileNotFoundError("missing FASTQ from {}: {}".format(record.run_id, source))
            compress_fastq(source, target)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2")
    args = parser.parse_args(argv)
    stage_run(args.manifest, args.run_id, args.r1, args.r2)


if __name__ == "__main__":
    main()

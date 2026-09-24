#!/usr/bin/env python3
"""Stage and run one reference-dependent mm10 MicroExonator smoke test."""

import argparse
import csv
import gzip
import json
import os
import pathlib
import shutil
import subprocess
import sys


PINNED_SNAKEMAKE = "7.32.4"


def read_simple_yaml(path):
    """Read the JSON-formatted YAML emitted by this runner."""
    with open(path) as handle:
        return json.load(handle)


def _absolute_local_manifest(source, destination):
    with open(source) as input_handle, open(destination, "w", newline="") as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.DictWriter(
            output_handle,
            fieldnames=("sample", "path"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            path = pathlib.Path(row["path"])
            if not path.is_absolute():
                path = pathlib.Path(source).parent / path
            writer.writerow({"sample": row["sample"], "path": str(path.resolve())})


def prepare_workdir(
    layout,
    genome_fasta,
    annotation_bed,
    workdir,
    fixture_root=None,
    repository=None,
):
    """Create the isolated run directory and return its concrete config path."""
    if layout not in ("single_end", "paired_end"):
        raise ValueError("layout must be single_end or paired_end")
    fixture_root = pathlib.Path(fixture_root or pathlib.Path(__file__).parent).resolve()
    repository = pathlib.Path(repository or fixture_root.parents[2]).resolve()
    genome_fasta = pathlib.Path(genome_fasta).resolve()
    annotation_bed = pathlib.Path(annotation_bed).resolve()
    workdir = pathlib.Path(workdir).resolve()
    if not genome_fasta.is_file():
        raise FileNotFoundError(genome_fasta)
    if not annotation_bed.is_file():
        raise FileNotFoundError(annotation_bed)
    if workdir.exists() and any(workdir.iterdir()):
        raise RuntimeError(
            "work directory is not empty; choose a new path so prior results are preserved: {}".format(
                workdir
            )
        )
    workdir.mkdir(parents=True, exist_ok=True)
    reference_dir = workdir / "reference"
    reference_dir.mkdir()
    staged_annotation = reference_dir / "GENCODE_VM25_comprehensive.bed"
    if annotation_bed.suffix == ".gz":
        with gzip.open(annotation_bed, "rb") as source, open(staged_annotation, "wb") as destination:
            shutil.copyfileobj(source, destination)
    else:
        shutil.copyfile(annotation_bed, staged_annotation)

    layout_root = fixture_root / layout
    _absolute_local_manifest(layout_root / "local_samples.tsv", workdir / "local_samples.tsv")
    with open(workdir / "local_samples.tsv") as handle:
        local_samples = list(csv.DictReader(handle, delimiter="\t"))
    download_dir = workdir / "download"
    download_dir.mkdir()
    for row in local_samples:
        script_path = download_dir / "{}.download.sh".format(row["sample"])
        with open(script_path, "w") as handle:
            handle.write("#!/bin/bash\n")
            handle.write("# FASTQ is staged by run_smoke_test.py\n")
    fastq_dir = workdir / "FASTQ"
    fastq_dir.mkdir()
    for row in local_samples:
        (fastq_dir / "{}.fastq.gz".format(row["sample"])).symlink_to(row["path"])
    staged_bulk = workdir / "bulk_samples.tsv"
    shutil.copyfile(layout_root / "bulk_samples.tsv", staged_bulk)
    paired_samples = "F"
    if layout == "paired_end":
        staged_pairs = workdir / "paired_samples.tsv"
        shutil.copyfile(layout_root / "paired_samples.tsv", staged_pairs)
        paired_samples = str(staged_pairs)

    source_link = workdir / "src"
    if not source_link.exists():
        source_link.symlink_to(repository / "src", target_is_directory=True)
    data_dir = workdir / "data"
    data_dir.mkdir(exist_ok=True)
    (data_dir / "Genome").symlink_to(genome_fasta)
    (workdir / "logs").mkdir(exist_ok=True)

    config = {
        "Genome_fasta": str(genome_fasta),
        "Gene_anontation_bed12": str(staged_annotation),
        "Gene_anontation_GTF": "NA",
        "GT_AG_U2_5": str(repository / "PWM" / "Mouse" / "mm10_GT_AG_U2_5.good.matrix"),
        "GT_AG_U2_3": str(repository / "PWM" / "Mouse" / "mm10_GT_AG_U2_3.good.matrix"),
        "ME_DB": str(fixture_root / "empty_microexons.bed"),
        "conservation_bigwig": "NA",
        "ME_len": "30",
        "max_read_len": "100",
        "min_reads_PSI": "5",
        "min_number_files_detected": "2",
        "min_detected_samples": "2",
        "filter_method": "robustness",
        "bulk_samples": str(staged_bulk),
        "paired_samples": paired_samples,
        "Single_Cell": "F",
        "Keep_fastq_gz": "T",
        "Optimize_hard_drive": "F",
        "working_directory": str(workdir) + os.sep,
    }
    config_path = workdir / "config.yaml"
    with open(config_path, "w") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return config_path


def _tool_version(command):
    result = subprocess.run(
        [str(command), "--version"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    return result.stdout.strip().splitlines()[0]


def _resolve_samtools(snakemake):
    configured = os.environ.get("MICROEXONATOR_SAMTOOLS")
    if configured:
        return configured
    discovered = shutil.which("samtools")
    if discovered:
        return discovered
    snakemake_path = pathlib.Path(snakemake).resolve()
    sibling = snakemake_path.parents[2] / "biotools" / "bin" / "samtools"
    if sibling.is_file():
        return str(sibling)
    raise RuntimeError(
        "samtools was not found; put it on PATH or set MICROEXONATOR_SAMTOOLS"
    )


def _build_temporary_fai(samtools, genome_fasta, workdir):
    index_path = pathlib.Path(workdir) / "reference" / "mm10.fa.fai"
    subprocess.run(
        [str(samtools), "faidx", "--fai-idx", str(index_path), str(genome_fasta)],
        check=True,
    )
    return index_path


def _execution_environment(snakemake):
    """Put the Conda base bin first so rule-environment activation stays first."""
    environment = os.environ.copy()
    snakemake_path = pathlib.Path(snakemake).resolve()
    if len(snakemake_path.parents) >= 4 and snakemake_path.parents[2].name == "envs":
        conda_bin = snakemake_path.parents[3] / "bin"
        path_entries = [
            entry
            for entry in environment.get("PATH", "").split(os.pathsep)
            if entry and pathlib.Path(entry) != conda_bin
        ]
        environment["PATH"] = os.pathsep.join([str(conda_bin)] + path_entries)
    return environment


def _run_and_log(command, workdir, environment=None):
    log_path = pathlib.Path(workdir) / "logs" / "snakemake.log"
    with open(log_path, "w") as log_handle:
        process = subprocess.Popen(
            command,
            cwd=workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            env=environment,
        )
        for line in process.stdout:
            sys.stdout.write(line)
            log_handle.write(line)
        return_code = process.wait()
    if return_code:
        raise RuntimeError(
            "Snakemake exited with status {}; failed work directory and log were preserved at {}".format(
                return_code, workdir
            )
        )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout", required=True, choices=("single_end", "paired_end"))
    parser.add_argument("--genome-fasta", required=True, type=pathlib.Path)
    parser.add_argument("--annotation-bed", required=True, type=pathlib.Path)
    parser.add_argument("--snakemake", required=True, type=pathlib.Path)
    parser.add_argument("--workdir", required=True, type=pathlib.Path)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument(
        "--fixture-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent,
        help="fixture laid out like this directory (default: the committed mm10 fixture)",
    )
    args = parser.parse_args(argv)

    fixture_root = args.fixture_dir.resolve()
    repository = pathlib.Path(__file__).resolve().parents[3]
    observed_version = _tool_version(args.snakemake)
    if observed_version != PINNED_SNAKEMAKE:
        raise RuntimeError(
            "this fixture requires Snakemake {}; observed {}".format(
                PINNED_SNAKEMAKE, observed_version
            )
        )
    config_path = prepare_workdir(
        args.layout,
        args.genome_fasta,
        args.annotation_bed,
        args.workdir,
        fixture_root,
        repository,
    )
    samtools = _resolve_samtools(args.snakemake)
    temporary_fai = _build_temporary_fai(samtools, args.genome_fasta.resolve(), args.workdir.resolve())
    print("Temporary FASTA index: {}".format(temporary_fai))
    command = [
        str(args.snakemake.resolve()),
        "--snakefile",
        str(repository / "MicroExonator.smk"),
        "--directory",
        str(args.workdir.resolve()),
        "--configfile",
        str(config_path),
        "--cores",
        str(args.cores),
        "--use-conda",
        "--conda-frontend",
        "mamba",
        "--conda-prefix",
        str(repository / ".snakemake" / "conda"),
        "--rerun-incomplete",
        "--printshellcmds",
        "quant",
    ]
    _run_and_log(
        command,
        args.workdir.resolve(),
        environment=_execution_environment(args.snakemake),
    )

    validator = fixture_root / "validate_smoke.py"
    subprocess.run(
        [
            sys.executable,
            str(validator),
            "--layout",
            args.layout,
            "--fixture-dir",
            str(fixture_root),
            "--run-dir",
            str(args.workdir.resolve()),
        ],
        check=True,
    )
    print("Smoke test passed; work directory retained at {}".format(args.workdir.resolve()))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(1)

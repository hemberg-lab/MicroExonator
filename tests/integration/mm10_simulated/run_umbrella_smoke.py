"""Umbrella-mode smoke run on the mm10 simulated fixture.

Six runs in two groups (A, B) and two batches: batch1 is single-end
(A1, A2, B1, B2), batch2 is paired-end (A3, B3), so both read dispatches, the
batch shards and their joins are exercised. Two comparisons: A_vs_B mixes
layouts (rMATS must refuse), A_vs_B_single_end leaves batch2 out.

The run never builds an index. --umbrella-reference is a JSON mapping for
config `umbrella_reference` and must name prebuilt full-mm10 indexes
(hisat2_index_prefix, whippet_index, salmon_index) plus the annotation files.
The fixture reads are small, so they are used whole against the full index.

  --dry-run   stage the work directory and run `snakemake -n quant_umbrella`
              (no tools needed; the reference_id may be a placeholder)
  default     full run with --use-conda, then validate_umbrella_smoke.py.
              rMATS and SUPPA2 conda packages are Linux-first: run this on
              Linux or in a Linux container.
"""

import argparse
import csv
import gzip
import importlib.util
import json
import os
import pathlib
import shutil
import subprocess
import sys


FIXTURE = pathlib.Path(__file__).resolve().parent
REPOSITORY = FIXTURE.parents[2]
PREBUILT = ("hisat2_index_prefix", "whippet_index", "salmon_index")
REQUIRED = ("genome_fasta", "annotation_gtf", "whippet_gtf", "me_db", "transcriptome_fasta",
            "decoys", "salmon_gtf", "splice_sites")


def _legacy_runner():
    spec = importlib.util.spec_from_file_location("mm10_smoke_runner", FIXTURE / "run_smoke_test.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_reference(reference):
    missing = [key for key in PREBUILT if not reference.get(key)]
    if missing:
        raise RuntimeError("the smoke run never builds an index; --umbrella-reference must name "
                           "prebuilt {}".format(", ".join(missing)))
    missing = [key for key in REQUIRED if not reference.get(key)]
    if missing:
        raise RuntimeError("--umbrella-reference lacks {}".format(", ".join(missing)))


def stage(workdir, genome_fasta, annotation_bed, reference, reference_id, fixture=FIXTURE,
          repository=REPOSITORY):
    """Write manifest, comparisons and config into workdir; return the config path."""
    check_reference(reference)
    workdir = pathlib.Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=False)
    with open(fixture / "umbrella_manifest.tsv") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    fields = list(rows[0])
    with open(workdir / "umbrella_manifest.tsv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            for key in ("source_1", "source_2"):
                if row[key]:
                    row[key] = str((fixture / row[key]).resolve())
            row["reference_id"] = reference_id
            writer.writerow(row)
    shutil.copyfile(fixture / "umbrella_comparisons.yaml", workdir / "umbrella_comparisons.yaml")

    staged_bed = workdir / "reference" / "annotation.bed"
    staged_bed.parent.mkdir()
    opener = gzip.open if str(annotation_bed).endswith(".gz") else open
    with opener(annotation_bed, "rb") as source, open(staged_bed, "wb") as destination:
        shutil.copyfileobj(source, destination)
    (workdir / "src").symlink_to(repository / "src", target_is_directory=True)
    (workdir / "data").mkdir()
    (workdir / "data" / "Genome").symlink_to(pathlib.Path(genome_fasta).resolve())
    (workdir / "logs").mkdir()

    config = {
        "Genome_fasta": str(pathlib.Path(genome_fasta).resolve()),
        "Gene_anontation_bed12": str(staged_bed),
        "Gene_anontation_GTF": reference["annotation_gtf"],
        "GT_AG_U2_5": str(repository / "PWM" / "Mouse" / "mm10_GT_AG_U2_5.good.matrix"),
        "GT_AG_U2_3": str(repository / "PWM" / "Mouse" / "mm10_GT_AG_U2_3.good.matrix"),
        "ME_DB": str(fixture / "empty_microexons.bed"),
        "conservation_bigwig": "NA",
        "ME_len": "30", "max_read_len": "100", "min_reads_PSI": "5",
        "min_number_files_detected": "2", "min_detected_samples": "2",
        "filter_method": "robustness", "Single_Cell": "F", "Keep_fastq_gz": "T",
        "Optimize_hard_drive": "F", "working_directory": str(workdir) + os.sep,
        "umbrella_manifest": str(workdir / "umbrella_manifest.tsv"),
        "umbrella_comparisons": str(workdir / "umbrella_comparisons.yaml"),
        "umbrella_reference": reference,
    }
    config_path = workdir / "config.yaml"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    return config_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome-fasta", required=True, type=pathlib.Path)
    parser.add_argument("--annotation-bed", required=True, type=pathlib.Path)
    parser.add_argument("--umbrella-reference", required=True, type=pathlib.Path,
                        help="JSON for config umbrella_reference, with prebuilt indexes")
    parser.add_argument("--reference-id",
                        help="content-derived reference_id (required unless --dry-run)")
    parser.add_argument("--snakemake", required=True, type=pathlib.Path)
    parser.add_argument("--workdir", required=True, type=pathlib.Path)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--genome-index-dir", type=pathlib.Path,
                        help="reuse the legacy data/Genome Bowtie + HISAT2 index")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    with open(args.umbrella_reference) as stream:
        reference = json.load(stream)
    if not args.reference_id and not args.dry_run:
        parser.error("--reference-id is required for a real run; compute it with "
                     "src/shard_guard.py reference_manifest or take it from a finished run")
    config_path = stage(args.workdir, args.genome_fasta, args.annotation_bed, reference,
                        args.reference_id or "dryrun")
    runner = _legacy_runner()
    if args.genome_index_dir:
        runner.link_genome_index(args.genome_index_dir, args.workdir)
    command = [str(args.snakemake.resolve()), "--snakefile", str(REPOSITORY / "MicroExonator.smk"),
               "--directory", str(args.workdir.resolve()), "--configfile", str(config_path),
               "--cores", str(args.cores), "--printshellcmds", "quant_umbrella"]
    if args.dry_run:
        command.append("-n")
    else:
        command += ["--use-conda", "--conda-prefix", str(REPOSITORY / ".snakemake" / "conda"),
                    "--rerun-incomplete"]
    runner._run_and_log(command, args.workdir.resolve(),
                        environment=runner._execution_environment(args.snakemake))
    if not args.dry_run:
        subprocess.run([sys.executable, str(FIXTURE / "validate_umbrella_smoke.py"),
                        "--workdir", str(args.workdir.resolve())], check=True)


if __name__ == "__main__":
    main()

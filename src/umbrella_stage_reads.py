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


# Every staged FASTQ is cleaned in one streaming pass while it is compressed:
#
# - The separator line becomes a bare "+". fasterq-dump repeats the read name
#   there ("+SRR1.1 ..."), and Whippet 1.6 then fails on any read that
#   contains N ("Cannot encode 78 to DNAAlphabet{2}"; its author: Whippet only
#   accepts the standard four-line FASTQ). Tested on the hg38 pilot: "+name"
#   with N reads crashes in single- and paired-end mode; a bare "+" with the
#   same N reads runs. Legacy fastq-dump --defline-qual '+' wrote it bare too.
# - Header lines keep printable ASCII only: any other byte (non-ASCII letters
#   in submitter read names, tabs, control characters) becomes "_". Such
#   names have crashed Bowtie and other tools downstream (the reason for the
#   legacy validate_fastq rule), and a tab breaks SAM files.
#
# awk runs in the C locale, so it works on bytes whatever the encoding.
CLEAN_FASTQ = ["awk", 'NR % 4 == 1 { gsub(/[^ -~]/, "_") } NR % 4 == 3 { print "+"; next } { print }']
CLEAN_ENV = dict(os.environ, LC_ALL="C")
# local .gz inputs are linked unless the first records need cleaning
SAMPLE_RECORDS = 10000


def _gzip_stream(incoming, destination, threads):
    """Write a binary FASTQ stream to destination as gzip, cleaned (CLEAN_FASTQ).

    pigz compresses when available, gzip otherwise.
    """
    pigz = shutil.which("pigz")
    compress = [pigz, "-p", str(max(1, threads)), "-c"] if pigz else ["gzip", "-c"]
    with open(destination, "wb") as outgoing:
        rewrite = subprocess.Popen(CLEAN_FASTQ, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   env=CLEAN_ENV)
        packer = subprocess.Popen(compress, stdin=rewrite.stdout, stdout=outgoing)
        rewrite.stdout.close()   # packer owns the read end now
        try:
            shutil.copyfileobj(incoming, rewrite.stdin, 1024 * 1024)
        finally:
            rewrite.stdin.close()
        for process, name in ((rewrite, "awk"), (packer, compress[0])):
            if process.wait() != 0:
                raise subprocess.CalledProcessError(process.returncode, name)


def compress_fastq(source, destination, threads=1):
    """gzip an uncompressed FASTQ, then delete it so the disk peak stays low."""
    with open(source, "rb") as incoming:
        _gzip_stream(incoming, destination, threads)
    Path(source).unlink()


def _is_clean(path, records=SAMPLE_RECORDS):
    """True when the first `records` records need no cleaning: bare "+" lines
    and printable-ASCII headers. Read names follow one pattern per run, so a
    sample is enough to decide; a file that fails is cleaned in full."""
    with gzip.open(path, "rb") as stream:
        for number in range(records * 4):
            line = stream.readline()
            if not line:
                break
            line = line.rstrip(b"\r\n")
            if number % 4 == 0 and any(byte < 0x20 or byte > 0x7E for byte in line):
                return False
            if number % 4 == 2 and line != b"+":
                return False
    return True


def stage_local_fastq(source, destination, threads=1):
    """Link gzip inputs that need no cleaning; clean and rewrite everything else."""
    source = Path(source)
    destination = Path(destination)
    if source.suffix == ".gz" and _is_clean(source):
        destination.symlink_to(source)
        return
    opener = {".gz": gzip.open, ".bz2": bz2.open}.get(source.suffix, open)
    with opener(source, "rb") as incoming:
        _gzip_stream(incoming, destination, threads)


def _sra_reads(work, accession, paired):
    """FASTQ files fasterq-dump --split-3 wrote for one accession.

    --split-3 keeps _1/_2 in step (spots with one mate go to <acc>.fastq,
    which a paired run ignores); a single-end run is <acc>.fastq.
    """
    mates = [work / (accession + "_1.fastq"), work / (accession + "_2.fastq")]
    single = work / (accession + ".fastq")
    if paired:
        missing = [str(path.name) for path in mates if not path.is_file()]
        if missing:
            raise FileNotFoundError("{} gave no {}: is the run paired-end in SRA?".format(
                accession, ", ".join(missing)))
        if single.is_file():
            single.unlink()   # unpaired spots; not used
        return mates
    if single.is_file():
        return [single]
    if mates[0].is_file() and not mates[1].is_file():
        return [mates[0]]
    raise FileNotFoundError("{} gave no single-end FASTQ: is the run single-end in SRA?".format(accession))


def stage_sra(accession, work, paired, threads=1):
    """Download with prefetch, dump with fasterq-dump; returns the FASTQ paths."""
    subprocess.run(["prefetch", "--max-size", "u", "-O", str(work), accession], check=True)
    downloaded = sorted((work / accession).glob(accession + ".sra*")) if (work / accession).is_dir() else []
    source = str(downloaded[0]) if downloaded else accession
    subprocess.run(["fasterq-dump", "--split-3", "--threads", str(max(1, threads)),
                    "--temp", str(work), "-O", str(work), source], check=True)
    if (work / accession).is_dir():
        shutil.rmtree(work / accession)   # the .sra is not needed once dumped
    return _sra_reads(work, accession, paired)


def processed_outputs(record):
    """Kept outputs that exist once a run is fully processed (paths from the workflow root)."""
    shard = "{}/{}/{}/{}".format(record.reference_id, record.project_id, record.group, record.batch_id)
    return [record.work_dir + "/reads.valid", record.work_dir + "/whippet/quant.psi.gz",
            "qc/{}.qc.tsv.gz".format(shard), "splicing/{}.microexonator.tsv.gz".format(shard),
            "splicing/{}.whippet.tsv.gz".format(shard)]


def stage_run(manifest_path, run_id, r1, r2=None, threads=1, tmpdir=None, allow_restage=False):
    manifest = load_umbrella_manifest(manifest_path)
    record = manifest.by_run[run_id]
    if not record.include:
        raise ValueError("umbrella run_id is excluded: {}".format(run_id))
    if not allow_restage and all(Path(path).exists() for path in processed_outputs(record)):
        # Every umbrella output of this run is kept, so nothing should need its
        # reads again: something upstream changed. Refuse before downloading.
        raise RuntimeError(
            "run {} is already fully processed, but a job asked for its reads again, which "
            "would {} it. Usually a rule, parameter or input upstream changed; "
            "`snakemake -n -r quant_umbrella` shows why. To reprocess on purpose, set "
            "umbrella_allow_restage: true".format(
                run_id, "download" if record.source_type == "sra" else "restage"))
    if (record.layout == "PE") != (r2 is not None):
        raise ValueError("output mate count does not match {} layout".format(record.layout))

    targets = [Path(r1)] + ([Path(r2)] if r2 is not None else [])
    targets[0].parent.mkdir(parents=True, exist_ok=True)
    if record.source_type == "fastq":
        sources = [record.source_1] + ([record.source_2] if r2 is not None else [])
        for source, target in zip(sources, targets):
            stage_local_fastq(source, target, threads)
        return

    # Uncompressed reads live here only until each is gzipped; point tmpdir at
    # node-local scratch to keep them off the shared file system.
    # tmpdir may name an environment variable, e.g. "$TMPDIR" for the per-job
    # node-local disk PBS provides; it is expanded here, on the node. Unset or
    # unexpanded, the run's work directory is used.
    tmpdir = os.path.expandvars(tmpdir) if tmpdir else None
    scratch = Path(tmpdir) if tmpdir and "$" not in tmpdir else targets[0].parent
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(scratch), prefix="stage_{}_".format(run_id)) as work_name:
        work = Path(work_name)
        if record.source_type == "sra":
            sources = stage_sra(record.source_1, work, r2 is not None, threads)
        elif r2 is not None:
            sources = [work / "R1.fastq", work / "R2.fastq"]
            sorted_bam = work / "name_sorted.bam"
            subprocess.run(["samtools", "sort", "-n", "-@", str(max(1, threads)),
                            "-T", str(work / "sort"), "-o", str(sorted_bam),
                            record.source_1], check=True)
            subprocess.run(["samtools", "fastq", "-@", str(max(1, threads)),
                            "-1", str(sources[0]), "-2", str(sources[1]),
                            "-0", os.devnull, "-s", os.devnull,
                            "-n", "-F", "0x900", str(sorted_bam)], check=True)
            sorted_bam.unlink()
        else:
            sources = [work / "R1.fastq"]
            with sources[0].open("wb") as stream:
                subprocess.run(["samtools", "fastq", "-@", str(max(1, threads)), "-F", "0x900",
                                record.source_1], stdout=stream, check=True)
        for source, target in zip(sources, targets):
            if not source.is_file():
                raise FileNotFoundError("missing FASTQ from {}: {}".format(record.run_id, source))
            compress_fastq(source, target, threads)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--tmpdir", help="scratch for uncompressed reads (default: next to the outputs)")
    parser.add_argument("--allow-restage", action="store_true",
                        help="stage a run even when all its outputs are kept already")
    args = parser.parse_args(argv)
    stage_run(args.manifest, args.run_id, args.r1, args.r2, args.threads, args.tmpdir,
              args.allow_restage)


if __name__ == "__main__":
    main()

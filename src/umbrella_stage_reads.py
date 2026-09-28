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


def _gzip_stream(incoming, destination, threads):
    """Write a binary stream to destination as gzip: pigz when available."""
    pigz = shutil.which("pigz")
    if pigz is None:
        with gzip.open(destination, "wb") as outgoing:
            shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
        return
    with open(destination, "wb") as outgoing:
        process = subprocess.Popen([pigz, "-p", str(max(1, threads)), "-c"],
                                   stdin=subprocess.PIPE, stdout=outgoing)
        try:
            shutil.copyfileobj(incoming, process.stdin, 1024 * 1024)
        finally:
            process.stdin.close()
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, "pigz")


def compress_fastq(source, destination, threads=1):
    """gzip an uncompressed FASTQ, then delete it so the disk peak stays low."""
    with open(source, "rb") as incoming:
        _gzip_stream(incoming, destination, threads)
    Path(source).unlink()


def stage_local_fastq(source, destination, threads=1):
    """Link gzip inputs; make a workflow-owned gzip for other FASTQs."""
    source = Path(source)
    destination = Path(destination)
    if source.suffix == ".gz":
        destination.symlink_to(source)
        return
    opener = bz2.open if source.suffix == ".bz2" else open
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


def stage_run(manifest_path, run_id, r1, r2=None, threads=1, tmpdir=None):
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
    args = parser.parse_args(argv)
    stage_run(args.manifest, args.run_id, args.r1, args.r2, args.threads, args.tmpdir)


if __name__ == "__main__":
    main()

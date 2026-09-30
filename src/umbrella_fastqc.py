"""Lightweight per-run QC kept next to the umbrella QC shards, for MultiQC.

fastqc   FastQC on the first N reads of each staged file (default 2M; 0 for
         all). A file's start shows base quality, GC, adapters, duplication
         and overrepresented sequences, but not tile-to-tile differences.
         Only fastqc_data.txt (tens of KB) is kept, in <out>/<run>_R<n>/,
         which MultiQC reads.
keep     Copies the small summaries MultiQC reads before their temporary
         per-run sources are deleted: the HISAT2 summary, the featureCounts
         summary (its BAM column renamed to the run ID, so MultiQC names the
         sample), Salmon's aux_info/meta_info.json and lib_format_counts.json
         (under salmon/<run>/, the layout MultiQC expects), and Whippet's .map.
"""

import argparse
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path


def subsample(source, destination, reads):
    """Write the first `reads` FASTQ records of a gzipped file (all when 0)."""
    limit = 4 * reads if reads else None
    with gzip.open(source, "rt") as stream, open(destination, "w") as out:
        for number, line in enumerate(stream):
            if limit is not None and number >= limit:
                break
            out.write(line)


def fastqc(reads, run_id, out, reads_per_file=2000000, threads=1, scratch=None):
    out = Path(out)
    with tempfile.TemporaryDirectory(prefix="fastqc-", dir=scratch) as directory:
        work = Path(directory)
        names = []
        for mate, source in enumerate(reads, start=1):
            name = "{}_R{}".format(run_id, mate)
            subsample(source, work / (name + ".fastq"), reads_per_file)
            names.append(name)
        subprocess.run(["fastqc", "--extract", "--quiet", "--threads", str(threads),
                        "-o", str(work)] + [str(work / (name + ".fastq")) for name in names],
                       check=True)
        for name in names:
            data = work / (name + "_fastqc") / "fastqc_data.txt"
            if not data.is_file():
                raise RuntimeError("FastQC wrote no fastqc_data.txt for {}".format(name))
            (out / name).mkdir(parents=True, exist_ok=True)
            shutil.copyfile(data, out / name / "fastqc_data.txt")


def keep(run_id, hisat2, featurecounts, salmon, whippet_map, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(hisat2, out / "{}.hisat2.summary.txt".format(run_id))
    with open(featurecounts) as stream:
        lines = stream.read().splitlines()
    header = lines[0].split("\t")
    lines[0] = "\t".join(header[:1] + [run_id])
    (out / "{}.featurecounts.txt.summary".format(run_id)).write_text("\n".join(lines) + "\n")
    salmon, target = Path(salmon), out / "salmon" / run_id
    (target / "aux_info").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(salmon / "aux_info" / "meta_info.json", target / "aux_info" / "meta_info.json")
    shutil.copyfile(salmon / "lib_format_counts.json", target / "lib_format_counts.json")
    shutil.copyfile(whippet_map, out / "{}.whippet.map.gz".format(run_id))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    run_fastqc = commands.add_parser("fastqc")
    run_fastqc.add_argument("--run-id", required=True)
    run_fastqc.add_argument("--reads", nargs="+", required=True)
    run_fastqc.add_argument("--out", required=True)
    run_fastqc.add_argument("--reads-per-file", type=int, default=2000000)
    run_fastqc.add_argument("--threads", type=int, default=1)
    run_fastqc.add_argument("--scratch")
    run_keep = commands.add_parser("keep")
    for name in ("run-id", "hisat2", "featurecounts", "salmon", "whippet-map", "out"):
        run_keep.add_argument("--" + name, required=True)
    args = parser.parse_args(argv)
    if args.command == "fastqc":
        fastqc(args.reads, args.run_id, args.out, args.reads_per_file, args.threads, args.scratch)
    else:
        keep(args.run_id, args.hisat2, args.featurecounts, args.salmon, args.whippet_map, args.out)


if __name__ == "__main__":
    main()

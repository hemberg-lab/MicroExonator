"""Keep circular output in every new Whippet quantification."""

import shlex


def quant_flags(flags=""):
    flags = str(flags or "").strip()
    return flags if "--circ" in shlex.split(flags) else (flags + " --circ").strip()


def read_paths(record):
    source = record.metadata.get("whippet_read_source", "shared") or "shared"
    if source not in ("shared", "backfill"):
        raise ValueError("whippet_read_source must be shared or backfill")
    directory = record.work_dir
    if source == "backfill":
        directory = directory.replace("umbrella/work/", "umbrella/whippet_backfill/", 1)
    return [directory + "/R1.fastq.gz"] + ([directory + "/R2.fastq.gz"] if record.layout == "PE" else [])

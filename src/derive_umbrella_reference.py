"""Small deterministic reference files derived from the supplied genome/GTF."""

import argparse
import gzip
from pathlib import Path

if __package__:
    from src.summarize_junctions import _gtf_introns
else:
    from summarize_junctions import _gtf_introns


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def decoy_names(genome_fasta):
    names = []
    seen = set()
    with open_text(genome_fasta) as stream:
        for line in stream:
            if not line.startswith(">"):
                continue
            name = line[1:].split()[0]
            if not name or name in seen:
                raise ValueError("empty or repeated genome FASTA sequence name: {}".format(name))
            seen.add(name)
            names.append(name)
    if not names:
        raise ValueError("genome FASTA has no sequence headers")
    return "".join(name + "\n" for name in names)


def splice_sites(annotation_gtf):
    introns = _gtf_introns(annotation_gtf)
    if not introns:
        raise ValueError("annotation GTF contains no transcript introns")
    return "".join("{}\t{}\t{}\t{}\n".format(chrom, start - 2, end, strand)
                   for chrom, start, end, strand in sorted(introns))


def validate_transcripts(transcript_fasta, tx2gene):
    expected = set()
    with open_text(tx2gene) as stream:
        next(stream)
        for line in stream:
            expected.add(line.split("\t", 1)[0])
    observed = set()
    with open_text(transcript_fasta) as stream:
        for line in stream:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in observed:
                    raise ValueError("duplicate transcript FASTA identifier: {}".format(name))
                observed.add(name)
    if not observed or observed != expected:
        raise ValueError("transcript FASTA / tx2gene mismatch: {} missing, {} extra".format(
            len(expected - observed), len(observed - expected)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("decoys", "splice-sites"):
        sub = commands.add_parser(command)
        sub.add_argument("--input", required=True)
        sub.add_argument("--output", required=True)
    validate = commands.add_parser("validate-transcripts")
    validate.add_argument("--fasta", required=True)
    validate.add_argument("--tx2gene", required=True)
    args = parser.parse_args(argv)
    if args.command == "validate-transcripts":
        validate_transcripts(args.fasta, args.tx2gene)
    else:
        content = decoy_names(args.input) if args.command == "decoys" else splice_sites(args.input)
        Path(args.output).write_text(content)


if __name__ == "__main__":
    main()

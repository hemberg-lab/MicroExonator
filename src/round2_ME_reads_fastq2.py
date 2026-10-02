"""Write the raw FASTQ records of the reads that aligned to microexon tags in Round2.

Usage: round2_ME_reads_fastq2.py <sample>.sam.pre_processed <sample>.fastq.gz
Writes <sample>.sam.pre_processed.fastq.

The original (untrimmed) reads are taken from the FASTQ, because Round2 aligns
reads trimmed to max_read_len. Each selected record is written as
"@<id>\\n<seq>\\n+\\n<qual>\\n", where <id> is the header up to the first
whitespace, which is what the earlier Biopython SeqIO version wrote.

The FASTQ is scanned as plain four-line records (no SeqRecord objects) and
decompressed by a separate pigz/gzip process, so this step is I/O bound rather
than parser bound.
"""

import csv
import shutil
import subprocess
import sys

csv.field_size_limit(100000000)


def ME_read_ids(alingment_pre_processed_round2):

    ME_reads = set([])

    with open(alingment_pre_processed_round2) as f:
        for row in csv.reader(f, delimiter = '\t'):
            if len(row)>=7:
                ME_reads.add(row[0].encode())

    return ME_reads


def decompress(row_fastq):

    tool = shutil.which("pigz") or "gzip"
    return subprocess.Popen([tool, "-dcf", row_fastq], stdout=subprocess.PIPE, bufsize=1 << 20)


def main(alingment_pre_processed_round2, row_fastq):

    ME_reads = ME_read_ids(alingment_pre_processed_round2)

    proc = decompress(row_fastq)
    readline = proc.stdout.readline

    with open(alingment_pre_processed_round2 + ".fastq", 'wb', buffering = 1 << 20) as fastq_out:

        n = 0
        while True:

            header = readline()

            if not header:
                break

            if not header.strip():  # tolerate blank lines between records, as SeqIO did
                continue

            seq = readline()
            plus = readline()
            qual = readline()
            n += 1

            if header[:1] != b"@" or plus[:1] != b"+" or not qual:
                sys.exit("{}: record {} is not a four-line FASTQ record".format(row_fastq, n))

            fields = header[1:].split(None, 1)
            read_id = fields[0] if fields else b""

            if read_id in ME_reads:

                seq = seq.rstrip()
                qual = qual.rstrip()

                if len(seq) != len(qual):
                    sys.exit("{}: read {} has sequence and quality of different lengths".format(
                        row_fastq, read_id.decode()))

                fastq_out.write(b"@" + read_id + b"\n" + seq + b"\n+\n" + qual + b"\n")

    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit("{}: decompression failed".format(row_fastq))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

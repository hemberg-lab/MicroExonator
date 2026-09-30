#!/bin/bash
# Concatenate the two mate files of a paired-end run into one FASTQ, appending
# _1 and _2 to the read names (not /1 and /2, which HISAT2 and Bowtie 2 strip
# from SAM read names). Without the suffixes both mates of a fragment
# share one name, and every step that works by read name treats them as one
# read: the Round2 genome blacklist then removes a junction read whenever its
# mate aligns unspliced to the genome.
#
# Usage: concat_mates.sh mate_1.fastq.gz mate_2.fastq.gz out.fastq.gz [threads]
# Compresses with pigz when it is on PATH, otherwise gzip; the output is the same.
#
# The separator line is written as a bare "+". fasterq-dump repeats the read
# name there ("+SRR1.1 ..."); left as it was, it would no longer match the
# renamed header, and Biopython (Round2 ME_reads) rejects such records with
# "Sequence and quality captions differ".

set -euo pipefail

mate1=$1
mate2=$2
out=$3
threads=${4:-1}

if command -v pigz > /dev/null; then
    compress=(pigz -p "$threads" -c)
else
    compress=(gzip -c)
fi

{
    gzip -dc "$mate1" | awk 'NR % 4 == 1 { $1 = $1 "_1" } NR % 4 == 3 { $0 = "+" } { print }'
    gzip -dc "$mate2" | awk 'NR % 4 == 1 { $1 = $1 "_2" } NR % 4 == 3 { $0 = "+" } { print }'
} | "${compress[@]}" > "$out"

#!/bin/bash
# Copy a local FASTQ (.fastq, .fastq.gz or .fastq.bz2) to a gzipped FASTQ,
# cleaning two things that break downstream tools:
#
# - The separator line is written as a bare "+". Whippet 1.6 fails with
#   "Cannot encode 78 to DNAAlphabet{2}" when the "+" line repeats the read
#   name (as fasterq-dump and some sequencing facilities write it) and reads
#   contain N.
# - Header bytes outside printable ASCII (non-ASCII read names, tabs, control
#   characters) become "_". Such names crash Bowtie, and a tab breaks SAM.
#
# A .gz input whose first 10,000 records need no cleaning is copied as it is;
# anything else is rewritten and cleaned in full. awk runs in the C locale, so
# it works on bytes whatever the encoding.
#
# Usage: stage_local_fastq.sh input.fastq[.gz|.bz2] output.fastq.gz [threads]
# Compresses with pigz when it is on PATH, otherwise gzip; the output is the same.

set -euo pipefail

input=$1
output=$2
threads=${3:-1}

export LC_ALL=C

clean='NR % 4 == 1 { gsub(/[^ -~]/, "_") } NR % 4 == 3 { print "+"; next } { print }'
check='NR % 4 == 1 && /[^ -~]/ { dirty = 1; exit } NR % 4 == 3 && $0 != "+" { dirty = 1; exit } NR >= 40000 { exit } END { exit dirty }'

if command -v pigz > /dev/null; then
    compress=(pigz -p "$threads" -c)
else
    compress=(gzip -c)
fi

case "$input" in
    *.gz)
        # awk stops reading early, so gzip may end on a closed pipe
        if { gzip -dc "$input" 2> /dev/null || true; } | awk "$check"; then
            cp "$input" "$output"
        else
            gzip -dc "$input" | awk "$clean" | "${compress[@]}" > "$output"
        fi
        ;;
    *.bz2)
        bzip2 -dc "$input" | awk "$clean" | "${compress[@]}" > "$output"
        ;;
    *.fastq|*.fq)
        awk "$clean" "$input" | "${compress[@]}" > "$output"
        ;;
    *)
        echo "stage_local_fastq.sh: only fastq, fastq.gz or fastq.bz2 are supported: $input" >&2
        exit 1
        ;;
esac

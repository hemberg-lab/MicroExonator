"""Summarise read lengths per sample (Report/read_lengths.tsv).

Splice-junction tags have max_read_len nt of exon on each side, and reads
longer than that are trimmed before they are aligned to the tags. This report
shows how many reads each sample had, their median and maximum length and the
share that was trimmed. Samples whose lengths look unusual for bulk RNA-seq
(very short reads, or most reads trimmed) are noted, because they are often
another library type, such as single-cell data, and should be checked.

Usage: read_length_report.py MAX_READ_LEN SAMPLE.tsv [SAMPLE.tsv ...]
"""

import csv
import os
import sys

SHORT_READS = 40


def summarise(path, max_read_len):
    hist = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hist[int(row["length"])] = hist.get(int(row["length"]), 0) + int(row["reads"])
    total = sum(hist.values())
    if total == 0:
        return total, 0, 0, 0.0, "no reads"
    running, median = 0, 0
    for length in sorted(hist):
        running += hist[length]
        if running * 2 >= total:
            median = length
            break
    trimmed = 100.0 * sum(n for length, n in hist.items() if length > max_read_len) / total
    notes = []
    if median < SHORT_READS:
        notes.append("median read shorter than {} nt".format(SHORT_READS))
    if trimmed > 50:
        notes.append("most reads longer than max_read_len and trimmed")
    return total, median, max(hist), trimmed, "; ".join(notes)


def main(max_read_len, paths):
    print("sample\treads\tmedian_length\tmax_length\tpercent_trimmed\tnote")
    for path in sorted(paths):
        sample = os.path.basename(path)[: -len(".tsv")]
        total, median, longest, trimmed, note = summarise(path, max_read_len)
        print("{}\t{}\t{}\t{}\t{:.1f}\t{}".format(sample, total, median, longest, trimmed, note))
        if note:
            print("WARNING: {}: {}".format(sample, note), file=sys.stderr)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2:])

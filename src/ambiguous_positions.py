"""List microexons whose genomic position cannot be resolved from reads.

A junction read shows the upstream exon, the microexon sequence and the
downstream exon, but not where in the intron the microexon sits. Very short
microexons (for example a 3 nt CAG) often match several positions with
canonical splice sites in the same intron. MicroExonator reports the position
with the best U2 splice-site score; this file lists every such event with its
alternative positions, so an annotated or expected coordinate that is missing
from the output can be traced to the event that carries its reads.

Usage: ambiguous_positions.py TOTAL.ME_centric.txt > ME_ambiguous_positions.txt
"""

import csv
import sys


def ambiguous_events(me_centric):
    csv.field_size_limit(1000000000)
    reported = {}
    for row in csv.reader(open(me_centric), delimiter="\t"):
        total_ME = row[11]
        if total_ME == "":
            continue
        candidates = [i.split("|") for i in total_ME.split(",")]
        if len({c[0] for c in candidates}) < 2:
            continue
        best = max(candidates, key=lambda item: float(item[1]))
        alternatives = sorted({c[0] for c in candidates} - {best[0]})
        reported[best[0]] = (row[6], alternatives)
    return reported


def main(me_centric):
    print("ME\tME_seq\tn_positions\talternative_positions")
    for microexon, (sequence, alternatives) in sorted(ambiguous_events(me_centric).items()):
        print("\t".join([microexon, sequence, str(len(alternatives) + 1), ",".join(alternatives)]))


if __name__ == "__main__":
    main(sys.argv[1])

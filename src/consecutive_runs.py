"""List microexons that sit in runs of adjacent microexons.

Two microexons are adjacent when one of them has a containing intron that
starts or ends exactly at the other (so an annotated transcript joins them
directly). Counts for adjacent microexons are corrected using the immediate
neighbour's tags (see neighbour_skip_counts.py). In runs of three or more, a
read can also reach a microexon two positions away; that evidence is not
recovered, so PSI for these events can be slightly underestimated.

Usage: consecutive_runs.py TOTAL.ME_centric.txt > ME_consecutive_runs.txt
"""

import csv
import sys

from neighbour_skip_counts import parse_intron, parse_microexon


def event_introns(me_centric):
    csv.field_size_limit(1000000000)
    introns = {}
    for row in csv.reader(open(me_centric), delimiter="\t"):
        total_SJs, total_ME = row[3], row[11]
        if total_ME == "":
            continue
        true_ME = max([i.split("|") for i in total_ME.split(",")], key=lambda item: float(item[1]))[0]
        introns.setdefault(true_ME, set()).update(total_SJs.split(","))
    return introns


def consecutive_runs(introns):
    """{microexon: (run members sorted, neighbours)} for runs of two or more."""
    by_start, by_end = {}, {}
    for microexon in introns:
        chrom, strand, start, end = parse_microexon(microexon)
        by_start.setdefault((chrom, strand, start), set()).add(microexon)
        by_end.setdefault((chrom, strand, end), set()).add(microexon)
    neighbours = {microexon: set() for microexon in introns}
    for microexon, microexon_introns in introns.items():
        for intron in microexon_introns:
            chrom, strand, start, end = parse_intron(intron)
            for other in by_end.get((chrom, strand, start), set()) | by_start.get((chrom, strand, end), set()):
                if other != microexon:
                    neighbours[microexon].add(other)
                    neighbours[other].add(microexon)
    runs = {}
    for microexon in introns:
        if microexon in runs or not neighbours[microexon]:
            continue
        members, stack = set(), [microexon]
        while stack:
            current = stack.pop()
            if current not in members:
                members.add(current)
                stack.extend(neighbours[current] - members)
        ordered = sorted(members, key=lambda m: parse_microexon(m)[2])
        for member in members:
            runs[member] = (ordered, sorted(neighbours[member]))
    return runs


def main(me_centric):
    print("ME\trun_size\trun_members\tadjacent_microexons")
    for microexon, (members, adjacent) in sorted(consecutive_runs(event_introns(me_centric)).items()):
        print("\t".join([microexon, str(len(members)), ",".join(members), ",".join(adjacent)]))


if __name__ == "__main__":
    main(sys.argv[1])

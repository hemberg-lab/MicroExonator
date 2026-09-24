"""Keep one Bowtie alignment per read, preferring alignments that cross a junction.

Tag flanks are cut from spliced transcripts, so when an exon is shorter than
the flank, a tag's flank also contains the next junction. A read from that
junction then aligns equally well inside the neighbouring tag's flank, where it
crosses no junction and is discarded. With a single reported alignment per
read (bowtie -k 1) such reads were lost at random, and more so with longer
flanks. Bowtie now reports several alignments per read (-k N); this filter
keeps, for each read, one alignment that crosses its tag's junction with at
least 8 nt on each side, choosing the fewest mismatches and breaking ties
with a seeded choice per read. Reads with no junction-crossing alignment keep
their first alignment and are discarded downstream as before.

Usage: bowtie ... -S | select_tag_alignments.py > sample.sam.raw
Alignments of one read must be consecutive, as Bowtie writes them.
"""

import random
import sys
import zlib

MIN_ANCHOR = 8
SEED = 123


def crosses_junction(row):
    """True if a SAM row crosses its tag's junction (or a microexon junction)."""
    anchors = row[2].split("|")[-1].split("_")
    start = int(row[3]) - 1
    end = start + len(row[9])
    up = int(anchors[0])
    boundaries = [up] if len(anchors) == 2 else [up, up + len(anchors[1])]
    return any(start <= b - MIN_ANCHOR and end >= b + MIN_ANCHOR for b in boundaries)


def mismatches(row):
    for field in row[11:]:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return 0


def choose(rows):
    crossing = [row for row in rows if crosses_junction(row)]
    if not crossing:
        return rows[0]
    fewest = min(mismatches(row) for row in crossing)
    best = [row for row in crossing if mismatches(row) == fewest]
    if len(best) == 1:
        return best[0]
    return random.Random(SEED + zlib.crc32(rows[0][0].encode())).choice(best)


def select(lines):
    group, name = [], None
    for line in lines:
        if line.startswith("@"):
            yield line
            continue
        row = line.rstrip("\n").split("\t")
        if row[0] != name and group:
            yield "\t".join(choose(group)) + "\n"
            group = []
        name = row[0]
        group.append(row)
    if group:
        yield "\t".join(choose(group)) + "\n"


if __name__ == "__main__":
    sys.stdout.writelines(select(sys.stdin))

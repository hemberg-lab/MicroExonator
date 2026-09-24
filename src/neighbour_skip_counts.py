"""Exclusion and inclusion evidence for consecutive microexons from their neighbours' tags.

With two or more microexons in a row (E1-m1-m2-E2), the junctions that skip
only one of them (E1|m2 skips m1, m1|E2 skips m2) are the neighbour's own
inclusion junctions. They never get an exon-exon tag, so skipping reads land
on the neighbour's inclusion tag and were only counted as inclusion for the
neighbour. Here, every read that crosses a junction of an inclusion tag with
8 nt anchors (the same rule as ME_up_count / ME_down_count) is also counted as
exclusion for any other microexon lying inside that junction. Only the
immediate junction is used, so this works for runs of any length.

The same sharing affects inclusion: a read on the m1|m2 junction fits both
m1's and m2's inclusion tags, and Bowtie (-k 1) reports only one, so the other
microexon loses it. A read on one microexon's tag that crosses the junction
shared with an adjacent microexon is therefore also counted as inclusion for
that neighbour (included_neighbours).
"""

import bisect
from collections import defaultdict

MIN_ANCHOR = 8


def parse_intron(intron):
    """chr1:200+1000 -> ("chr1", "+", 200, 1000)."""
    chrom, coords = intron.rsplit(":", 1)
    strand = "+" if "+" in coords else "-"
    start, end = coords.split(strand)
    return chrom, strand, int(start), int(end)


def parse_microexon(microexon):
    """chr1_+_500_512 -> ("chr1", "+", 500, 512); chromosome names may contain _."""
    fields = microexon.split("_")
    return "_".join(fields[:-3]), fields[-3], int(fields[-2]), int(fields[-1])


def junction_id(chrom, strand, start, end):
    return "{}:{}{}{}".format(chrom, start, strand, end)


class MicroexonIndex:
    """Find microexons that lie entirely inside a genomic interval."""

    def __init__(self, microexons):
        by_strand = defaultdict(list)
        for microexon in set(microexons):
            chrom, strand, start, end = parse_microexon(microexon)
            by_strand[(chrom, strand)].append((start, end, microexon))
        self._entries = {key: sorted(values) for key, values in by_strand.items()}
        self._starts = {key: [e[0] for e in values] for key, values in self._entries.items()}
        self._by_start = {}
        self._by_end = {}
        for (chrom, strand), values in self._entries.items():
            for start, end, microexon in values:
                self._by_start.setdefault((chrom, strand, start), []).append((microexon, end - start))
                self._by_end.setdefault((chrom, strand, end), []).append((microexon, end - start))

    def starting_at(self, chrom, strand, position):
        return self._by_start.get((chrom, strand, position), [])

    def ending_at(self, chrom, strand, position):
        return self._by_end.get((chrom, strand, position), [])

    def inside(self, chrom, strand, start, end):
        key = (chrom, strand)
        if key not in self._entries:
            return []
        entries, starts = self._entries[key], self._starts[key]
        found = []
        for index in range(bisect.bisect_left(starts, start), len(entries)):
            me_start, me_end, microexon = entries[index]
            if me_start >= end:
                break
            if me_end <= end:
                found.append(microexon)
        return found


def crossed_junctions(intron, microexon, read_start, matches, anchor_up, me_len):
    """Genomic junctions of an inclusion tag that a read crosses with 8 nt anchors.

    Tag coordinates run in transcript order: anchor_up nt of upstream sequence,
    the microexon, then the downstream sequence.
    """
    chrom, strand, istart, iend = parse_intron(intron)
    _, _, me_start, me_end = parse_microexon(microexon)
    read_end = read_start + matches
    if strand == "+":
        upstream, downstream = (istart, me_start), (me_end, iend)
    else:
        upstream, downstream = (me_end, iend), (istart, me_start)
    junctions = []
    if read_start <= anchor_up - MIN_ANCHOR and read_end >= anchor_up + MIN_ANCHOR:
        junctions.append(upstream)
    if read_start <= anchor_up + me_len - MIN_ANCHOR and read_end >= anchor_up + me_len + MIN_ANCHOR:
        junctions.append(downstream)
    return [(chrom, strand, start, end) for start, end in junctions]


def skipped_by_read(index, intron, microexon, read_start, matches, anchor_up, me_len):
    """{skipped microexon: junction ID} for one read on `microexon`'s inclusion tag."""
    skipped = {}
    for chrom, strand, start, end in crossed_junctions(
        intron, microexon, read_start, matches, anchor_up, me_len
    ):
        for other in index.inside(chrom, strand, start, end):
            if other != microexon:
                skipped.setdefault(other, junction_id(chrom, strand, start, end))
    return skipped


def included_neighbours(index, intron, microexon, read_start, matches, anchor_up, me_len):
    """{adjacent microexon: weight} for one read on `microexon`'s inclusion tag.

    A read crossing the junction between this microexon and an adjacent one
    is inclusion evidence for the neighbour too. Weights are in the half-read
    units that correct_quant.py uses for inclusion reads: 1 when the read
    crosses one of the neighbour's junctions, 2 when it spans the neighbour
    (both junctions with 8 nt anchors).
    """
    chrom, strand, istart, iend = parse_intron(intron)
    read_end = read_start + matches
    tag_down = anchor_up + me_len
    # Genomic coordinate where an upstream / downstream neighbour would end / start.
    if strand == "+":
        upstream = index.ending_at(chrom, strand, istart)
        downstream = index.starting_at(chrom, strand, iend)
    else:
        upstream = index.starting_at(chrom, strand, iend)
        downstream = index.ending_at(chrom, strand, istart)
    found = {}
    if read_start <= anchor_up - MIN_ANCHOR and read_end >= anchor_up + MIN_ANCHOR:
        for other, other_len in upstream:
            if other != microexon:
                found[other] = 2 if read_start <= anchor_up - other_len - MIN_ANCHOR else 1
    if read_start <= tag_down - MIN_ANCHOR and read_end >= tag_down + MIN_ANCHOR:
        for other, other_len in downstream:
            if other != microexon:
                found[other] = 2 if read_end >= tag_down + other_len + MIN_ANCHOR else 1
    return found


def twin_intron_index(introns, microexon):
    """Index of the neighbour's intron whose tag also contains `microexon`'s side.

    For a neighbour X of microexon Y, that is X's containing intron ending at
    Y's start or starting at Y's end; the read could equally have aligned to
    X's tag on that intron. Falls back to the first intron.
    """
    _, _, me_start, me_end = parse_microexon(microexon)
    for position, intron in enumerate(introns):
        _, _, start, end = parse_intron(intron)
        if end == me_start or start == me_end:
            return position
    return 0

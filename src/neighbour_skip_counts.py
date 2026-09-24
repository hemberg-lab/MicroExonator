"""Exclusion evidence for consecutive microexons from their neighbours' tags.

With two or more microexons in a row (E1-m1-m2-E2), the junctions that skip
only one of them (E1|m2 skips m1, m1|E2 skips m2) are the neighbour's own
inclusion junctions. They never get an exon-exon tag, so skipping reads land
on the neighbour's inclusion tag and were only counted as inclusion for the
neighbour. Here, every read that crosses a junction of an inclusion tag with
8 nt anchors (the same rule as ME_up_count / ME_down_count) is also counted as
exclusion for any other microexon lying inside that junction. Only the
immediate junction is used, so this works for runs of any length.
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

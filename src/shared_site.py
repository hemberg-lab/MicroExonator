"""Microexons that share one splice site with a longer annotated exon.

When a microexon M shares its acceptor (or donor) with a longer exon L, reads
that cross the shared junction and end inside M's length fit both M's
inclusion tag and L, so they cannot be assigned. Reads crossing M's other
junction come only from M, and in an M molecule both junctions yield the same
number of single-junction reads. M's inclusion is therefore counted as its
full-spanning reads plus the single-junction reads of the unshared side
(the "mirror" count), and L enters exclusion once, through the junction on
its far side, which M never uses.

Tags are written in transcript orientation, so a shared acceptor is the "up"
junction of M's inclusion tag and a shared donor the "down" junction.
"""

import re


def shared_sides(shares_acceptor, shares_donor):
    """Inclusion-tag sides ("up"/"down") that a longer exon also uses."""
    sides = set()
    if shares_acceptor:
        sides.add("up")
    if shares_donor:
        sides.add("down")
    return sides


def unique_side_half_reads(up_only, down_only, sides):
    """Single-junction reads that can only come from the microexon.

    With both sides shared there is no clean side; the smaller count is the
    one less contaminated by the longer exons.
    """
    if sides == {"up"}:
        return down_only
    if sides == {"down"}:
        return up_only
    return min(up_only, down_only)


def exon_bounds(exon):
    """(chrom, start, end) from the annotation key 'chrom:start<strand>end'."""
    chrom, start, end = re.match(r"^(.*):(\d+)[+-](\d+)$", exon).groups()
    return chrom, int(start), int(end)


def far_and_shared_introns(exon, strand, shares_acceptor, estart_introns, eend_introns):
    """Introns on the longer exon's unshared side, and on its shared side.

    estart_introns maps 'chrom_coord' to introns ending at that coordinate
    (joining the exon starting there); eend_introns to introns starting at it.
    """
    chrom, start, end = exon_bounds(exon)
    low = estart_introns.get(chrom + "_" + str(start), set())
    high = eend_introns.get(chrom + "_" + str(end), set())
    # acceptor is the low coordinate on "+", the high one on "-"
    shared_is_low = shares_acceptor == (strand == "+")
    return (high, low) if shared_is_low else (low, high)


def long_exon_units(alternatives, strand, shares_acceptor, estart_introns, eend_introns,
                    coverage, tagged_introns, credited):
    """Exclusion reads for each longer exon, counting every intron once.

    Exons are counted through their tagged far-side introns. Only when none
    of them has one do they fall back to their shared-side introns (mixing
    the two would count the same molecules twice). Introns already in
    `credited` (shared with an exon counted earlier for the same microexon)
    add nothing. Returns one count per exon, in the order given.
    """
    sides = [far_and_shared_introns(exon, strand, shares_acceptor, estart_introns, eend_introns)
             for exon in alternatives]
    any_far = any(i in tagged_introns for far, _ in sides for i in far)
    counts = []
    for far, shared in sides:
        use = sorted(i for i in (far if any_far else shared) if i in tagged_introns)
        total = 0
        for intron in use:
            if intron not in credited:
                credited.add(intron)
                total += coverage[intron]
        counts.append(total)
    return counts

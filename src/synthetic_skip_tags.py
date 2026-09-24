"""Synthetic skipping tags for microexons whose skipping junction has no tag.

Skipping (exon-exon) tags normally come from annotated junctions only
(Round1/ME_TAGs.fa). A microexon whose containing intron is not an annotated
junction, or whose flanking exons are too short for that library, then has no
exclusion evidence and its PSI is pinned at 1. For each such intron we write
one tag joining the same flanks used for the inclusion tag, without the
microexon. Introns that start or end at another microexon are skipped: their
skipping junction is the neighbour's own inclusion junction, and a synthetic
tag would duplicate the neighbour's inclusion tag.
"""

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def load_tag_introns(path):
    """Return the intron IDs of two-anchor (exon-exon) tags in a tag FASTA."""
    introns = set()
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            fields = line[1:].strip().split("|")
            if len(fields) == 3 and len(fields[2].split("_")) == 2:
                introns.add(fields[0])
    return introns


def intron_id(chrom, start, strand, end):
    """MicroExonator's junction ID, e.g. chr1:200+1000."""
    return "{}:{}{}{}".format(chrom, start, strand, end)


def microexon_boundaries(microexons):
    """Index (chrom, strand, start) and (chrom, strand, end) of each microexon.

    microexons: iterable of (chrom, strand, start, end), zero-based half-open.
    """
    starts, ends = set(), set()
    for chrom, strand, start, end in microexons:
        starts.add((chrom, strand, int(start)))
        ends.add((chrom, strand, int(end)))
    return starts, ends


def parse_microexon_id(microexon):
    """Split chrom_strand_start_end, allowing underscores in the chromosome."""
    strand, start, end = microexon.split("_")[-3:]
    chrom = "_".join(microexon.split("_")[:-3])
    return chrom, strand, int(start), int(end)


def needs_synthetic_tag(chrom, strand, start, end, existing, written, boundaries):
    """True when intron start-end needs a synthetic skipping tag.

    existing: intron IDs that already have a skipping tag.
    written: intron IDs that already received a synthetic tag.
    boundaries: (starts, ends) from microexon_boundaries().
    """
    junction = intron_id(chrom, start, strand, end)
    if junction in existing or junction in written:
        return False
    starts, ends = boundaries
    if (chrom, strand, int(end)) in starts or (chrom, strand, int(start)) in ends:
        return False
    return True


def skip_tag_record(chrom, strand, start, end, transcript, up_sequence, down_sequence):
    """FASTA header and sequence for the skipping tag of one intron.

    up_sequence / down_sequence are the genomic-orientation flanks ending at the
    intron start and starting at the intron end. Anchors follow the inclusion
    tags: tag-orientation lengths, reversed for minus-strand genes.
    """
    sequence = str(up_sequence) + str(down_sequence)
    anchors = [len(up_sequence), len(down_sequence)]
    if strand == "-":
        sequence = reverse_complement(sequence)
        anchors = anchors[::-1]
    header = "|".join(
        [intron_id(chrom, start, strand, end), transcript, "_".join(map(str, anchors))]
    )
    return header, sequence.upper()

"""Annotated introns that contain each microexon, from one bedtools call.

Get_annotated_microexons*.py used to run one `bedtools intersect` per
microexon (a temporary file and a full scan of the intron BED each time). This
runs the same intersect once for all microexons. The query is unchanged: the
introns (A) are intersected with the microexon interval (B) written as
`chrom start end-1 ME 0 strand`, same strand, B fully covered (-F 1). Each
microexon's hits come back as the intron lines it got before, in the same
(intron file) order.
"""

from pybedtools import BedTool


def containing_introns(intron_bed, microexons):
    """{(chrom, start, end, strand, ...) key: [intron line, ...]} for each microexon key.

    `microexons` yields keys whose first four fields are chrom, start, end and
    strand (0-based start, end exclusive). Microexons with no containing intron
    are absent from the result.
    """
    keys = list(microexons)
    if not keys:
        return {}
    query = "\n".join(
        " ".join([key[0], str(key[1]), str(key[2] - 1), str(index), "0", key[3]])
        for index, key in enumerate(keys)
    )
    hits = {}
    for interval in intron_bed.intersect(BedTool(query, from_string=True),
                                         wa=True, wb=True, s=True, F=1, nonamecheck=True):
        fields = str(interval).rstrip("\n").split("\t")
        hits.setdefault(keys[int(fields[9])], []).append("\t".join(fields[:6]))
    return hits

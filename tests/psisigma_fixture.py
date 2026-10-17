"""A tiny HISAT2-like paired-end fixture for the PSI-Sigma module.

One gene on chr1 (+ strand):
  T1: E1 101-300, E2 1001-1060, E3 2001-2200        (cassette exon E2)
  T2: 101-1060, 2001-2200                           (retains intron 301-1000)
  T3: E1, E3                                        (skips E2)
Reads are written as SAM with HISAT2's conventions (MAPQ 60 and NH:i:1 for
unique alignments, MAPQ 1 and NH:i:2 with a secondary for multi-mappers) and
cover the cases PSI-Sigma's counting distinguishes: overlapping mates on one
junction, multi-mappers, soft-clipped, deleted and inserted bases, single
reads in intron-retention positions.
"""

import random

LENGTH = 3000
EXON = {"E1": (101, 300), "E2": (1001, 1060), "E3": (2001, 2200)}


def genome(seed=7):
    rng = random.Random(seed)
    bases = [rng.choice("ACGT") for _ in range(LENGTH)]
    # canonical motifs for every annotated intron (+ strand GT...AG)
    for start, end in ((301, 1000), (1061, 2000), (301, 2000)):
        bases[start - 1:start + 1] = list("GT")
        bases[end - 2:end] = list("AG")
    # a non-canonical junction 1500-1599 used by a few reads
    bases[1499:1501] = list("CA")
    bases[1597:1599] = list("TT")
    return "".join(bases)


def gtf(gene_lines=True, transcript_lines=True):
    lines = []
    attrs = 'gene_id "G1"; gene_name "GENE1";'
    if gene_lines:
        lines.append("chr1\ttest\tgene\t101\t2200\t.\t+\t.\t" + attrs)
    structures = {"T1": [EXON["E1"], EXON["E2"], EXON["E3"]], "T2": [(101, 1060), EXON["E3"]],
                  "T3": [EXON["E1"], EXON["E3"]]}
    for transcript, exons in structures.items():
        extra = ' transcript_id "{}"; transcript_type "protein_coding";'.format(transcript)
        if transcript_lines:
            lines.append("chr1\ttest\ttranscript\t{}\t{}\t.\t+\t.\t{}{}".format(
                exons[0][0], exons[-1][1], attrs, extra))
        for start, end in exons:
            lines.append("chr1\ttest\texon\t{}\t{}\t.\t+\t.\t{}{}".format(start, end, attrs, extra))
    return "\n".join(lines) + "\n"


def _record(name, flag, pos, cigar, sequence, mapq, nh, mate_pos, tlen, extra=()):
    tags = ["NH:i:{}".format(nh)] + list(extra)
    return "\t".join([name, str(flag), "chr1", str(pos), str(mapq), cigar, "=", str(mate_pos), str(tlen),
                      sequence, "I" * len(sequence)] + tags)


def _seq(genome_, blocks):
    return "".join(genome_[start - 1:start - 1 + length] for start, length in blocks)


def sam(sample_psi, seed=1, depth=60):
    """SAM text for one sample; `sample_psi` (0-1) sets E2 inclusion and IR."""
    rng = random.Random(seed)
    g = genome()
    records = []

    def pair(name, first, second, nh=1, secondary=None):
        mapq = 60 if nh == 1 else 1
        (p1, c1, s1), (p2, c2, s2) = first, second
        span = p2 + len(s2) - p1
        records.append(_record(name, 99, p1, c1, s1, mapq, nh, p2, span))
        records.append(_record(name, 147, p2, c2, s2, mapq, nh, p1, -span))
        if secondary:
            records.append(_record(name, 355, secondary[0], secondary[1], secondary[2], mapq, nh, p1, 0))

    count = 0
    for _ in range(depth):
        count += 1
        name = "r{}".format(count)
        draw = rng.random()
        if draw < sample_psi:
            # inclusion: mate 1 E1->E2, mate 2 E2->E3; sometimes both cross E1->E2 (overlap)
            first = (276, "25M700N25M", _seq(g, [(276, 25), (1001, 25)]))
            if rng.random() < 0.3:
                second = (281, "20M700N30M", _seq(g, [(281, 20), (1001, 30)]))
            else:
                second = (1036, "25M940N25M", _seq(g, [(1036, 25), (2001, 25)]))
            pair(name, first, second)
        else:
            first = (276, "25M1700N25M", _seq(g, [(276, 25), (2001, 25)]))
            second = (2010, "50M", _seq(g, [(2010, 50)]))
            pair(name, first, second)
    # intronic reads (retained intron 301-1000), unique, some with D/I or clipping
    for index in range(int(depth * sample_psi * (0.7 + 0.15 * (seed % 4)))):
        count += 1
        name = "r{}".format(count)
        start = 301 + rng.randrange(0, 600)
        kind = index % 5
        if kind == 1:
            first = (start, "20M2D28M", _seq(g, [(start, 20), (start + 22, 28)]))
        elif kind == 2:
            first = (start, "20M2I28M", _seq(g, [(start, 20)]) + "AA" + _seq(g, [(start + 20, 28)]))
        elif kind == 3:
            first = (start, "5S45M", "NNNNN" + _seq(g, [(start, 45)]))
        else:
            first = (start, "50M", _seq(g, [(start, 50)]))
        second = (start + 100, "50M", _seq(g, [(start + 100, 50)]))
        pair(name, first, second)
    # multi-mappers: one spliced (E1->E3), one intronic, each with a secondary
    for _ in range(5):
        count += 1
        pair("r{}".format(count), (276, "25M1700N25M", _seq(g, [(276, 25), (2001, 25)])),
             (2010, "50M", _seq(g, [(2010, 50)])), nh=2,
             secondary=(276, "25M1700N25M", _seq(g, [(276, 25), (2001, 25)])))
        count += 1
        pair("r{}".format(count), (500, "50M", _seq(g, [(500, 50)])), (600, "50M", _seq(g, [(600, 50)])), nh=2,
             secondary=(700, "50M", _seq(g, [(700, 50)])))
    # a few reads over a non-canonical junction 1500-1599
    for _ in range(12):
        count += 1
        pair("r{}".format(count), (1475, "25M100N25M", _seq(g, [(1475, 25), (1600, 25)])),
             (1700, "50M", _seq(g, [(1700, 50)])))
    records.sort(key=lambda line: int(line.split("\t")[3]))
    header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{}\n".format(LENGTH)
    return header + "\n".join(records) + "\n"


def mapq255(sam_text):
    """The same SAM as STAR would mark it: MAPQ 255 for NH:i:1, else 0."""
    lines = []
    for line in sam_text.splitlines():
        if line.startswith("@"):
            lines.append(line)
            continue
        fields = line.split("\t")
        fields[4] = "255" if "NH:i:1" in fields[11:] else "0"
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"

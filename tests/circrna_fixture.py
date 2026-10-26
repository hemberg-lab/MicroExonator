"""A tiny genome with two genes, two circRNAs and simulated paired-end reads.

chr1, gene G1 (+): exons 201-400, 1001-1150, 1601-1800, 2601-2800 (1-based, inclusive)
    circRNA C1 = exons 2-3, back-spliced junction joins 1800 -> 1001
chr2, gene G2 (-): exons 3201-3400, 2201-2400, 1201-1400 (transcript order)
    circRNA C2 = exon 2201-2400 alone (single-exon circle), BSJ 2400 -> 2201
Introns carry canonical GT-AG motifs in the transcript's direction.

Reads are unstranded (mates in random orientation), 2 x 75 nt from 220-260 nt
fragments. `expected_bsj` counts the fragments of each circle in which a mate crosses
the back-spliced junction with at least `anchor` bases on each side; these are the
reads every tool should be able to see.
"""

import gzip
import random

CHROM_LENGTH = 4000
READ = 75
G1_EXONS = [(201, 400), (1001, 1150), (1601, 1800), (2601, 2800)]
G2_EXONS = [(3201, 3400), (2201, 2400), (1201, 1400)]           # transcript order on -
C1 = ("chr1", 1001, 1800, "+", [(1001, 1150), (1601, 1800)])
C2 = ("chr2", 2201, 2400, "-", [(2201, 2400)])
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def genome(seed=11):
    rng = random.Random(seed)
    chroms = {}
    for chrom, exons, strand in (("chr1", G1_EXONS, "+"), ("chr2", sorted(G2_EXONS), "-")):
        bases = [rng.choice("ACGT") for _ in range(CHROM_LENGTH)]
        for (_, left_end), (right_start, _) in zip(exons, exons[1:]):
            # intron = left_end+1 .. right_start-1 (1-based); motif in transcript direction
            donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
            bases[left_end:left_end + 2] = list(donor)
            bases[right_start - 3:right_start - 1] = list(acceptor)
        chroms[chrom] = "".join(bases)
    # flanking motifs for the circRNAs' outer splice sites (already intronic above for C1,
    # and for C2 its neighbours are exons, so the introns around it carry the motifs)
    return chroms


def write_genome(path, chroms):
    with open(path, "w") as stream:
        for name, sequence in chroms.items():
            stream.write(">{}\n".format(name))
            for i in range(0, len(sequence), 60):
                stream.write(sequence[i:i + 60] + "\n")
    lines = []
    offset = 0
    with open(path) as stream:
        text = stream.read()
    for name, sequence in chroms.items():
        header = ">{}\n".format(name)
        offset = text.index(header, offset) + len(header)
        lines.append("{}\t{}\t{}\t60\t61\n".format(name, len(sequence), offset))
    with open(str(path) + ".fai", "w") as stream:
        stream.write("".join(lines))


def gtf():
    lines = []
    for gene, chrom, strand, exons in (("G1", "chr1", "+", G1_EXONS), ("G2", "chr2", "-", G2_EXONS)):
        ordered = sorted(exons)
        attrs = 'gene_id "{0}"; transcript_id "{0}.1"; gene_name "{0}"; transcript_type "protein_coding";'.format(gene)
        lines.append("\t".join([chrom, "test", "gene", str(ordered[0][0]), str(ordered[-1][1]), ".", strand, ".",
                                'gene_id "{0}"; gene_name "{0}";'.format(gene)]))
        lines.append("\t".join([chrom, "test", "transcript", str(ordered[0][0]), str(ordered[-1][1]), ".",
                                strand, ".", attrs]))
        for start, end in ordered:
            lines.append("\t".join([chrom, "test", "exon", str(start), str(end), ".", strand, ".", attrs]))
    return "\n".join(lines) + "\n"


def transcript_sequence(chroms, chrom, strand, exons):
    sequence = "".join(chroms[chrom][start - 1:end] for start, end in sorted(exons))
    return revcomp(sequence) if strand == "-" else sequence


def catalogue_bed():
    """BED6 (0-based start): both true circRNAs, an unexpressed decoy, and rows to drop."""
    return "".join([
        "chr1\t1000\t1800\tcircAtlas_C1\t0\t+\n",
        "chr1\t1000\t1800\tcircBase_C1\t0\t+\n",          # same BSJ from a second source
        "chr2\t2200\t2400\tcircAtlas_C2\t0\t-\n",
        "chr1\t200\t1150\tdecoy_E1E2\t0\t+\n",             # not expressed as a circle
        "chrUn\t10\t500\tunknown_contig\t0\t+\n",
        "chr1\t3900\t4200\toutside\t0\t+\n",
        "chr2\t1200\t1400\tboth_strands\t0\t+\n",
        "chr2\t1200\t1400\tboth_strands\t0\t-\n",
    ])


def simulate(chroms, seed=5, linear_fragments=600, circle_fragments=(160, 120), anchor=10):
    rng = random.Random(seed)
    pairs = []
    expected = {}

    def emit(name, fragment):
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        pairs.append((name, fragment[:READ], revcomp(fragment[-READ:])))

    t1 = transcript_sequence(chroms, "chr1", "+", G1_EXONS)
    t2 = transcript_sequence(chroms, "chr2", "-", G2_EXONS)
    for i in range(linear_fragments):
        transcript = t1 if i % 2 == 0 else t2
        size = rng.randint(220, 260)
        start = rng.randint(0, len(transcript) - size)
        emit("f{:06d}".format(len(pairs)), transcript[start:start + size])
    for (chrom, start, end, strand, exons), count in zip((C1, C2), circle_fragments):
        circle = transcript_sequence(chroms, chrom, strand, exons)
        length = len(circle)
        doubled = circle * 3
        spanning = 0
        for _ in range(count):
            size = rng.randint(150, min(200, length - 10))
            offset = rng.randint(0, length - 1)
            fragment = doubled[offset:offset + size]
            # the junction sits at multiples of `length` in transcript coordinates
            mates = [(offset, offset + READ), (offset + size - READ, offset + size)]
            if any(any(left + anchor <= junction <= right - anchor for junction in (length, 2 * length))
                   for left, right in mates):
                spanning += 1
            emit("f{:06d}".format(len(pairs)), fragment)
        expected["{}:{}|{}".format(chrom, start, end)] = spanning
    rng.shuffle(pairs)
    return pairs, expected


def star_chimeric_lines(reads_per_circle=(12, 9), seed=3, read=READ):
    """Emulated STAR Chimeric.out.junction lines (legacy 14-column format, single-end style).

    Used only where STAR itself cannot run (STAR's macOS builds cannot spawn
    readFilesCommand); the cluster smoke test runs real STAR. Column semantics follow the
    STAR manual: brkpt_donorA / brkpt_acceptorB are the first intron bases next to the
    donor and acceptor (1-based), junction type 1 = GT/AG, start_aln 1-based, CIGARs on
    the + strand of the genome."""
    rng = random.Random(seed)
    lines = []
    for (chrom, start, end, strand, _), count in zip((C1, C2), reads_per_circle):
        for i in range(count):
            x = rng.randint(20, read - 20)          # bases on the circle's 3' side of the BSJ
            name = "chim{}_{:03d}".format(chrom, i)
            if strand == "+":
                # donor = last base `end`, acceptor = first base `start`
                fields = [chrom, end + 1, "+", chrom, start - 1, "+", 1, 0, 0, name,
                          end - x + 1, "{}M{}S".format(x, read - x), start, "{}S{}M".format(x, read - x)]
            else:
                # transcript runs end -> start: donor side at `start`, acceptor at `end`
                fields = [chrom, start - 1, "-", chrom, end + 1, "-", 2, 0, 0, name,
                          start, "{}S{}M".format(read - x, x), end - (read - x) + 1, "{}M{}S".format(read - x, x)]
            lines.append("\t".join(str(value) for value in fields) + "\n")
    return "".join(lines), {"{}:{}|{}".format(c[0], c[1], c[2]): n for c, n in zip((C1, C2), reads_per_circle)}


def write_fastq(pairs, r1, r2):
    with gzip.open(r1, "wt") as one, gzip.open(r2, "wt") as two:
        for name, first, second in pairs:
            one.write("@{}\n{}\n+\n{}\n".format(name, first, "I" * len(first)))
            two.write("@{}\n{}\n+\n{}\n".format(name, second, "I" * len(second)))


def build(directory, seed=5):
    """Write genome.fa(+fai), annotation.gtf, catalogue.bed6, R1/R2; return paths and truth."""
    from pathlib import Path
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    chroms = genome()
    write_genome(directory / "genome.fa", chroms)
    (directory / "annotation.gtf").write_text(gtf())
    (directory / "catalogue.bed6").write_text(catalogue_bed())
    pairs, expected = simulate(chroms, seed=seed)
    write_fastq(pairs, directory / "R1.fastq.gz", directory / "R2.fastq.gz")
    return {"genome": directory / "genome.fa", "gtf": directory / "annotation.gtf",
            "catalogue": directory / "catalogue.bed6", "r1": directory / "R1.fastq.gz",
            "r2": directory / "R2.fastq.gz", "expected_bsj": expected, "chroms": chroms}

"""Flag microexon junctions whose sequence also exists unspliced in the genome.

Reads that align to the genome without splicing are removed before counting
(Filter1_round2.py). When a junction's sequence has an exact copy elsewhere,
typically a processed pseudogene (retrocopy) of one isoform, reads from that
junction are indistinguishable from the copy and are removed too. A copy of the
skipping junction pushes PSI up; a copy of the inclusion junction pushes it
down. This script lists those events so their PSI can be read with care.

Usage:
  junction_genomic_copies.py cores TAGS_FASTA ME_CENTRIC > cores.fa
  junction_genomic_copies.py report CORES_FASTA BOWTIE_HITS > report.tsv
"""

import csv
import sys

FLANK = 30


def read_fasta(path):
    name, chunks = None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def event_junctions(me_centric):
    """(ME, intron, microexon sequence) for every containing intron of every event."""
    csv.field_size_limit(1000000000)
    for row in csv.reader(open(me_centric), delimiter="\t"):
        total_SJs, micro_exon_seq_found, total_ME = row[3], row[6], row[11]
        if total_ME == "":
            continue
        true_ME = max([i.split("|") for i in total_ME.split(",")], key=lambda item: float(item[1]))[0]
        for intron in total_SJs.split(","):
            yield true_ME, intron, micro_exon_seq_found


def junction_cores(tags_fasta, me_centric):
    """Core sequences (30 nt each side of the junction) for skipping and inclusion tags."""
    wanted = {}
    for microexon, intron, me_seq in event_junctions(me_centric):
        wanted.setdefault((intron, ""), set()).add(microexon)
        wanted.setdefault((intron, me_seq), set()).add(microexon)
    cores = {}
    for name, seq in read_fasta(tags_fasta):
        intron, _, anchors = name.split("|")
        fields = anchors.split("_")
        me_seq = fields[1] if len(fields) == 3 else ""
        if (intron, me_seq) not in wanted:
            continue
        anchor_up = int(fields[0])
        if anchor_up < FLANK or len(seq) < anchor_up + len(me_seq) + FLANK:
            continue
        core = seq[anchor_up - FLANK:anchor_up + len(me_seq) + FLANK]
        kind = "inclusion" if me_seq else "skipping"
        for microexon in wanted[(intron, me_seq)]:
            cores["{}|{}|{}".format(microexon, kind, intron)] = core
    return cores


def report(cores_fasta, bowtie_hits):
    """Group bowtie (default output) hits by event and junction type."""
    names = [name for name, _ in read_fasta(cores_fasta)]
    copies = {}
    for row in csv.reader(open(bowtie_hits), delimiter="\t"):
        if len(row) < 4:
            continue
        copies.setdefault(row[0], []).append("{}:{}:{}".format(row[2], int(row[3]) + 1, row[1]))
    lines = ["ME\tjunction\tintron\tgenomic_copies"]
    for name in names:
        if name in copies:
            microexon, kind, intron = name.split("|")
            lines.append("\t".join([microexon, kind, intron, ",".join(sorted(copies[name]))]))
    return lines


if __name__ == "__main__":
    if sys.argv[1] == "cores":
        for name, core in junction_cores(sys.argv[2], sys.argv[3]).items():
            print(">{}\n{}".format(name, core))
    elif sys.argv[1] == "report":
        print("\n".join(report(sys.argv[2], sys.argv[3])))
    else:
        sys.exit(__doc__)

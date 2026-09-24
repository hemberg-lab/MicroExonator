import sys
import csv
from collections import defaultdict

csv.field_size_limit(100000000)


def load_tag_sequences(tags_fasta):
    tags = {}
    name = None
    chunks = []
    for line in open(tags_fasta):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                tags[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line)
    if name is not None:
        tags[name] = "".join(chunks)
    return tags


def tag_mismatches(tags, tag, start, seq):
    """Mismatches of a read (tag orientation, 0-based start) against its tag."""
    reference = tags.get(tag, "")[int(start):int(start) + len(seq)]
    if len(reference) != len(seq):
        return len(seq)
    return sum(1 for a, b in zip(reference.upper(), seq.upper()) if a != b)


def main(pre_processed, genome_sam, tags_fasta="NA"):

    # By default a read is removed if it aligns anywhere in the genome
    # without splicing. With a tag FASTA, it is removed only when its genome
    # alignment has no more mismatches than its tag alignment, so reads that
    # match a splice junction better than a retrocopy of it are kept.
    tags = load_tag_sequences(tags_fasta) if tags_fasta != "NA" else None
    genome_mismatches = {}

    read_SJ = defaultdict(set)
    black_list = set([])

    # for row in csv.reader(open(dust), delimiter = '>'):

    #   black_list.add(row[1])

    # for row in csv.reader(open(repbase), delimiter = '\t'):

    #   black_list.add(row[9])

    for row in csv.reader(open(genome_sam), delimiter = '\t'):
        try:
            if row[1]=="0" or row[1]=="16":
                black_list.add(row[0])
                nm = [int(f[5:]) for f in row[11:] if f.startswith("NM:i:")]
                genome_mismatches[row[0]] = min(genome_mismatches.get(row[0], 99), nm[0] if nm else 0)
        except (IndexError, ValueError):
            pass

    for row in csv.reader(open(pre_processed), delimiter = '\t'):
        try:
            read, flag, tag, start, cigar, seq, qual = row

            SJ = tag.split("|")[0]
            read_SJ[read].add(SJ)
        except ValueError:
            pass

    for row in csv.reader(open(pre_processed), delimiter = '\t'):
        try:
            read, flag, tag, start, cigar, seq, qual = row

            #if (read in black_list)==False and len(read_SJ[read])==1:
            if read not in black_list:
                print(("\t".join(row)))
            elif tags is not None and genome_mismatches[read] > tag_mismatches(tags, tag, start, seq):
                print(("\t".join(row)))
        except ValueError:
            pass
        #print black_list


main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "NA") #, sys.argv[3], sys.argv[4])

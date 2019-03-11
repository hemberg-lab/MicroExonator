import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}

def Genomictabulator(fasta):

	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	print >> sys.stderr, "OK"

	f.close()


def main(bed12):

    transcripts_seq = defaultdict(str)

    for row in csv.reader(open(bed12), delimiter = '\t'):
        start = int(row[1])
        end = int(row[2])
        strand = row[5]
        bn = int(row[9])
        chrom = row[0]
        transcript = row[3]
        blocksizes = map(int, row[10].strip(",").split(","))
        qstarts = map (int, row[11].strip(",").split(","))

        for q, b in zip(qstarts, blocksizes):


            estart = start + q
            eend = start + q + b
            elength = eend - estart

            exon_seq = Genome[chrom][estart:eend]

            transcripts_seq[(transcript, strand)] += exon_seq


    for key_value in transcripts_seq.items():

        transcript_strand, seq = key_value
        transcript, strand = transcript_strand

        if strand=="-":
            seq = seq.reverse_complement()

        seq = str(seq).upper()

        print(">" + transcript)
        print(seq)


if __name__ == '__main__':
    Genomictabulator(sys.argv[1])
    main(sys.argv[2])

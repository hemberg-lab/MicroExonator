
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main(fastq1, fastq2, fastq12):

    r1 = open(fastq1)
    r2 = open(fastq2)
    out = open(fastq12, 'w')

    for read in SeqIO.parse(r1, "fastq"):

        rd_id = read.id + "_1"
        rd_Q = read.letter_annotations["phred_quality"]
        reads_rd = SeqRecord( read.seq, rd_id, description = "" )
        reads_rd.letter_annotations["phred_quality"] = rd_Q
        out.write(reads_rd.format("fastq"))

    for read in SeqIO.parse(r2, "fastq"):

        rd_id = read.id + "_2"
        rd_Q = read.letter_annotations["phred_quality"]
        reads_rd = SeqRecord( read.seq, rd_id, description = "" )
        reads_rd.letter_annotations["phred_quality"] = rd_Q
        out.write(reads_rd.format("fastq"))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])  
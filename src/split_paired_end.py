import sys
import csv
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord



def main (fastq):
	""" Splits the reads at N position. It generate two files rd1 and rd2 """
		
	#N = 80	

	f = open(fastq)
	# rd1 = open(fastq + '.rd1', 'w')
	# rd2 = open(fastq + '.rd2', 'w')
		
	for read in SeqIO.parse(f, "fastq"):
		
		N = len(read.seq)/2

		# if len(read.seq) == 155:
		# 	N = 80 
		
		rd1_seq = read.seq[:N]
		rd2_seq = read.seq[N:]
		
		rd1_id = read.id + "_1"
		rd2_id = read.id + "_2"
				
		rd1_Q = read.letter_annotations["phred_quality"][:N]
		rd2_Q = read.letter_annotations["phred_quality"][N:]
		
		reads_rd1 = SeqRecord( rd1_seq, rd1_id, description = "" )
		reads_rd1.letter_annotations["phred_quality"] = rd1_Q	

		reads_rd2 = SeqRecord( rd2_seq, rd2_id, description = "" )
		reads_rd2.letter_annotations["phred_quality"] = rd2_Q	
		
		# rd1.write(reads_rd1.format("fastq"))
		# rd2.write(reads_rd2.format("fastq"))

		print reads_rd1.format("fastq")
		print reads_rd2.format("fastq")


if __name__ == '__main__':
	main(sys.argv[1])

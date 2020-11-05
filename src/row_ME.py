import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from random import randint, sample
from operator import itemgetter
import re

Genome = {}

def Genomictabulator(fasta):

	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	print >> sys.stderr, "OK"

	f.close()


def main(sam_pre_processed):

	fastq_out = open( ".".join(sys.argv[2].split(".")[:-1]) + ".row_ME.fastq", 'w')



	for row in csv.reader(open(sam_pre_processed), delimiter = '\t'):

		if len(row)==14: #To avoid rare errors (like SRR2138604.sam.pre_processed)
			read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found = row
			intron_tag, transcript_ID, anchors = tag.split("|")

			chr = "_".join(re.findall(r"[\w']+", intron_tag)[:-2])
			istart, iend = re.findall(r"[\w']+", intron_tag)[-2:]

			istart = int(istart)
			iend = int(iend)
			
			try:

				intron_seq = str(Genome[chr][istart:iend]).upper()

				micro_exons_coords = []

				island = "AG" + DR_corrected_micro_exon_seq_found + "GT"
				rev_island = str(Seq(island).reverse_complement())

				strand = "+"

				if "-" in intron_tag:
					strand = "-"

				if strand == "+" and island in intron_seq:

					for i in [i for i in range(len(intron_seq)) if intron_seq.startswith(island, i)]:

						ME_start = i + 2 + istart
						ME_end = ME_start + len(DR_corrected_micro_exon_seq_found)
						ME_chr = chr
						ME_strand = strand

						micro_exons_coords.append("_".join((map(str, [ME_chr, ME_strand, ME_start, ME_end]))))


				elif strand == "-" and rev_island in intron_seq:


					for i in [i for i in range(len(intron_seq)) if intron_seq.startswith(rev_island, i)]:

						ME_start = i + 2 + istart
						ME_end = ME_start + len(DR_corrected_micro_exon_seq_found)
						ME_chr = chr
						ME_strand = strand

						micro_exons_coords.append("_".join((map(str, [ME_chr, ME_strand, ME_start, ME_end]))))

				micro_exons_coords = ",".join(micro_exons_coords)

				if micro_exons_coords!="":
					print "\t".join(row) + "\t" + micro_exons_coords
                    
                    if flag==16:
                        seq = str(Seq(island).reverse_complement())
                        qual = qual[::-1]
					
					if len(seq)==len(qual):

                        fastq_out.write("@" + read + "\n")
                        fastq_out.write(seq + "\n")
                        fastq_out.write("+" + "\n")
                        fastq_out.write(qual + "\n")
                        
					elif len(seq)>len(qual):  ## preventing errors with hisat
                        
                        qual2 = qual + qual[ -(len(seq) - len(qual)) : ]
                            
                        fastq_out.write("@" + read + "\n")
                        fastq_out.write(seq + "\n")
                        fastq_out.write("+" + "\n")
                        fastq_out.write(qual2 + "\n")
                        
                   elif len(seq)<len(qual):
                    
                        qual2 = qual[:len(seq)]
                        
                        fastq_out.write("@" + read + "\n")
                        fastq_out.write(seq + "\n")
                        fastq_out.write("+" + "\n")
                        fastq_out.write(qual2 + "\n")
                        
					# ME_fastq = SeqRecord( seq, id = read, description = "" )
					# ME_fastq.letter_annotations["phred_quality"] = qual

					# fastq_out.write(ME_fastq.format("fastq"))


					# print "@" + read.qname
					# print seq
					# print "+"
					# print q

					
			except KeyError:
				pass 

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2])

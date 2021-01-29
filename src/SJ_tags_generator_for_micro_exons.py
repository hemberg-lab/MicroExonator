import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from random import randint, sample
from operator import itemgetter
from collections import defaultdict
from operator import itemgetter


Transcriptome = {}
Genome = {}

def Genomictabulator(fasta):

	print >> sys.stderr, "Loading the genome into RAM memory ...",

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	print >> sys.stderr, "OK"

	f.close()



def Transcriptometabulator(genecode_fasta):

	print >> sys.stderr, "Loading the genome into RAM memory ...",

	for record in SeqIO.parse(genecode_fasta, "fasta"):
		id = str(record.id).split("|")[0].split(" ")[0]
		Transcriptome[id] = record.seq

		#print(id)

	print >> sys.stderr, "OK"


def main(bed12, ME_len, max_read_len):

	n = max_read_len

	transcript_intron_info = defaultdict(list)

	min_intron_lenght = 80

	for row in csv.reader(open(bed12), delimiter = '\t'):

		try:


			qName = row[3]
			seq = Transcriptome[qName]

			qstarts = map (int, row[11].strip(",").split(","))
			blocksizes = map(int, row[10].strip(",").split(","))

			start = int(row[1])
			strand = row[5]
			bn = int(row[9])
			chr = row[0]
			qstart = 0

			for q1, q2, b, b2 in zip(qstarts, qstarts[1:], blocksizes, blocksizes[1:]):

				qstart = qstart + b
				tag_start = qstart - n
				tag_end = qstart + n

				#if tag_start <= 0:
				#	print tag_start, qstart, tag_end, strand

				istart = start + q1 + b
				iend = start + q2
				ilen = iend - istart
				intron = row[0] + ":" +  str(istart) + row[5] + str(iend)
				intron = chr + ":" + str(istart) + strand + str(iend)
				ilength = iend - istart

				block_up = n
				block_down = n
				dn = str(Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]).upper()


				if strand == '+' :                          #Para los que aliniean en la hebra +

					if tag_start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
						tag_start = 0
						block_up = qstart

					if tag_end>len(seq):
						tag_end=len(seq)
						block_down = tag_end - qstart


					tag = seq[tag_start:tag_end]


				if strand == '-' :

					dn = str((Genome[chr][istart:(istart+2)] + Genome[chr][(iend-2):iend]).reverse_complement()).upper()

					if tag_end>len(seq):                 #Para los que alinian en la hebra - es todo al inverso
						tag_end=len(seq)
						block_up = tag_end - qstart

					tag = seq[-tag_end:-tag_start]

					if tag_start<=0:

						tag = seq[-tag_end:]
						block_down = qstart


				if b > ME_len and b2 > ME_len and ilength >= min_intron_lenght and (dn=="GTAG" or dn=="GCAG" or dn=="ATAC"):  # hay que agregarle el filtro de los micro exones!!


					info = qName, tag, chr, istart, iend, strand, block_up, block_down, block_up + block_down
					transcript_intron_info[intron].append(info)


		except KeyError:
			pass


	for i in transcript_intron_info.items():

		infos = i[1]
		intron = i[0]

		qName, tag, chr, istart, iend, strand, block_up, block_down, sum_blocks = max(infos, key=itemgetter(8))


		ID = ">" + intron + "|" + qName + "|" + str(block_up) + "_" + str(block_down)

		print ID
		print tag


#>chr12:3701518+3702264|ENST00000562877.1|100_19
#AGCTTTCTGTTTAGTTGTGTCAATCGCAGGCCACTCTGCTGAGCATCTTCTCCCAGGAGTACCAGAAACACATTAAAAGAACACATGCCAAACATCATACTTCGGAAGCAATTGAAAGT


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Transcriptometabulator(sys.argv[2])
	main (sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))



#El filtro del los intrones canonicos fue anadido despues

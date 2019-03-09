import sys
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from random import randint, sample
from operator import itemgetter
import re


tags = {}


def Tags_indexer(tags_fasta):
	
	print >> sys.stderr, "Cargando a fasta en la ram ...",
	
	for record in SeqIO.parse(genecode_fasta, "tags_fasta"):
		id = str(record.id).split("|")[0]
		tags[id] = record.seq
		
	print >> sys.stderr, "OK"


def ascii_classifier(c): 

	ascii = ord(c)
	if ascii >= 48 and ascii <= 57:
		return 'number'
	else:
		return 'letter'

def cigar_parser(cigar):

	aux_str = ''
	cigar_vars = []

	for c in cigar:
		c_type = ascii_classifier(c)

		if c_type == 'number':
			aux_str += c

		elif c_type == 'letter':
			cigar_vars.append((c, int(aux_str)))
			aux_str = ''

	return cigar_vars

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def DR_counter(ME5U, ME, ME3D):
	
	L = len(ME)
			
	DRU = 0
	DRD = 0

	ME3U = ME5U + ME
	ME5D = ME + ME3D			
		
	try:
		while ME5U[L-1-DRU]==ME3U[L-1-DRU]:
			DRU += 1
			
			if  ME5U[L-1-DRU]!=ME3U[L-1-DRU]: 
				break
	except IndexError:
		pass 
	try:
		while ME5D[DRD]==ME3D[DRD]:
			DRD += 1

			if ME5D[DRD]!=ME3D[DRD]:
				break
	except IndexError:
		pass
	
	return DRU, DRD


def main(sam, stranded):

	micro_exons_found = defaultdict(set)
	total_micro_exons = defaultdict(set)

	reads_micro_exon = []
	reads_tags_intron = defaultdict(set)
	
	for row in csv.reader(open(sam), delimiter = '\t'):

		if row!=[]:

			try:


				if row[0][0] != "@":
					
					read = row[0]
					flag = int(row[1])
					tag = row[2]
					start = int(row[3]) - 1             #Sam es 1 referenciado y es mas comodo trabajar en cordeneadas 0 refereciadas
					cigar = row[5]
					seq = row[9]
					qual = row[10]

					
					alinments = []

					if (flag==0):

						alinments.append((read, flag, tag, start, cigar, seq, qual))

						if row[13].strip("AS:i:") == row[14].strip("XS:i:") and len(row)>=16 and "XA:Z:" in row[-1]:

							c = 0

							for x in row[-1].strip("XA:Z:").split(";"):

								if x!='':
									tag, start, cigar, n = x.split(",")
									if "+" in start:
										start = int(start.strip("+")) - 1
										alinments.append((read, flag, tag, start, cigar, seq, qual))

					elif stranded=="F" and flag==16:

						alinments.append((read, flag, tag, start, cigar, seq, qual))

						seq = str(Seq(seq).reverse_complement())

						cigar_list = cigar_parser(cigar)

						cigar = ""

						for block in cigar_list[::-1]:
							cigar+= "".join(map(str, block[::-1]))


						if row[13].strip("AS:i:") == row[14].strip("XS:i:") and len(row)>=16 and "XA:Z:" in row[-1]:

							c = 0

							for x in row[-1].strip("XA:Z:").split(";"):

								if x!='':
									tag, start, cigar, n = x.split(",")
									if "+" in start:
										start = int(start.strip("+")) - 1
										alinments.append((read, flag, tag, start, cigar, seq, qual))


					for i in alinments:

						read, flag, tag, start, cigar, seq, qual = i
						intron_tag, transcript_ID, anchors = tag.split("|")

						anchor_up, anchor_down = anchors.split("_")
						anchor_up = int(anchor_up)
						anchor_down = int(anchor_down)

						cigar_list = cigar_parser(cigar)

						q_block_starts = [0]
						q_block_ends = []			
						q_block = 0
						var_index = 0
					
						for var in cigar_list:
							var_type = var[0]
							var_value = var[1]
							var_index += 1
							
							if var_type == 'M':
								q_block += var_value						
								
							if var_type == 'D':
								q_block += var_value
													
							if var_type == 'I':
								q_block_ends.append(q_block_starts[-1] + q_block)
								q_block_starts.append(q_block_ends[-1] + var_value)
								q_block += 0
								
							if var_type == 'N':
								q_block = 0
								
							if var_index == len(cigar_list):
								q_block_ends.append(q_block_starts[-1] + q_block)

						if len(q_block_ends)==2: #deveria tener solo 2 q_blocks ya que deberia encontrar una insercion asocciada a un micro-exon


							micro_exon_seq_found = seq[q_block_ends[0]:q_block_starts[1]].upper()
							micro_exon_len_found = len(micro_exon_seq_found)

							micro_exon_seq_up = seq[:q_block_ends[0]].upper()
							micro_exon_seq_down = seq[q_block_starts[1]:].upper()

							if micro_exon_seq_up!="" and micro_exon_seq_down!="":  #Esto evita insersiones terminales

								DRU, DRD = DR_counter(micro_exon_seq_up, micro_exon_seq_found, micro_exon_seq_down)

								I_pos_tag = start + q_block_ends[0]

								
								if (I_pos_tag - DRU <= anchor_up <= I_pos_tag + DRD):

									DR_corrected_micro_exon_seq_found = seq[q_block_ends[0] + (anchor_up-I_pos_tag):q_block_starts[1] + (anchor_up-I_pos_tag)].upper()

									info = (read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found)

									reads_micro_exon.append(info)
									reads_tags_intron[read].add(intron_tag)

			except IndexError:
				pass


	for i in reads_micro_exon:

		read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found = i
		intron_tag, transcript_ID, anchors = tag.split("|")
		chr, istart, iend = re.findall(r"[\w']+", intron_tag)

		q_block_starts = ",".join(map(str, q_block_starts))
		q_block_ends = ",".join(map(str, q_block_ends))

		istart = int(istart)
		iend = int(iend)

		number_of_host_introns = len(reads_tags_intron[read])

		strand = "+"

		if "-" in intron_tag:
			strand = "-"

		if number_of_host_introns == 1:

			print "\t".join(map(str, [read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found]))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])



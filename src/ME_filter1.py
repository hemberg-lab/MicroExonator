import sys
import csv
from collections import defaultdict
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import pyBigWig

csv.field_size_limit(1000000000)


flag_dict = {'73':[1,1], '89':[1,1], '121':[1,-1], '153':[-1,-1], '185':[-1,-1], '137':[-1,1], '99':[1,1], '147':[-1,-1], '83':[1,-1], '163':[-1,1], '67':[1,1], '115':[1,-1], '179':[-1,-1], '81':[1,-1], "161":[-1,1], '97':[1,1], '145':[-1,-1], '65':[1,1], '129':[-1,1], '113':[1,-1], '177':[-1,-1] }

Genome = {}

def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0


def Genomictabulator(fasta):

	print >> sys.stderr, "Cargando genoma en la memoria RAM ...",

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	print >> sys.stderr, "OK"

	f.close()

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

def PWM_to_dict(file):
	reader = csv.reader(open(file), delimiter = '\t')
	header = reader.next()
	header_dict = {}
	col = 0

	matrix = {}

	for name in header:
		header_dict[name] = col
		col += 1

	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	N_freq = []

	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]

		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
		N_freq.append(0)

	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	matrix["N"] = N_freq

	return matrix


#def main(row_ME, reads_genome, dust, repbase, U2_GTAG_5_file, U2_GTAG_3_file, phylop_vertebrates, phylop_primates): old


def main(row_ME, reads_genome, U2_GTAG_5_file, U2_GTAG_3_file, phylop, ME_len):




	phylop_bw = pyBigWig.open(phylop)


	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])

	TOTAL_U2_max_score = U2_GTAG_5_max_score + U2_GTAG_3_max_score


	black_list = set([])
	exons_reads_genome = defaultdict(list)

	# for row in csv.reader(open(dust), delimiter = '>'):

	# 	black_list.add(row[1])

	# for row in csv.reader(open(repbase), delimiter = '\t'):

	# 	black_list.add(row[9])

	forward = "Rd1"

	pair_ori = 0
	if forward == "Rd1":
		pair_ori = 1
	elif forward == "Rd2":
		pair_ori = -1

	for row in csv.reader(open(reads_genome), delimiter = '\t'):


		if row[0][0] != "@":
			#if "N" in row[5]:
			read = row[0]
			flag = row[1]
			chr = row[2]
			start = int(row[3]) - 1             #Sam es 1 referenciado y es mas comodo trabajar en cordeneadas 0 refereciadas
			cigar = row[5]
			seq = row[9]

			pair_ori = 0
			if forward == "Rd1":
				pair_ori = 1
			elif forward == "Rd2":
				pair_ori = -1

			self_strand = 1
			pair_strand = '+'

			#Si no se tiene el flag XS:A:- se tienen que implementar las operaciones a nivel de bits:

			if (1 & int(flag)):    #paired end
				pair_number = flag_dict[flag][0]
				self_strand = flag_dict[flag][1]
				if pair_ori*self_strand*pair_number==-1:
					pair_strand = '-'

			elif (16 & int(flag)):   #single end
				self_strand = -1
				pair_strand = '-'

			if self_strand == -1:
				seq = str(Seq(seq).reverse_complement())

			if flag == "0" or flag == "16":

				aux_str = ''
				cigar_vars = []

				strand = "+"

				if flag == "16":
					strand = "-"

				for c in cigar:
					c_type = ascii_classifier(c)

					if c_type == 'number':
						aux_str += c

					elif c_type == 'letter':
						cigar_vars.append((c, int(aux_str)))
						aux_str = ''

				Exon_starts = [start]
				Exon_ends = []

				block = 0
				var_index = 0

				score = 0

				for var in cigar_vars:
					var_type = var[0]
					var_value = var[1]
					var_index += 1

					if var_type == 'M':
						block += var_value
						score += var_value

					if var_type == 'D':
						block += var_value
						score -= 1

					if var_type == 'I':
						block += 0
						score -= var_value

					if var_type == 'N':
						Exon_ends.append(Exon_starts[-1] + block)
						Exon_starts.append(Exon_ends[-1] + var_value)
						block = 0

					if var_index == len(cigar_vars):
						Exon_ends.append(Exon_starts[-1] + block)

				exons_reads_genome[read].append((chr, strand, start, cigar, score, Exon_starts, Exon_ends))

				# if len(Exon_ends)>=2:

				# 	print chr, start, cigar, matches, Exon_starts, Exon_ends



	for row in csv.reader(open(row_ME), delimiter = '\t'):


		read, flag, tag, start, cigar, seq, qual, q_block_starts, q_block_ends,  micro_exon_seq_found, I_pos_tag, DRU, DRD, DR_corrected_micro_exon_seq_found, micro_exons_coords = row
		intron_tag, transcript_ID, anchors = tag.split("|")
		chr, istart, iend = re.findall(r"[\w']+", intron_tag)
		up, down = anchors.split("_")
		up = int(up)
		down = int(down)

		start = int(start)


		t_score = 0

		aux_str = ''
		cigar_vars = []

		for c in cigar:
			c_type = ascii_classifier(c)

			if c_type == 'number':
				aux_str += c

			elif c_type == 'letter':
				cigar_vars.append((c, int(aux_str)))
				aux_str = ''


		var_index = 0

		for var in cigar_vars:
			var_type = var[0]
			var_value = var[1]
			var_index += 1

			if var_type == 'M':
				t_score += var_value

			elif var_type == 'I':
				t_score += var_value

			elif var_type == 'D':
				t_score -= 1


		micro_exons = []

		for i in micro_exons_coords.split(","):
			micro_exons.append((i.split("_")))

		U2_scores = []
		TOTAL_mean_conservation = []




		for i in micro_exons:



			ME_strand, ME_start, ME_end = i[-3:]
			ME_chr = "_".join(i[:-3])

			ME_start = int(ME_start)
			ME_end = int(ME_end)



			mean_conservation= phylop_bw.stats(ME_chr, ME_start-2, ME_end+2, type="mean")[0]




			if mean_conservation==None:

				mean_conservation=0





			TOTAL_mean_conservation.append(mean_conservation)



			ME5 = str(Genome[ME_chr][ME_start-14:ME_start+3]).upper()
			ME3 = str(Genome[ME_chr][ME_end-3:ME_end+10]).upper()


			if ME_strand == "-":

				ME5 = str(Genome[ME_chr][ME_end-3:ME_end+14].reverse_complement()).upper()
				ME3 = str(Genome[ME_chr][ME_start-10:ME_start+3].reverse_complement()).upper()




			U2_score = 0

			i = 0

			for N in ME5:
				U2_score += U2_GTAG_3[N][i]
				i += 1

			i = 0

			for N in ME3:
				U2_score += U2_GTAG_5[N][i]
				i += 1

			U2_score = percent(U2_score, TOTAL_U2_max_score)

			U2_scores.append(U2_score)


		tag_alingment = "|".join([tag, cigar])
		genome_alingment = "None"
		g_score = "None"

		same_ME = False

		if start <= (up-8) and cigar.count("I") == 1 and cigar.count("D") == 0 and cigar.count("S") == 0 and (read in black_list)==False and len(DR_corrected_micro_exon_seq_found) <= ME_len:


			if read in exons_reads_genome:

				g_chr, g_strand, g_start, g_cigar, g_score, g_Exon_starts, g_Exon_ends = max(exons_reads_genome[read],key=lambda item:item[4])

				genome_alingment = "|".join([g_chr, g_strand, str(g_start), g_cigar])

				for s, e in zip(g_Exon_starts, g_Exon_ends):

					if "_".join(map(str, [g_chr, g_strand, s, e])) in micro_exons_coords:

						same_ME = True
#

				if same_ME:
					print read, seq, qual, tag_alingment, t_score, genome_alingment, g_score, same_ME, len(DR_corrected_micro_exon_seq_found), DR_corrected_micro_exon_seq_found, len(micro_exons), max(U2_scores), max(TOTAL_mean_conservation), micro_exons_coords, ",".join(map(str, U2_scores)), ",".join(map(str, TOTAL_mean_conservation))


				elif t_score > g_score:
					print read, seq, qual, tag_alingment, t_score, genome_alingment, g_score, same_ME, len(DR_corrected_micro_exon_seq_found), DR_corrected_micro_exon_seq_found, len(micro_exons), max(U2_scores), max(TOTAL_mean_conservation), micro_exons_coords, ",".join(map(str, U2_scores)), ",".join(map(str, TOTAL_mean_conservation))

			else:
				print read, seq, qual, tag_alingment, t_score, genome_alingment, g_score, same_ME, len(DR_corrected_micro_exon_seq_found), DR_corrected_micro_exon_seq_found, len(micro_exons), max(U2_scores), max(TOTAL_mean_conservation), micro_exons_coords, ",".join(map(str, U2_scores)), ",".join(map(str, TOTAL_mean_conservation))


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	#main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9] ) # Old
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]))

#python ~/my_src/ME/Pipeline/ME_filter1.py _clip10.trim.sam.row_ME _clip10.trim.sam.row_ME.hg19.sam _clip10.trim.sam.row_ME.fastq.dust _clip10.trim.sam.row_ME.fastq.Repbase

#python ~/my_src/ME/Pipeline/ME_filter1.py _clip1.trim.sam.row_ME _clip1.trim.sam.row_ME.hg19.sam _clip1.trim.sam.row_ME.fastq.dust _clip1.trim.sam.row_ME.fastq.Repbase

#python ~/my_src/ME/Pipeline/ME_filter1.py _clip1.trim.sam.row_ME _clip1.trim.sam.row_ME.hg19.sam _clip1.trim.sam.row_ME.fastq.dust _clip1.trim.sam.row_ME.fastq.Repbase ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix

#python ~/my_src/ME/Pipeline/ME_filter1.py ~/db/genome/hg19.fa _clip1.trim.sam.row_ME _clip1.trim.sam.row_ME.hg19.sam _clip1.trim.sam.row_ME.fastq.dust _clip1.trim.sam.row_ME.fastq.Repbase ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix

#python ~/my_src/ME/Pipeline/ME_filter1.py ~/db/genome/hg19.fa _clip1.trim.sam.row_ME _clip1.trim.sam.row_ME.hg19.sam _clip1.trim.sam.row_ME.fastq.dust _clip1.trim.sam.row_ME.fastq.Repbase ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw

#python ~/my_src/ME/Pipeline/ME_filter1.py ~/db/genome/hg19.fa _clip1.trim.sam.row_ME _clip1.trim.sam.row_ME.hg19.sam _clip1.trim.sam.row_ME.fastq.dust _clip1.trim.sam.row_ME.fastq.Repbase ~/db/PWM/hg19_GT_AG_U2_5.good.matrix ~/db/PWM/hg19_GT_AG_U2_3.good.matrix ~/db/hg19.100way.phyloP100way.bw ~/db/hg19.46way.phyloP46way.primates.bw

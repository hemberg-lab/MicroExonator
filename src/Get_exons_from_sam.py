import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict


SeqTable = []

#dicionario de fag [Rd1/Rd2, +/-]
flag_dict = {'73':[1,1], '89':[1,1], '121':[1,-1], '153':[-1,-1], '185':[-1,-1], '137':[-1,1], '99':[1,1], '147':[-1,-1], '83':[1,-1], '163':[-1,1], '67':[1,1], '115':[1,-1], '179':[-1,-1], '81':[1,-1], "161":[-1,1], '97':[1,1], '145':[-1,-1], '65':[1,1], '129':[-1,1], '113':[1,-1], '177':[-1,-1] }



def ascii_classifier(c): 

	ascii = ord(c)
	if ascii >= 48 and ascii <= 57:
		return 'number'
	else:
		return 'letter'
			

def main(sam):          #hay que indicar si forward es Rd1 o Rd2
	reader = csv.reader(open(sam), delimiter = '\t')

	exon_count = defaultdict(int)
	
	forward = "Rd1"
	
	pair_ori = 0
	if forward == "Rd1":
		pair_ori = 1
	elif forward == "Rd2":
		pair_ori = -1
	
	
	for row in reader:
		if row[0][0] != "@":
			if "N" in row[5]:
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
						
				aux_str = ''
				cigar_vars = [] 
							
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
						
				for var in cigar_vars:
					var_type = var[0]
					var_value = var[1]
					var_index += 1
					
					if var_type == 'M':
						block += var_value						
						
					if var_type == 'D':
						block += var_value
											
					if var_type == 'I':
						block += 0
						
					if var_type == 'N':
						Exon_ends.append(Exon_starts[-1] + block)
						Exon_starts.append(Exon_ends[-1] + var_value)
						block = 0
						
					if var_index == len(cigar_vars):
						Exon_ends.append(Exon_starts[-1] + block)

				if len(Exon_starts)>=3:

				
					#for e5s, e5e, e3s, e3e  in zip(Exon_starts, Exon_ends, Exon_starts[1:], Exon_ends[1:]):
					for estart, eend in zip(Exon_starts[1:-1], Exon_ends[1:-1]):

						exon = "_".join(map(str, [chr, estart, eend]))

						exon_count[exon] += 1


	for exon, count in 	exon_count.items():

		print exon, count
							
			
					
			



if __name__ == '__main__':
	main(sys.argv[1])

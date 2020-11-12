import sys
import csv
from collections import defaultdict

csv.field_size_limit(100000000)

def main(pre_processed, genome_sam):

	read_SJ = defaultdict(set)
	black_list = set([])

	# for row in csv.reader(open(dust), delimiter = '>'):

	# 	black_list.add(row[1])

	# for row in csv.reader(open(repbase), delimiter = '\t'):

	# 	black_list.add(row[9])

	for row in csv.reader(open(genome_sam), delimiter = '\t'):
                
            try:
	        if row[1]=="0" or row[1]=="16":
	            black_list.add(row[0])
            except ValueError:
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
		if (read in black_list)==False:
		    print "\t".join(row)
            except ValueError:
                pass
		#print black_list


main(sys.argv[1], sys.argv[2]) #, sys.argv[3], sys.argv[4])

import sys
import csv


def main(alingment_pre_processed_round2):
	
	for row in csv.reader(open(alingment_pre_processed_round2), delimiter = '\t'):

		read, flag, tag, start, cigar, seq, qual = row

		print "@" + read
		print seq
		print "+"
		print qual


if __name__ == '__main__':
	main(sys.argv[1])

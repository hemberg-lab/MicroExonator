import sys
import csv


def main(alingment_pre_processed_round2):
	
	for row in csv.reader(open(alingment_pre_processed_round2), delimiter = '\t'):

		try:
			read, flag, tag, start, cigar, seq, qual = row

			if len(seq)>len(qual):
				qual = qual + qual[ -(len(seq) - len(qual)) : ]

			elif len(seq)<len(qual):
				qual = qual[:len(seq)]

			print "@" + read
			print seq
			print "+"
			print qual

		except ValueError: #minor fraction of lines dont have 7 columns due to errors.
			pass


if __name__ == '__main__':
	main(sys.argv[1])

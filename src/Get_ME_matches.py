
import sys
import csv


def main():

    for row in csv.reader(open(ME_centric), delimiter = '\t'):

        ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_scores, mean_conservations, P_MEs, total_ME = row

		ME_strand, ME_start, ME_end = ME.split("_")[-3:]
		ME_chrom =  "_".join(ME.split("_")[:-3])
    
    
if __name__ == '__main__':
        Genomictabulator(sys.argv[1])
        main (sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]), sys.argv[8])

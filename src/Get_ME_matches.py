
import sys
import csv


def main(ME_centric):

    for row in csv.reader(open(ME_centric), delimiter = '\t'):

        ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_scores, mean_conservations, P_MEs, total_ME = row

		ME_strand, ME_start, ME_end = ME.split("_")[-3:]
		ME_chrom =  "_".join(ME.split("_")[:-3])
		
		total_ME.split(",")
    
    
if __name__ == '__main__':
        main (sys.argv[1])

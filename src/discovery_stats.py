
import sys
import csv


##### Use as #####

    #python bed12 MicroExonator_high_qual ME_len

    #Compatible with python2 and python3

######


def main(bed12, ME_high_cual , ME_len):

    ME_annotated = set([])
    ME_MicroExonator = set([])

    with csv.reader(open(bed12) as  db:

        reader = csv.reader( high_qual_filters , delimiter = '\t') )

		for row in csv.reader(open(bed12), delimiter = '\t'):


			blocksizes = list(map(int, row[10].strip(",").split(",")))
			qstarts = list(map(int, row[11].strip(",").split(",")))


			start = int(row[1])
			end = int(row[2])
			strand = row[5]
			bn = int(row[9])
			chrom = row[0]


			for q1, q2, q3, b1, b2, b3 in zip(qstarts, qstarts[1:] , qstarts[2:], blocksizes, blocksizes[1:], blocksizes[2:]):


				estart = start + q2
				eend = start + q2 + b2
				elength = eend - estart
				exon = (chrom, strand, str(estart), str(eend))
                transcript = row[3]
				SJ_start = start + q1 + b1
				SJ_end = start + q3
				ME_intron = " ".join([chrom, str(SJ_start), str(SJ_end), "SJ", "0", strand])

				if elength <= ME_len:

                    ME_annotated.add(exon)


    with csv.reader(open(ME_high_cual) as  high_qual_filters:

        reader = csv.DictReader( high_qual_filters , delimiter = '\t') )

        for row in reader:

            # Hearder is "ME", "transcript", "sum_total_coverage", "total_SJs", "total_coverages", "len_micro_exon_seq_found", "micro_exon_seq_found", "total_number_of_micro_exons_matches", "U2_scores",  "mean_conservations_vertebrates", "P_MEs", "total_ME",   "ME_P_value", "ME_type"

            ME = reader["ME"]

            chrom = "_".join(ME.split("_")[:-3])  # To avoid errors with chrom_
            starnd, start, end = ME.split("_")[-3:]

            ME = (chrom, strand, start, end)

            ME_MicroExonator.add(ME)


    ######## TO DO #######

    # Operate sets of microexos found at annotation (bed12) -> stored as a set in ME_annotated
    # and ME found by MicroExonator -> stored as a set in ME_MicroExonator

    ## Hint lock for set operation in python

    ### https://docs.python.org/2/library/sets.html


    ### generate percent of overlap, and do plots (maybe ven diagram)

    ######################

if __name__ == '__main__':
	main (sys.argv[1], sys.argv[2], int(sys.argv[3]) )

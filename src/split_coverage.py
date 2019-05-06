import sys
import csv

def main(ME_SJ_coverage):

    with open(ME_SJ_coverage) as F:

        reader = csv.reader(F, delimiter="\t")

        #print( "\t".join( [ "ME", "Coord", "PSI", "CI_Lo", "CI_Hi", "is_alternative_5", "alternatives_5",  "is_alternative_3", "alternatives_3" ]) )

        for row in reader:

            FILE, ME, total_SJs, ME_SJ_coverages, sum_ME_coverage, sum_ME_SJ_coverage_up_down_uniq, sum_ME_SJ_coverage_up, sum_ME_SJ_coverage_down, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3,  total_cov_alternatives_3 = row
            
            with open( FILE + ".sam.pre_processed.filter1.ME_SJ_coverage", "a") as out:
                out.write( "\t".join(row) + "\n")

if __name__ == '__main__':
	main(sys.argv[1]  )

from snakemake.utils import min_version
import csv
from collections import defaultdict
import sys




#with open(snakemake.output[0]) as out
with open("test.txt") as out:

    min_read_per_sample = int(snakemake.params[0])
    ME_n_samples = defaultdict(int)

    #for file in snakemake.input:
    for file in sys.argv[1:]:

        with open(file) as cov_file:

            reader = csv.reader(file)

            for row in reader:

                FILE_NAME, ME, total_SJs, ME_SJ_coverages, sum_ME_coverage, sum_ME_SJ_coverage_up_down_uniq, sum_ME_SJ_coverage_up, sum_ME_SJ_coverage_down, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3, total_cov_alternatives_3 = row

                if int(sum_ME_SJ_coverage_up_down_uniq)>=min_read_per_sample:
                    ME_n_samples[ME] += 1
                    
    out.write("\t".join(["ME", "N_samples" ] + "\n")
    
    for ME, n in ME_n_samples.items():
        
        out.write("\t".join([ME, str(n)] + "\n")
                    
                    
#if __name__ == '__main__':
#	main(sys.argv[1:])
                

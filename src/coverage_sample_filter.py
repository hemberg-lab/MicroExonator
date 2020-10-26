from snakemake.utils import min_version
import csv
from collections import defaultdict
import sys




with open(snakemake.output[0]) as out
#with open("test.txt", "w") as out:

    min_read_per_sample = int(snakemake.params[0])
    #min_read_per_sample = 5
    ME_n_samples = defaultdict(int)
    all_ME = set([])

    #for file in snakemake.input:
    for file in sys.argv[1:]:

        with open(file) as cov_file:

            reader = csv.reader(cov_file, delimiter="\t")
           

            for row in reader:


                FILE_NAME, ME, total_SJs, ME_SJ_coverages, sum_ME_coverage, sum_ME_SJ_coverage_up_down_uniq, sum_ME_SJ_coverage_up, sum_ME_SJ_coverage_down, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3, total_cov_alternatives_3 = row

                all_ME.add(ME) 
                if int(sum_ME_SJ_coverage_up_down_uniq)>=min_read_per_sample:
                    ME_n_samples[ME] += 1
                    
    out.write("\t".join(["ME", "N_samples" ]) + "\n")
    
    for ME in all_ME:
        
        n = ME_n_samples[ME]
        
        out.write("\t".join([ME, str(n)]) + "\n")

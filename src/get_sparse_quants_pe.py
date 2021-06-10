import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import random

#This script filter the correted quant files to get sparce quant files


paired_sum = defaultdict(int)
with open( gzip.open(snakemake.output["corrected_sparse"], "wt"), "wt") as out, gzip.open(snakemake.input["corrected_quant_rd1"], "rt") as rd1, gzip.open(snakemake.input["corrected_quant_rd2"], "rt") as rd2:

    header =  "\t".join(["sample", "ME", "ME_coverages", "excluding_covs",  "PSI", "CI_Lo", "CI_Hi"])

    out.write(header + "\n")  
    
    paired_sum = defautdict(int)

    reader1 = csv.reader(rd1, delimiter="\t")
    reader2 = csv.reader(rd2, delimiter="\t")

    for row in reader1:

        sample, ME, ME_coverages, excluding_covs = row
        ME_coverages = float(ME_coverages)
        excluding_covs = float(excluding_covs)
        pair_index = paired_dic[sample]
        paired_sum[(pair_index, ME, "ME_coverages")] += ME_coverages
        paired_sum[(pair_index, ME, "excluding_covs")] += excluding_covs  
            
    for row in reader2:

        sample, ME, ME_coverages, excluding_covs = row
        ME_coverages = float(ME_coverages)
        excluding_covs = float(excluding_covs)
        pair_index = paired_dic[sample]
        paired_sum[(pair_index, ME, "ME_coverages")] += ME_coverages
        paired_sum[(pair_index, ME, "excluding_covs")] += excluding_covs
        
        
    for info, count in paired_sum.items():

        sample, ME, count_type = info

        if count_type=="ME_coverages":

            ME_coverages = count
            excluding_covs = paired_sum[(sample, ME, "excluding_covs")]

            if ME_coverages + excluding_covs >=5:

                PSI = ME_coverages/(excluding_covs+ME_coverages)
                CI_Lo, CI_Hi = calcBin(ME_coverages,  ME_coverages+excluding_covs)


                outrow = "\t".join(map(str, [sample, ME, ME_coverages, excluding_covs,  PSI, CI_Lo, CI_Hi], ))

                out.write(outrow + "\n")

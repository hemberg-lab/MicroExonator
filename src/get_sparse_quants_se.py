import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import random

#This script filter the correted quant files to get sparce quant files


paired_sum = defaultdict(int)
with open( "/home/jovyan/Hannes/TOTAL.corrected.light.psi.hannes.2.tsv" , "wt") as out:

    header =  "\t".join(["sample", "ME", "ME_coverages", "excluding_covs",  "PSI", "CI_Lo", "CI_Hi"])

    out.write(header + "\n")    

    for f in filenames:
        
        sample = f.split("/")[-1].split(".")[0]
        
        
        if sample not in single_cell and sample not in completed:
            
            print(sample)
        
            if sample in paired_dic:

                with gzip.open(f, "rt") as file:

                    reader = csv.reader(file, delimiter="\t")

                    for row in reader:

                        sample, ME, ME_coverages, excluding_covs = row
                        ME_coverages = float(ME_coverages)
                        excluding_covs = float(excluding_covs)

                        pair_index = paired_dic[sample]

                        paired_sum[(pair_index, ME, "ME_coverages")] += ME_coverages
                        paired_sum[(pair_index, ME, "excluding_covs")] += excluding_covs

            else:
                with gzip.open(f, "rt") as file:

                    reader = csv.reader(file, delimiter="\t")

                    for row in reader:

                        sample, ME, ME_coverages, excluding_covs = row

                        ME_coverages = float(ME_coverages)
                        excluding_covs = float(excluding_covs)

                        if ME_coverages + excluding_covs >=5:

                            PSI = ME_coverages/(excluding_covs+ME_coverages)
                            CI_Lo, CI_Hi = calcBin(ME_coverages,  ME_coverages+excluding_covs)


                            outrow = "\t".join(map(str, [sample, ME, ME_coverages, excluding_covs,  PSI, CI_Lo, CI_Hi], ))

                            out.write(outrow + "\n")

import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import random

#This script filter the correted quant files to get sparce quant files


paired_sum = defaultdict(int)
with gzip.open( snakemake.input["corrected_quant"] , "rt") as f, gzip.open(snakemake.output["corrected_sparse"] , "wt") as out:

    header =  "\t".join(["sample", "ME", "ME_coverages", "excluding_covs",  "PSI", "CI_Lo", "CI_Hi"])

    out.write(header + "\n")    
    reader = csv.reader(f, delimiter="\t")

    for row in reader:

        sample, ME, ME_coverages, excluding_covs = row

        ME_coverages = float(ME_coverages)
        excluding_covs = float(excluding_covs)

        if ME_coverages + excluding_covs >=5:

            PSI = ME_coverages/(excluding_covs+ME_coverages)
            CI_Lo, CI_Hi = calcBin(ME_coverages,  ME_coverages+excluding_covs)

            outrow = "\t".join(map(str, [sample, ME, ME_coverages, excluding_covs,  PSI, CI_Lo, CI_Hi], ))
            out.write(outrow + "\n")

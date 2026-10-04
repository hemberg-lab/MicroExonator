import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import random

#This script filter the correted quant files to get sparce quant files for pseudopools


def binP(N, p, x1, x2):
    p = float(p)
    q = p/(1-p)
    k = 0.0
    v = 1.0
    s = 0.0
    tot = 0.0

    while(k<=N):
            tot += v
            if(k >= x1 and k <= x2):
                    s += v
            if(tot > 10**30):
                    s = s/10**30
                    tot = tot/10**30
                    v = v/10**30
            k += 1
            v = v*q*(N+1-k)/k
    return s/tot

def calcBin(vx, vN, vCL = 95):
    '''
    Calculate the exact confidence interval for a binomial proportion
    Usage:
    >>> calcBin(13,100)
    (0.07107391357421874, 0.21204372406005856)
    >>> calcBin(4,7)
    (0.18405151367187494, 0.9010086059570312)
    '''
    vx = float(vx)
    vN = float(vN)
    #Set the confidence bounds
    vTU = (100 - float(vCL))/2
    vTL = vTU

    vP = vx/vN
    if(vx==0):
            dl = 0.0
    else:
            v = vP/2
            vsL = 0
            vsH = vP
            p = vTL/100

            while((vsH-vsL) > 10**-5):
                    if(binP(vN, v, vx, vN) > p):
                            vsH = v
                            v = (vsL+v)/2
                    else:
                            vsL = v
                            v = (v+vsH)/2
            dl = v

    if(vx==vN):
            ul = 1.0
    else:
            v = (1+vP)/2
            vsL =vP
            vsH = 1
            p = vTU/100
            while((vsH-vsL) > 10**-5):
                    if(binP(vN, v, 0, vx) < p):
                            vsH = v
                            v = (vsL+v)/2
                    else:
                            vsL = v
                            v = (v+vsH)/2
            ul = v
    return (dl, ul)
  
  
  
# One pseudo-bulk per output file: Report/quant/corrected/PSI_sparse/single_cell/<cell type>-<n>.corrected.PSI.gz
pseudo_pool_ID = snakemake.output["corrected_sparse"].split("/")[-1]
if pseudo_pool_ID.endswith(".corrected.PSI.gz"):
    pseudo_pool_ID = pseudo_pool_ID[:-len(".corrected.PSI.gz")]
cell_type = pseudo_pool_ID.rsplit("-", 1)[0].replace("_", " ")

# PSI is reported only with at least min_reads_PSI corrected reads, as in the
# bulk files; a microexon with exclusion reads only is kept with PSI 0.
min_reads = float(snakemake.params["min_reads"])

pool_ME_coverages = defaultdict(float)
pool_excluding_covs = defaultdict(float)

for f in snakemake.input["cells"]:

    with gzip.open(f, "rt") as file:

        reader = csv.reader(file, delimiter="\t")

        for row in reader:

            sample, ME, ME_coverages, excluding_covs = row
            pool_ME_coverages[ME] += float(ME_coverages)
            pool_excluding_covs[ME] += float(excluding_covs)


with gzip.open(snakemake.output["corrected_sparse"], "wt") as out:

    header = "\t".join(["ME", "pseudo_pool", "cell_type",  "ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi"])

    out.write(header + "\n")

    for ME in sorted(pool_ME_coverages):

        ME_coverages = pool_ME_coverages[ME]
        excluding_covs = pool_excluding_covs[ME]
        total = ME_coverages + excluding_covs

        if total >= min_reads and total > 0:

            CI_Lo, CI_Hi = calcBin(ME_coverages, total)
            PSI = ME_coverages / total

            out_line = "\t".join(map(str, [ME, pseudo_pool_ID, cell_type,  ME_coverages, excluding_covs, PSI, CI_Lo, CI_Hi]))
            out.write(out_line + "\n")

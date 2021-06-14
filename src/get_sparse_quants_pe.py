import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import random

#This script filter the correted quant files to get sparce quant files

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

paired_dic = dict()

with open(snakemake.config["paired_samples"]) as file:
    
    reader = csv.reader(file, delimiter="\t")
    
    for row in reader:
        
        paired_dic[row[0]] = row[1]
        paired_dic[row[1]] = row[1]
        

paired_sum = defaultdict(int)

with gzip.open(snakemake.output["corrected_sparse"], "wt") as out, gzip.open(snakemake.input["corrected_quant_rd1"], "rt") as rd1, gzip.open(snakemake.input["corrected_quant_rd2"][0], "rt") as rd2:

    header =  "\t".join(["sample", "ME", "ME_coverages", "excluding_covs",  "PSI", "CI_Lo", "CI_Hi"])

    out.write(header + "\n")  
    
    paired_sum = defaultdict(int)

    reader1 = csv.reader(rd1, delimiter="\t")
    #reader2 = csv.reader(rd2, delimiter="\t")

    for row in reader1:

        sample, ME, ME_coverages, excluding_covs = row
        ME_coverages = float(ME_coverages)
        excluding_covs = float(excluding_covs)
        pair_index = paired_dic[sample]
        paired_sum[(pair_index, ME, "ME_coverages")] += ME_coverages
        paired_sum[(pair_index, ME, "excluding_covs")] += excluding_covs  
    
    reader2 = csv.reader(rd2, delimiter="\t")
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

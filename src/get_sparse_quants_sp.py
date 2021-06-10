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
  
  
  
pool_ME_coverages = defaultdict(int)
pool_excluding_covs = defaultdict(int)



for f in open(snakemake.input["cells"]):
    
    sample = f.split("/")[-1].split(".")[0]
    
    if sample in pseudo_pool_dict:

        with gzip.open(f, "rt") as file:

            reader = csv.reader(file, delimiter="\t")
            
            for row in reader:
                
                sample, ME, ME_coverages, excluding_covs = row       
                pseudo_pool_ID = pseudo_pool_dict[sample]

                ME = ME
                
                if float(ME_coverages) > 0:
                    pool_ME_coverages[(pseudo_pool_ID, ME)] += float(ME_coverages)
                    
                if float(excluding_covs) > 0:    
                    pool_excluding_covs[(pseudo_pool_ID, ME)] += float(excluding_covs)
                    
                    
pseudo_pool_set = set(pool_ME_coverages.keys())                   
                    
with gzip.open(snakemake.input["corrected_sparse"], "wt") as out:
  
    header = "\t".join(["ME", "pseudo_pool", "cell_type",  "ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi"])
    
    out.write(header + "\n")
    
    for key in pseudo_pool_set:

        if len(key)==2:
            
            pseudo_pool, ME = key

            cell_type = pseudo_pool.split("-")[0].replace("_", " ")
            ME_coverages = pool_ME_coverages[key]
            excluding_covs = pool_excluding_covs[key]
            CI_Lo, CI_Hi = calcBin(float(ME_coverages),  ME_coverages+excluding_covs)
            PSI = ME_coverages/(ME_coverages+excluding_covs)

            out_line = "\t".join(map(str, [ME, pseudo_pool, cell_type,  ME_coverages, excluding_covs, PSI, CI_Lo, CI_Hi]))
            out.write(out_line + "\n")

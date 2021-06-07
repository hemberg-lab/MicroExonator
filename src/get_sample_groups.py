import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import glob
import itertools
import csv
import gzip
from collections import defaultdict
import random


#This script gets the sample groups (output 1) and compute quantitative filters to get a list of reliable detected microexons.

csv.field_size_limit(100000000)

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
  
def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
  

primary_clusters = defaultdict(list)



        
        
with open(snakemake.params["run_metadata_sc"]) as file:
  
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
        
        primary_clusters[snakemake.params["cell_type"]].append(snakemake.params["cells"])



def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
    

    
    
sample_group = dict()
pseudo_pool_dict =  dict()

for cluster in primary_clusters:
    
    c = 0
    cells = primary_clusters[cluster]
    pseudo_pools = partition( cells , 3 )
    pseudo_pool_ID = ""
    
    for pool in pseudo_pools:
        
        c+=1
        pseudo_pool_ID =cluster.replace(" ", "_") + "-" + str(c)
        
        
        for cell in pool:
            pseudo_pool_dict[cell] = pseudo_pool_ID
            sample_group[cell] = pseudo_pool_ID  

 with open(snakemake.params["bulk_samples"]) as file:
    
    reader = csv.DicReader(file, delimiter="\t")
    
    for row in reader:
        
        sample_group[row["sample"]] = row["condition"]  
        

with open(snakemake.params["pseudo_pool_membership"]) as out_pseudo_pool_membership, open(snakemake.params["sample_groups"]) as out_sample_groups:
    
    for cell, sp in pseudo_pool_dict.itemps():
        
        out = "\t".join([cell, sp])
        
        out_pseudo_pool_membership.write(out + "\n")
        
 
    
    for sample, group in sample_group.itemps():
        
        out = "\t".join([sample, group])
        
        sample_group.write(out + "\n")
        
    
    
# input : 
# params : run_metadata_sc, cell_type, cells
# out : pseudo_pool_membership, sample_groups

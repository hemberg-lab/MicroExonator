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


  
def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
  

primary_clusters = defaultdict(list)

        
with open(config["cluster_metadata"]) as file:
  
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
        primary_clusters[config["cluster_name"]].append(config["file_basename"])



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
            pseudo_pool_dict[pseudo_pool_ID] = cell
            sample_group[pseudo_pool_ID] = cell  
           
 with open(config["bulk_samples"]) as file:
    
    reader = csv.DicReader(file, delimiter="\t")
    for row in reader:
        sample_group[row["condition"]] = row["sample"]  
        

with open("pseudo_pool_membership.txt", "w") as out_pseudo_pool_membership, open("sample_groups.txt", "w") as out_sample_groups:
    
    for sp, cell  in pseudo_pool_dict.items():
        
        out = "\t".join([cell, sp])  
        out_pseudo_pool_membership.write(out + "\n")
        
    
    for  group, sample in sample_group.items():
        
        out = "\t".join([sample, group])
        sample_group.write(out + "\n")
        
    
    
rule quant_pseudo_pools:
    input:
        cells = lambda w: pseudo_pool_dict[w.cluster],
    output:
        "Report/quant/corrected/pseudo_pools/{cluster}.corrected.sparce.psi",
    priority: 10
    script:
        "../src/get_sparse_quants_sp.py"  
    
# input : 
# params : run_metadata_sc, cell_type, cells
# out : pseudo_pool_membership, sample_groups

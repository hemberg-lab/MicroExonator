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

pe_samples = set([])

if "paired_samples" in config:
    
    with  open(config["paired_samples"]) as file:
        
        reader = reader(file, delimiter="\t")
        for row in reader:
            
            pe_samples.add(row[0])
            pe_samples.add(row[1])


with open(config["cluster_metadata"]) as file:
  
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
        primary_clusters[config["cluster_name"]].append(config["file_basename"])



def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
    

sample_group_se = dict()
sample_group_pe = dict()
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
           
 with open(config["bulk_samples"]) as file:
    
    reader = csv.DicReader(file, delimiter="\t")
    for row in reader:
        
        if row["sample"] in pe_samples:  
            sample_group_se[row["condition"]] = row["sample"]
        else:
            sample_group_pe[row["condition"]] = row["sample"]  
        

with open("pseudo_pool_membership.txt", "w") as out_pseudo_pool_membership, open("sample_groups.txt", "w") as out_sample_groups:
    
    for sp, cell  in pseudo_pool_dict.items():
        
        out = "\t".join([cell, sp])  
        out_pseudo_pool_membership.write(out + "\n")
        
    
    for  group, sample in sample_group_se.items():
        
        out = "\t".join([sample, group])
        sample_group.write(out + "\n")
        
    for  group, sample in sample_group_pe.items():
        
        out = "\t".join([sample, group])
        sample_group.write(out + "\n")
        
        
def get_cell_sp(cluster):
    return(expand( "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz", sample = pseudo_pool_dict[cluster])
    
    
rule get_sparse_quants_sp:
    input:
        cells = get_cell_sp(w.cluster)
    output:
        corrected_sparse = "Report/quant/sparse/single_cell/{cluster}.corrected.PSI.gz"
    priority: 10
    script:
        "../src/get_sparse_quants_sp.py" 
          
        
rule get_sparse_quants_se:
    input:
        corrected_quant = "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz"
    output:
        corrected_sparse = "Report/quant/sparse/blulk/{sample}.corrected.PSI.gz"
    priority: 10
    script:
        "../src/get_sparse_quants_se.py"
           
           
rule get_sparse_quants_pe:
    input:
        corrected_quant = "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz"
    output:
        corrected_sparse = "Report/quant/sparse/blulk/{sample}.corrected.PSI.gz"
    priority: 10
    script:
        "../src/get_sparse_quants_pe.py"
           
         
  rule detection_filter:
    input : 
        bulk_se = expand("Report/quant/sparse/blulk/{sample}.corrected.PSI.gz", sample = sample_group_se.keys()),
        bulk_pe = expand("Report/quant/sparse/blulk/{sample}.corrected.PSI.gz", sample = sample_group_pe.keys()),
        bulk_ME_reads_se = expand( "Round2/ME_reads/{sample}.counts.tsv",  sample = sample_group_se.keys()),
        bulk_ME_reads_pe = expand( "Round2/ME_reads/{sample}.counts.tsv",  sample = sample_group_pe.keys()),
        single_cell = expand("Report/quant/sparse/single_cell/{cluster}.corrected.PSI.gz", cluster = pseudo_pool_dict.keys()),
        single_cell_reads = expand( "Round2/ME_reads/{cluster}.counts.tsv",  cluster = pseudo_pool_dict.keys())
    out:
       detected_list = "Report/filter/detected_ME.txt"
    priority: 10
    script:
        "../src/get_sparse_quants_sp.py"  


import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version


#This script gets the sample groups (output 1) and compute quantitative filters to get a list of reliable detected microexons.

csv.field_size_limit(100000000)


  

  

primary_clusters = defaultdict(list)


with open("/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Runs/Paper/single_cell/MicroExonator/Tasic_clustering.txt") as file:
    
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
        
        primary_clusters[row["primary_type"]].append(row["Run_s"])
        
        
with open(snakemake.params["run_metadata_sc"]) as file:
  
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
        
        primary_clusters[snakemake.params["cell_type"]].append(snakemake.params["cells"])



def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
    
    
sample_group = 

with open(snakemake.params["run_metatda_sc"]) as file:
    
    reader = csv.DicReader(file, delimiter="\t")
    
    for row in reader:
        

 with open(snakemake.params["run_metatda_bulk"]) as file:
    
    reader = csv.DicReader(file, delimiter="\t")
    
    for row in reader:

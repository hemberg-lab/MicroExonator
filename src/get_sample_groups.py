import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version


#This script gets the sample groups (output 1) and compute quantitative filters to get a list of reliable detected microexons.

csv.field_size_limit(100000000)


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

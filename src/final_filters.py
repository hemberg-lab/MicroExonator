
import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version

#This gets a high quality list of microexons that were reliable discovered

se_dict = snakemake.params["se_dict"]
pe_dict = snakemake.params["pe_dict"]
pseudo_dict = snakemake.params["pseudo_dict"]

sample_group_CI_Lo = defaultdict(list)

for input in snakemake.input["bulk_se"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            sample_group = se_dict[row["sample"]]
            sample_group_CI_Lo[(sample_group, row["ME"])].append(row["CI_Lo"])

for input in snakemake.input["bulk_pe"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            sample_group = pe_dict[row["sample"]]
            sample_group_CI_Lo[(sample_group, row["ME"])].append(row["CI_Lo"])
             
for input in snakemake.input["single_cell"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            sample_group = pe_dict[row["sample"]]
            sample_group_CI_Lo[(sample_group, row["ME"])].append(row["CI_Lo"])
            
total_reads_files = snakemake.input["bulk_ME_reads_se"] + snakemake.input["bulk_ME_reads_pe"] + snakemake.input["single_cell_reads"]


spanning_reads_ME = defaultdict(int)

for input in total_reads_files:
    with open(input) as file:
        reader = csv.DictReader(file, delimiter="\t")
        
        for row in reader:
            spanning_reads_ME[row["ME"]] += int(row["Spanning_cov"]]
                                                

for key, value in sample_group_CI_Lo.items():

                                                
                                                

        
        

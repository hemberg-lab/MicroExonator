
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
            try:
                sample_group = se_dict[row["sample"]]
                sample_group_CI_Lo[(sample_group, row["ME"])].append(row["CI_Lo"])
            except KeyError:
                pass

for input in snakemake.input["bulk_pe"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            try:
                sample_group = pe_dict[row["sample"]]
                sample_group_CI_Lo[(sample_group, row["ME"])].append(row["CI_Lo"])
            except KeyError:
                pass
             
for input in snakemake.input["single_cell"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            try:
                sample_group = pe_dict[row["sample"]]
                sample_group_CI_Lo[(sample_group, row["ME"])].append(row["CI_Lo"])
            except KeyError:
                pass
            
total_reads_files = snakemake.input["bulk_ME_reads_se"] + snakemake.input["bulk_ME_reads_pe"] + snakemake.input["single_cell_reads"]


spanning_reads_ME = defaultdict(int)


for input in total_reads_files:
    with open(input) as file:
        reader = csv.DictReader(file, delimiter="\t")
        
        for row in reader:
            spanning_reads_ME[row["ME"]] += int(row["Spanning_cov"])
            


with open(snakemake.output["robustly_detected_ME"], "wt") as out:

    header =  "\t".join(["ME", "sample_group", "total_mesurements", "detected_samples"])
    out.write(header + "\n")
    
    for key, CI_lo_list in sample_group_CI_Lo.items():
        sample_group, ME = key

        total_mesurements = len(CI_lo_list)
        detected_samples = [x for x in CI_lo_list_len if x >= 0.1]

        if detected_samples/total_mesurements > 0.5:

            outrow = map(str, [ME, sample_group, total_mesurements, detected_samples])
            out.write(outrow + "\n")   
        

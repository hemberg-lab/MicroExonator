import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version


sample_group_CI_Lo = defaultdict(list)

for input in snakemake.input["files"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            sample_group_CI_Lo[row["ME"]].append(row["CI_Lo"])
            
    
spanning_reads_ME = defaultdict(int)

for input in total_reads_files:
    with open(input) as file:
        reader = csv.DictReader(file, delimiter="\t")
        
        for row in reader:
            spanning_reads_ME[row["ME"]] += int(row["Spanning_cov"])

with open(snakemake.output["detected"], "wt") as out:            
    for ME in sample_group_CI_Lo:

        CI_lo_list =  sample_group_CI_Lo[ME]

        total_mesurements = len(CI_lo_list)
        detected_samples = [x for x in CI_lo_list_len if x >= 0.1]

        if detected_samples/total_mesurements > 0.5:

            outrow = map(str, [ME, sample_group, total_mesurements, detected_samples])
            out.write(outrow + "\n")   
    

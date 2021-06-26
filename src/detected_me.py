import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version


sample_group_CI_Lo = defaultdict(list)

for input in snakemake.input["PSI_files"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            sample_group_CI_Lo[row["ME"]].append(float(row["CI_Lo"]))
            
    
spanning_reads_ME = defaultdict(int)

for input in snakemake.input["ME_reads"]:
    with open(input) as file:
        reader = csv.DictReader(file, delimiter="\t")
        
        for row in reader:
            spanning_reads_ME[row["ME"]] += int(row["Spanning_cov"])

with open(snakemake.output["detected"], "wt") as out:            
    header = "\t".join(['ME', 'total_mesurements', 'detected_samples', 'spanning_reads'])
    out.write(header + "\n")
    
    for ME in sample_group_CI_Lo:

        #sample_group = snakemake.params["sample_group"]

        CI_lo_list =  sample_group_CI_Lo[ME]

        total_mesurements = len(CI_lo_list)
        detected_samples = len([x for x in  CI_lo_list if x >= 0.1])

        if detected_samples/total_mesurements > 0.5 and spanning_reads_ME[ME]>=3:

            outrow = "\t".join(map(str, [ME, total_mesurements, detected_samples, spanning_reads_ME[ME]]))
            out.write(outrow + "\n")

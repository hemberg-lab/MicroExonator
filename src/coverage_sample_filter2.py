from snakemake.utils import min_version
import csv
from collections import defaultdict
import gzip

min_read_per_sample = int(snakemake.params[0])

with open(snakemake.output[0], 'w', newline='') as csvfile:

    csvfile.write("\t".join(["ME", "N_samples" ]) + "\n")
    
    files = [csv.DictReader(gzip.open(i,'rt'), delimiter="\t") for i in snakemake.input["PSI_files"] ]
    
    for rows in zip(*files):

        ME = ""
        detected_files = 0
        high_PSI = False

        for row in rows:
            
            ME = row["ME_coords"]
            
            if int(row["Unique_ME_reads"])>=min_read_per_sample:
                detected_files += 1
     
        csvfile.write("\t".join([ME, str(detected_files)]) + "\n")

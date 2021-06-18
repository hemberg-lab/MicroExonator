mport sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version

#This gets a high quality list of microexons that were reliable discovered

se_dict = snakemake.params["se_dict"]


ME_CI_Lo = defaultdict(list)

for input in snakemake.input["files"]:
    with gzip.open(input, "rt") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            ME_CI_Lo[row["ME"]].append(int(row["CI_Lo"]))
  


with open(snakemake.output["detected"], "wt") as out:

    header =  "\t".join(["ME", "total_mesurements", "detected_samples"])
    out.write(header + "\n")
    
    for ME, CI_lo_list in sample_group_CI_Lo.items():

        total_mesurements = len(CI_lo_list)
        detected_samples = [x for x in CI_lo_list_len if x >= 0.1]

        if detected_samples/total_mesurements > 0.5:

            outrow = map(str, [ME, total_mesurements, detected_samples])
            out.write(outrow + "\n")   

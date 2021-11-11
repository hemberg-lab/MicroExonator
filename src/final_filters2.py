import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version

#This gets a high quality list of microexons that were reliable discovered

min_detected_samples = int(snakemake.params["min_detected_samples"])


total_filter_files = snakemake.input["bulk_se"] + snakemake.input["bulk_pe"] + snakemake.input["single_cell"]


detected_ME = set([])

for f in total_filter_files:
    with open(f) as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:

            detected_samples = int(row['detected_samples'])

            if detected_samples>=min_detected_samples:

                detected_ME.add(row["ME"])

with open(snakemake.input["ME_centric"]) as file, open(snakemake.output["detected_list"], "w") as out:

    reader = csv.reader(file, delimiter="\t")

    header = ['ME', 'Transcript', 'Total_coverage', 'Total_SJs', 'ME_coverages', 'ME_length', 'ME_seq', 'ME_matches', 'U2_score', 'Mean_conservation', 'P_MEs', 'Total_ME']
    out.write("\t".join(header) + "\n")

    for row in reader:

        ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found , total_number_of_micro_exons_matches, U2_scores, mean_conservations_vertebrates, P_MEs, total_ME = row

        out_row = [ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, U2_scores, mean_conservations_vertebrates, P_MEs, total_ME]

        if ME in detected_ME:
            
            out.write("\t".join(out_row) + "\n")


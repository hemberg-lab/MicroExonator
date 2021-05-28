import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version

csv.field_size_limit(100000000)

tags_hit_dict = dict()

with open(snakemake.input["ME_centric"]) as file:
    
    reader = csv.reader(file, delimiter="\t")
    
    for row in reader:
        
        ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_scores, mean_conservations, P_MEs, total_ME = row
        
        for SJ in total_SJs.split(","):
            
            tags_hit_dict[(SJ, micro_exon_seq_found)] = ME
            

            
quant_ME_spanning_cov = defaultdict(int)

with open(snakemake.input["spanning_ME_reads"]) as file:
    
    reader = csv.reader(file, delimiter="\t")
    
    for row in reader:
        
        read, flag, tag, start, cigar, seq, qual = row
        
        SJ = tag.split("|")[0]
        micro_exon_seq_found = tag.split("|")[-1].split("_")[1]
        ME = tags_hit_dict[(SJ, micro_exon_seq_found)]
        
        quant_ME_spanning_cov[ME] +=1
        


with gzip.open(snakemake.input["quant"], "rt") as file, gzip.open(snakemake.output["corrected_quant"], "wt") as out: 
    
    reader = csv.DictReader(file, delimiter="\t")
    sample = snakemake.input["quant"].split("/")[-1].split(".")[0]
    
    for row in reader:
        
        ME = row["ME_coords"]
        
            
        full_counts =  quant_ME_spanning_cov[ME]


        ME_coverages = sum(map(float, row["ME_coverages"].split(",") ))
        excluding_covs = sum(map(float, row["SJ_coverages"].split(",") ))

        TOTAL_Alt5_3 = 0

        if row["Alt5"]!='None':
            excluding_covs += sum(map(float, row["Alt5_coverages"].split(",") ))
            TOTAL_Alt5_3 += sum(map(float, row["Alt5_coverages"].split(",") )) 

        if row["Alt3"]!='None':
            excluding_covs += sum(map(float, row["Alt3_coverages"].split(",") ))
            TOTAL_Alt5_3 += sum(map(float, row["Alt3_coverages"].split(",") ))
            
        half_counts = (ME_coverages - full_counts - TOTAL_Alt5_3)/2
        
        if half_counts <0:
            half_counts = 0
            
        ME_coverages_corrected = full_counts + half_counts
        
        outrow = "\t".join(map(str, [sample, ME, ME_coverages_corrected, excluding_covs], ))
        out.write(outrow + "\n")


with open(snakemake.output["count_spanning_ME_reads"], "wt") as out:
    
    header =  "\t".join(["ME", "Spanning_cov"])
    out.write(header + "\n")  

    for ME, cov in quant_ME_spanning_cov.items():
        
        outrow = "\t".join(map(str, [ME, cov], ))
        out.write(outrow + "\n")


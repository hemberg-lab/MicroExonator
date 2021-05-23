import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version

csv.field_size_limit(100000000)
csv.field_size_limit()

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
        
        





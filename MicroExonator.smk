#version 0.9.0

import yaml
from collections import defaultdict
import csv

configfile : "config.yaml"
DATA = set([])

def str2bool(v):
  if v==True:
    return True
  elif v==False:
    return False
  else:
    return v.lower() in ("yes", "true", "t", "1")

rule quant:
    input:
        "Report/out.high_quality.txt",
        expand("Report/quant/{sample}.out_filtered_ME.PSI.txt", sample=DATA)
        #"Report/out_filtered_ME.PSI.txt",        
        #"Report/stats/Microexons.not_consensus",
        #"Report/stats/Microexons.annotation.stats"        
        
        #"Report/out_filtered_ME.txt"
        #expand("Genome_aligments/{Software}/TOTAL.exons.{Software}", Software=["Hisat2", "STAR", "Olego"])
        # expand("Genome_aligments/{Software}/{sample}.sam.SJ_count", sample=DATA, Software=["Hisat2", "STAR"]),
        #expand("Whippet/Quant/{sample}.psi.gz", sample=DATA),
        #expand("Ground_Truth/{sample}.GT.SJ_count", sample=DATA)




if 'cluster_metadata' in config:

    cluster_files = defaultdict(list)
    cluster_files_metadata = defaultdict(list)
    single_cell_files = set([])

    with open(config["cluster_metadata"]) as Single_Cell:

        Single_Cell_clustering = csv.DictReader(Single_Cell, delimiter="\t")

        for row in Single_Cell_clustering:

            cluster_files[row[config["cluster_name"]].replace(" ", "_")].append(row[config["file_basename"]])
            single_cell_files.add(row[config["file_basename"]])



#### MicroExonator ####

if ("deletion_penalty" in config)==False:
    config["deletion_penalty"]="6"
    
if ("insertion_penalty" in config)==False:
    config["insertion_penalty"]="2"

config["indel_penalty"] = ",".join([str(config["deletion_penalty"]), str(config["insertion_penalty"])])

if ("ME_DB" in config)==False:
    config["ME_DB"]="touch/VastDb.bed12"

if ("paired_samples" in config)==False:
    config["paired_samples"]="F"

if ("min_reads_PSI" in config)==False:
    config["min_reads_PSI"]="5"


include : "rules/init.smk"
include : "rules/Get_data.smk"


rule bamfiles:
    input:
        expand("Whippet/BAM/{samples}.bam", samples=DATA), 
        expand("Whippet/BAM/{samples}.bam.bai", samples=DATA)


if str2bool(config.get("downstream_only", False)):
    pass
elif str2bool(config.get("skip_discovery_and_quant", False)):
    include : "rules/Round2_post_processing.smk"
elif str2bool(config.get("skip_discovery", False)):
    include : "rules/Round2.smk"
    include : "rules/Round2_post_processing.smk"
else:
    include : "rules/Round1.smk"
    include : "rules/Round1_post_processing.smk"
    include : "rules/Round2.smk"
    include : "rules/Round2_post_processing.smk"

rule discovery:
    input:
        expand("Round1/{sample}.sam.row_ME.filter1", sample=DATA )
#        "Round2/ME_canonical_SJ_tags.de_novo.fa"

##### Downstream Analysis ####

if "whippet_bin_folder" in config:
   include : "rules/Whippet_quant.smk"

if "whippet_delta" in config:
   with open(config["whippet_delta"], 'r') as stream:
      whippet_delta = yaml.safe_load(stream)
   include : "rules/Whippet_delta.smk"


#### Single Cell ###

if not "Single_Cell" in config:
   config["Single_Cell"]="F"

if str2bool(config["Single_Cell"]):
   include : "rules/Snakepool.py"
   include : "rules/pseudo_pool.smk"

#### Benchmark ####

#include : "rules/Benchmark.smk



#### Re-run incomplete round1 ####

import os

round1_incomplete = []

for file in DATA:
    if os.path.isfile('./Round1/' + file  + '.sam.row_ME.filter1')!=True:
        round1_incomplete.append(file)
  
rule rerun_incomplete_round1:
    input:
        expand("Round1/{sample}.sam.row_ME.filter1", sample=round1_incomplete )
        
        
        
round2_incomplete = []

for file in DATA:
    if os.path.isfile('./Round2/' + file  + '.sam.pre_processed.filter1.ME_SJ_coverage')!=True:
        round2_incomplete.append(file)
  
rule rerun_incomplete_round2:
    input:
        expand("Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage", sample=round2_incomplete )
    
    
include : "rules/sashimi.smk"
  

rule correct_quant:
    input:
        quant = "Report/quant/{sample}.out_filtered_ME.PSI.gz",
        ME_centric = "Round2/TOTAL.ME_centric.txt",
        spanning_ME_reads =  "Round2/ME_reads/{sample}.ME_spanning_reads.tsv"
    output:
        corrected_quant = "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz",
        count_spanning_ME_reads = "Round2/ME_reads/{sample}.tsv"
    script:
        "src/correct_quant.py"
		

rule get_all_corrected_quant:
    input:
        expand("Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz", sample=DATA)

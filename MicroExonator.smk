#version 0.9.0

import yaml
from collections import defaultdict
import csv
import random
import warnings
from snakemake.exceptions import WorkflowError
from src.filtering_config import (
    load_filter_groups,
    paired_read_samples,
    resolve_filter_method,
    selected_microexon_output,
)
random.seed(123)

configfile : "config.yaml"
DATA = set([])
to_validate = set([])
EBWT = [
    "1.ebwt",
    "2.ebwt",
    "3.ebwt",
    "4.ebwt",
    "rev.1.ebwt",
    "rev.2.ebwt",
]

warnings.simplefilter("default", DeprecationWarning)
try:
    FILTER_METHOD = resolve_filter_method(config)
except ValueError as error:
    raise WorkflowError(str(error))
FILTERED_ME_OUTPUT = selected_microexon_output(FILTER_METHOD)

if "validate_fastq_list" in config:

    with open(config["validate_fastq_list"]) as fastq_list:
        reader = csv.reader(fastq_list, delimiter="\t")
        for row in reader:
            to_validate.add(row[0])

def str2bool(v):
  if v==True:
    return True
  elif v==False:
    return False
  else:
    return v.lower() in ("yes", "true", "t", "1")

def hard_drive_behavior(wildcards):
    
    fastq = wildcards.sample

    if config.get("Optimize_hard_drive", False)=="T":

        if "validate_fastq_list" in config:

            if fastq in to_validate:
                return("FASTQ/round2/valid/" + fastq + ".valid.fastq.gz")
            else:
                return(  "FASTQ/round2/" + fastq + ".fastq.gz")

        else:
            return(  "FASTQ/round2/" + fastq + ".fastq.gz")
    else:

        if "validate_fastq_list" in config:

            if fastq in to_validate:
                return("FASTQ/valid/" + fastq + ".valid.fastq.gz")
            else:
                return(  "FASTQ/" + fastq + ".fastq.gz")
        else:

            return("FASTQ/" + fastq + ".fastq.gz")


cluster_files = defaultdict(list)
cluster_files_metadata = defaultdict(list)
single_cell_files = set([])



csv.field_size_limit(100000000)

pe_samples = set([])
paired_dict = dict()

if "paired_samples" in config:
    
    if config["paired_samples"]!="F":
    
        with open(config["paired_samples"]) as file:

            reader = csv.reader(file, delimiter="\t")
            for row in reader:

                pe_samples.add(row[0])
                pe_samples.add(row[1])
                paired_dict[row[0]] = row[1]
            
            
            
            
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

if "Single_Cell" not in config:
    config["Single_Cell"]="F"


include : "rules/init.smk"
include : "rules/Get_data.smk"


try:
    filter_groups = load_filter_groups(config, DATA, paired_dict=paired_dict)
except ValueError as error:
    raise WorkflowError(str(error))

sample_group_se = defaultdict(list, filter_groups["bulk_se"])
sample_group_pe = defaultdict(list, filter_groups["bulk_pe"])
sample_group_se_set = set(
    sample for samples in sample_group_se.values() for sample in samples
)
sample_group_pe_set = set(
    sample for samples in sample_group_pe.values() for sample in samples
)

primary_clusters = defaultdict(list, filter_groups["single_cell"])
cluster_files = defaultdict(
    list,
    {
        cluster.replace(" ", "_"): cells
        for cluster, cells in filter_groups["single_cell"].items()
    },
)
single_cell_files = set(
    sample for samples in primary_clusters.values() for sample in samples
)

def partition(list_in, n):
    randomized = list(list_in)
    random.shuffle(randomized)
    return [randomized[i::n] for i in range(n)]

pseudo_pool_dict = defaultdict(list)
pseudo_pool_dict_simple = dict()
cluster_pseudo_pools = defaultdict(list)
cluster_cells = defaultdict(list)

for cluster, cells in primary_clusters.items():
    cluster_id = cluster.replace(" ", "_")
    number_of_pools = min(3, len(cells))
    if cluster_id and number_of_pools:
        for pool_number, pool in enumerate(
            partition(cells, number_of_pools), start=1
        ):
            pseudo_pool_id = "{}-{}".format(cluster_id, pool_number)
            cluster_pseudo_pools[cluster_id].append(pseudo_pool_id)
            cluster_cells[cluster_id].extend(pool)
            for cell in pool:
                pseudo_pool_dict[pseudo_pool_id].append(cell)
                pseudo_pool_dict_simple[cell] = pseudo_pool_id

with open("pseudo_pool_membership.txt", "w") as pool_membership, open(
    "sample_groups.txt", "w"
) as sample_groups:
    for pseudo_pool, cells in pseudo_pool_dict.items():
        for cell in cells:
            pool_membership.write("{}\t{}\n".format(cell, pseudo_pool))
    for group, samples in sample_group_se.items():
        for sample in samples:
            sample_groups.write("{}\t{}\n".format(sample, group))
    for group, samples in sample_group_pe.items():
        for sample in samples:
            sample_groups.write("{}\t{}\n".format(sample, group))


rule quant:
    input:
        FILTERED_ME_OUTPUT,
        "Report/ME_junction_genomic_copies.txt",
        expand(
            "Report/quant/{sample}.out_filtered_ME.PSI.uncorrected.gz",
            sample=DATA,
        ),
        expand(
            "Report/quant/corrected/counts/{sample}.ME.adj_counts.gz",
            sample=DATA,
        )


rule bamfiles:
    input:
        expand("Whippet/BAM/{samples}.bam", samples=DATA), 
        expand("Whippet/BAM/{samples}.bam.bai", samples=DATA)


if str2bool(config.get("downstream_only", False)):
    pass
elif str2bool(config.get("skip_discovery_and_quant", False)):
    include : "rules/Round2_post_processing.smk"
    if FILTER_METHOD == "robustness":
        include : "rules/Round2_quant_filters.smk"
elif str2bool(config.get("skip_discovery", False)):
    include : "rules/Round2.smk"
    include : "rules/Round2_post_processing.smk"
    include : "rules/Round2_quant_filters.smk"
else:
    include : "rules/Round1.smk"
    include : "rules/Round1_post_processing.smk"
    include : "rules/Round2.smk"
    include : "rules/Round2_post_processing.smk"
    include : "rules/Round2_quant_filters.smk"
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

#include : "rules/Round2_quant_filters.smk"

#### Single Cell ###

if str2bool(config["Single_Cell"]) and "whippet_bin_folder" in config:
#   include : "rules/Snakepool.py"
    include : "rules/pseudo_pool.smk"
    ruleorder: quant_pool_pb > whippet_quant
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
    

rule get_whippet_psi:
    input:
        expand("Whippet/Quant/{sample}.psi.gz", sample=DATA)

    
include : "rules/sashimi.smk"
  
#ruleorder: quant_pool_pb > whippet_quant

import sys
import csv
import gzip
from collections import defaultdict
from snakemake.utils import min_version
import glob
import itertools
import csv
import gzip
from collections import defaultdict
import random
random.seed(123)

#This script gets the sample groups (output 1) and compute quantitative filters to get a list of reliable detected microexons.

# csv.field_size_limit(100000000)


  
# def partition (list_in, n):  # Function to do random pooling
#     random.shuffle(list_in)
#     return [list_in[i::n] for i in range(n)]
  

# primary_clusters = defaultdict(list)

# pe_samples = set([])
# paired_dict = dict()

# if "paired_samples" in config:
    
#     if config["paired_samples"]!="F":
    
#         with open(config["paired_samples"]) as file:

#             reader = csv.reader(file, delimiter="\t")
#             for row in reader:

#                 pe_samples.add(row[0])
#                 paired_dict[row[0]] = row[1]

if 'cluster_metadata' in config:
    with open(  config["cluster_metadata"]) as file:
  
        reader = csv.DictReader(file, delimiter="\t")
    
        for row in reader:
            DATA.add(row[config["file_basename"]])

#print(DATA)
    
#         primary_clusters[row[config["cluster_name"]]].append(row[config["file_basename"]])



# def partition (list_in, n):  # Function to do random pooling
#     random.shuffle(list_in)
#     return [list_in[i::n] for i in range(n)]
    

# sample_group_se = defaultdict(list)
# sample_group_pe = defaultdict(list)
# sample_group_se_set = set()
# sample_group_pe_set = set()

# pseudo_pool_dict =  defaultdict(list)
# pseudo_pool_dict_simple = dict()
# cluster_pseudo_pools = defaultdict(list)
# cluster_cells = defaultdict(list)

# for cluster in primary_clusters:
    
#     c = 0
#     cells = primary_clusters[cluster]
#     pseudo_pools = partition( cells , 3 )
#     pseudo_pool_ID = ""

#     if cluster !="":    

#         for pool in pseudo_pools:
        
#             c+=1
#             pseudo_pool_ID =cluster.replace(" ", "_") + "-" + str(c)
#             cluster_pseudo_pools[cluster.replace(" ", "_")].append(pseudo_pool_ID)	
 
#             for cell in pool:
#                 pseudo_pool_dict[pseudo_pool_ID].append(cell)
#                 cluster_cells[cluster.replace(" ", "_")].append(cell)
#                 pseudo_pool_dict_simple[cell] = pseudo_pool_ID	 
            
# if "bulk_samples" in config:
           
#     with open(config["bulk_samples"]) as file:

#         reader = csv.DictReader(file, delimiter="\t")
#         for row in reader:
#             if row["sample"] in pe_samples:
#                 if row["sample"] in paired_dict:
#                     sample_group_pe[row["condition"]].append(row["sample"])
#                     #sample_group_pe[row["sample"]] = row["condition"]
#                     sample_group_pe_set.add(row["sample"])
#             else:
#                 sample_group_se[row["condition"]].append(row["sample"]),
#                 #sample_group_se[row["sample"]] = row["condition"]
#                 sample_group_se_set.add(row["sample"])
        

# with open("pseudo_pool_membership.txt", "w") as out_pseudo_pool_membership, open("sample_groups.txt", "w") as out_sample_groups:
    
#     for sp, cells  in pseudo_pool_dict.items():  
#         for cell in cells:
#             out = "\t".join([cell, sp])  
#             out_pseudo_pool_membership.write(out + "\n")
        
    
#     for  group, samples in sample_group_se.items():
#         for sample in samples:
#             out = "\t".join([sample, group])
#             out_sample_groups.write(out + "\n")
        
#     for  group, sample in sample_group_pe.items():
#         for sample in samples:
#             out = "\t".join([sample, group])
#             out_sample_groups.write(out + "\n")
        
if str2bool(config.get("optimise_disk", False)):
    rule correct_quant:
        input:
            quant = "Report/quant/{sample}.out_filtered_ME.PSI.uncorrected.gz",
            ME_centric = "Round2/TOTAL.ME_centric.txt",
            spanning_ME_reads =  "Round2/ME_reads/{sample}.ME_spanning_reads.tsv"
        output:
            corrected_quant = temp("Report/quant/corrected/counts/{sample}.ME.adj_counts.gz"),
            count_spanning_ME_reads = "Round2/ME_reads/{sample}.counts.tsv"
        priority: 300
        script:
            "../src/correct_quant.py"
else:
    rule correct_quant:
        input:
            quant = "Report/quant/{sample}.out_filtered_ME.PSI.uncorrected.gz",
            ME_centric = "Round2/TOTAL.ME_centric.txt",
            spanning_ME_reads =  "Round2/ME_reads/{sample}.ME_spanning_reads.tsv"
        output:
            corrected_quant = protected("Report/quant/corrected/counts/{sample}.ME.adj_counts.gz"),
            count_spanning_ME_reads = "Round2/ME_reads/{sample}.counts.tsv"
        priority: 300
        script:
            "../src/correct_quant.py"
		

rule get_all_corrected_quant:
    input:
        expand("Report/quant/corrected/counts/{sample}.ME.adj_counts.gz", sample=DATA)        

#print(expand("Report/quant/corrected/counts/{sample}.ME.adj_counts.gz", sample=DATA))        
        
def get_cell_sp(cluster):
    return(expand( "Report/quant/corrected/counts/{sample}.ME.adj_counts.gz", sample = pseudo_pool_dict[cluster]))
    
rule get_PSI_sparse_quants_sp:
    input:
        cells = lambda w : get_cell_sp(w.cluster)
    output:
        corrected_sparse = protected("Report/quant/sparse/single_cell/{cluster}.corrected.PSI.gz")
    priority: 10
    script:
        "../src/get_sparse_quants_sp.py" 
                
rule get_PSI_sparse_quants_se:
    input:
        corrected_quant = "Report/quant/corrected/counts/{sample}.ME.adj_counts.gz"
    output:
        corrected_sparse = protected("Report/quant/corrected/PSI_sparse/bulk/se/{sample}.corrected.PSI.gz")
    priority: 10
    script:
        "../src/get_sparse_quants_se.py"

          
def get_pair(rd1):
    return(paired_dict[rd1])
           
rule get_PSI_sparse_quants_pe:
    input:
        corrected_quant_rd1 = "Report/quant/corrected/counts/{sample}.ME.adj_counts.gz",
        corrected_quant_rd2 = lambda w : expand( "Report/quant/corrected/counts/{rd2}.ME.adj_counts.gz", rd2=paired_dict[w.sample])
    output:
        corrected_sparse = protected("Report/quant/corrected/PSI_sparse/bulk/pe/{sample}.corrected.PSI.gz")
    priority: 10
    script:
        "../src/get_sparse_quants_pe.py"

rule detection_filter_se:
    input:
        PSI_files = lambda w : expand("Report/quant/corrected/PSI_sparse/bulk/se/{sample}.corrected.PSI.gz", sample = sample_group_se[w.sample_group] ),
        ME_reads = lambda w : expand("Round2/ME_reads/{sample_group}.counts.tsv", sample_group = sample_group_se[w.sample_group] )
    output:
        detected = "Report/filter/se/{sample_group}.detected.txt"
    script:
        "../src/detected_me.py"
        
rule detection_filter_pe:
    input:
        PSI_files = lambda w : expand("Report/quant/corrected/PSI_sparse/bulk/pe/{sample}.corrected.PSI.gz", sample = sample_group_pe[w.sample_group] ),
        ME_reads = lambda w : expand("Round2/ME_reads/{sample_group}.counts.tsv", sample_group = sample_group_pe[w.sample_group] )
    output:
        detected = "Report/filter/pe/{sample_group}.detected.txt"
    script:
        "../src/detected_me.py"
	
rule detection_filter_sp:
    input:
        PSI_files = lambda w : expand("Report/quant/corrected/PSI_sparse/single_cell/{pseudo_pool}.corrected.PSI.gz", pseudo_pool = cluster_pseudo_pools[w.cluster] ),
	ME_reads = lambda w : expand("Round2/ME_reads/{pseudo_pool}.counts.tsv", pseudo_pool = cluster_cells[w.cluster] )
    output:
        detected = "Report/filter/sc/{cluster}.detected.txt"
    script:
        "../src/detected_me.py"	

#print(cluster_pseudo_pools)
#print(cluster_pseudo_pools)

rule detection_filter_total:
    input:
        expand("Report/filter/se/{sample_group}.detected.txt",  sample_group=sample_group_se.keys()),
        expand("Report/filter/pe/{sample_group}.detected.txt",  sample_group=sample_group_pe.keys()),                
        expand("Report/filter/sc/{cluster}.detected.txt",  cluster=cluster_pseudo_pools.keys())
        
        
rule detection_filter:
    input : 
        bulk_se = expand("Report/filter/se/{sample_group}.detected.txt",  sample_group=sample_group_se.keys()),
        bulk_pe = expand("Report/filter/pe/{sample_group}.detected.txt",  sample_group=sample_group_pe.keys()),
        single_cell = expand("Report/filter/sc/{cluster}.detected.txt",  cluster=cluster_pseudo_pools.keys()),
        ME_centric = "Round2/TOTAL.ME_centric.txt"   
    params:
        min_detected_samples = int(config.get("min_detected_samples", 1))
    output:
       detected_list = "Report/out.robustly_detected.txt"
    priority: 10
    script:
        "../src/final_filters2.py"  


## GENE EXPRESSION


rule get_transcriptome:
    input:
        genome = config["Genome_fasta"],
        gtf = config["Gene_anontation_GTF"]
    output:
        "Genome/transcriptome.fa"
    conda: "../envs/core.yaml"
    shell:
        "gffread -w {output} -g {input}"

           
rule salmon_index:
    input:
        "Genome/transcriptome.fa"
    output:
        directory("salmon/transcriptome_index")
    log:
        "logs/salmon/transcriptome_index.log"
    threads: 2
    params:
        extra=""
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i {output} --threads {threads} "
           
 
rule salmon_quant_se:
    input:
        r = "FASTQ/{sample}.fastq.gz",
        #r = "FASTQ/{sample}.fastq.gz",
        index = "salmon/transcriptome_index"
    output:
        quant = 'salmon/SE/{sample}/quant.sf'
    log:
        'logs/salmon/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        out_dir = 'salmon/SE/{sample}',
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra=""
    threads: 2
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -i {input.index} -l {params.libtype} -r {input.r}  -p {threads}  -o {params.out_dir}"        


rule salmon_quant_reads_pe:
    input:
        rd1 = "FASTQ/{sample}.fastq.gz",
        rd2 = lambda w : expand( "FASTQ/{rd2}.fastq.gz", rd2=paired_dict[w.sample]), 
        index = "salmon/transcriptome_index"
    output:
        quant = 'salmon/PE/{sample}/quant.sf'
    log:
        'logs/salmon/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        out_dir = 'salmon/PE/{sample}',
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra=""
    threads: 2
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -i {input.index} -l {params.libtype}  -1 {input.rd1} -2 {input.rd2} -p {threads}  -o {params.out_dir}"


rule salmon_all_quant:
    input:
        expand( 'salmon/SE/{sample}/quant.sf', sample=sample_group_se_set), 
        expand( 'salmon/PE/{sample}/quant.sf', sample=sample_group_pe_set) 

           

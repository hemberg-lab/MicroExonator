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

csv.field_size_limit(100000000)


  
def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
  

primary_clusters = defaultdict(list)

pe_samples = set([])
paired_dict = dict()

if "paired_samples" in config:
    
    with open(config["paired_samples"]) as file:
        
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            
            pe_samples.add(row[0])
            paired_dict[row[0]] = row[1]


with open(config["cluster_metadata"]) as file:
  
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
        primary_clusters[row[config["cluster_name"]]].append(row[config["file_basename"]])



def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]
    

sample_group_se = defaultdict(set)
sample_group_pe = defaultdict(set)
sample_group_se_set = set()
sample_group_pe_set = set()

pseudo_pool_dict =  dict()


for cluster in primary_clusters:
    
    c = 0
    cells = primary_clusters[cluster]
    pseudo_pools = partition( cells , 3 )
    pseudo_pool_ID = ""
    
    for pool in pseudo_pools:
        
        c+=1
        pseudo_pool_ID =cluster.replace(" ", "_") + "-" + str(c)
        
        for cell in pool:
            pseudo_pool_dict[pseudo_pool_ID] = cell
           
with open(config["bulk_samples"]) as file:
    
    reader = csv.DictReader(file, delimiter="\t")
    for row in reader:
        if row["sample"] in pe_samples:
            if row["sample"] in paired_dict:
                sample_group_pe[row["condition"]].add(row["sample"])
                sample_group_pe_set.add(row["sample"])
        else:
            sample_group_se[row["condition"]].add(row["sample"])
            sample_group_se_set.add(row["sample"])
        

with open("pseudo_pool_membership.txt", "w") as out_pseudo_pool_membership, open("sample_groups.txt", "w") as out_sample_groups:
    
    for sp, cell  in pseudo_pool_dict.items():
        
        out = "\t".join([cell, sp])  
        out_pseudo_pool_membership.write(out + "\n")
    
    for  group, samples in sample_group_se.items():
        for sample in samples:
            out = "\t".join([sample, group])
            out_sample_groups.write(out + "\n")
        
    for  group, sample in sample_group_pe.items():
        for sample in samples:
            out = "\t".join([sample, group])
            out_sample_groups.write(out + "\n")
        
        
def get_cell_sp(cluster):
    return(expand( "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz", sample = pseudo_pool_dict[cluster]))
    
rule get_sparse_quants_sp:
    input:
        cells = lambda w : get_cell_sp(w.cluster)
    output:
        corrected_sparse = "Report/quant/sparse/single_cell/{cluster}.corrected.PSI.gz"
    priority: 10
    script:
        "../src/get_sparse_quants_sp.py" 
                
rule get_sparse_quants_se:
    input:
        corrected_quant = "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz"
    output:
        corrected_sparse = "Report/quant/sparse/bulk/se/{sample}.corrected.PSI.gz"
    priority: 10
    script:
        "../src/get_sparse_quants_se.py"

          
def get_pair(rd1):
    return(paired_dict[rd1])
           
rule get_sparse_quants_pe:
    input:
        corrected_quant_rd1 = "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz",
        corrected_quant_rd2 = lambda w : expand( "Report/quant/corrected/{rd2}.out_filtered_ME.PSI.gz", rd2=paired_dict[w.sample])
    output:
        corrected_sparse = "Report/quant/sparse/bulk/pe/{sample}.corrected.PSI.gz"
    priority: 10
    script:
        "../src/get_sparse_quants_pe.py"
           
         
rule detection_filter:
    input : 
        bulk_se = expand("Report/quant/sparse/bulk/se/{sample}.corrected.PSI.gz", sample = sample_group_se_set),
        bulk_pe = expand("Report/quant/sparse/bulk/pe/{sample}.corrected.PSI.gz", sample = sample_group_pe_set),
        single_cell = expand("Report/quant/sparse/single_cell/{cluster}.corrected.PSI.gz", cluster = pseudo_pool_dict.keys()),   
        bulk_ME_reads_se = expand( "Round2/ME_reads/{sample}.counts.tsv",  sample = sample_group_se_set),
        bulk_ME_reads_pe = expand( "Round2/ME_reads/{sample}.counts.tsv",  sample = sample_group_pe_set),
        single_cell_reads = expand( "Round2/ME_reads/{cluster}.counts.tsv",  cluster = set(pseudo_pool_dict.values()))
    params:
        se_dict = sample_group_se,
        pe_dict = sample_group_pe,
        pseudo_dict = pseudo_pool_dict
    output:
       detected_list = "Report/filter/robustly_detected_ME.txt"
    params:
       min_PSI = config["min_PSI"]
    priority: 10
    script:
        "../src/final_filters.py"  


## GENE EXPRESSION
           
rule salmon_index:
    input:
        #"gffcompare/extended_ref_annotation.fa"
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
        quant = 'salmon/{sample}/quant.sf',
        lib = 'salmon/{sample}/lib_format_counts.json'
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
        r1 = "FASTQ/{sample}.fastq.gz",
        r2 = lambda w : "FASTQ/" + get_pair(w.sample) + ".fastq.gz", 
        index = "salmon/transcriptome_index"
    output:
        quant = temp('salmon/{sample}/quant.sf'),
        lib = 'salmon/{sample}/lib_format_counts.json'
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
        expand( 'salmon/SE/{sample}/quant.sf', sample=sample_group_se.keys()), 
        expand( 'salmon/PE/{sample}/quant.sf', sample=sample_group_pe.keys()) 

           

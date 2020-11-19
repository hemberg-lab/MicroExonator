import glob, os
import random
import csv
import gzip
from collections import defaultdict


def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

n_sb = 5

cluster_files_pb = dict()

for cluster, files in cluster_files.items():
    sb = 1
    for pool in partition(files, n_sb):
        cluster_files_pb[(cluster, sb)] = pool
        sb += 1

def get_files_by_cluster_pb(cluster, pool_ID):
    ext = ".fastq.gz"
    path="FASTQ/"
    return([path + x + ext for x in cluster_files_pb[(cluster, int(pool_ID))]])

rule quant_pool_pb:
    input:
        fastq = lambda w: get_files_by_cluster_pb(w.cluster, w.pool_ID),
	index = "Whippet/Index/whippet.jls"
    output:
        "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.gene.tpm.gz",
	"Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.isoform.tpm.gz",
        "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.jnc.gz",
        "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.map.gz",
        "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz"
    params:
        bin = config["whippet_bin_folder"],
        output = "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}"
    priority: 10
    shell:
        "julia {params.bin}/whippet-quant.jl <( cat {input.fastq} ) --force-gz -x {input.index}  -o {params.output}"
        
#print(expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))))
        
rule get_pseudo_pools:
    input:
        expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  )))
	
	
rule collapse_pseudo_pools:
  input: 
      gene = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.gene.tpm.tsv",
      isoform = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.isoform.tpm.tsv"


rule merge_quant_by_cluster_gene_sp:
    input:
        files = expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.gene.tpm.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))),
        jnc =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.jnc.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))),
        mapf =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.map.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))),
        psi =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  )))
    params:
        feature = "Gene"
    output:
        merged = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.gene.tpm.tsv"
    script:
        "../src/merge_quant.py"


rule merge_quant_by_cluster_isoform_sp:
    input:
        files = expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.gene.tpm.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))),
        jnc =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.jnc.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))),
        mapf =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.map.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))),
        psi =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  )))
    params:
        feature = "Isoform"
    output:
        merged = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.isoform.tpm.tsv"
    script:
        "../src/merge_quant.py"
        

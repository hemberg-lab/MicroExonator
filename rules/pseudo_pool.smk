import glob, os
import random
import csv
import gzip
from collections import defaultdict


def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

#n_sb = 5
if "cells_pseudobulks" in config:
    n_cells = int(config["cells_pseudobulks"])
else:
    n_cells = 15 

cluster_files_pb = dict()
sb_IDs = set()

for cluster, files in cluster_files.items():
    sb = 1
    
    if "n_pseudobulks" in config:
        n_sb = int(config["n_pseudobulks"])
    
    else:
        n_sb = round(len(files)/n_cells)
        if n_sb<3:
            n_sb=3
        
    for pool in partition(files, n_sb):
        cluster_files_pb[(cluster, sb)] = pool
        sb_IDs.add(cluster + "_" + str(sb))
        sb += 1

rule get_pseudo_bulk_membership:
    output:
        table = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulk_membership.tsv"
    run:
        with open(output.table, "w") as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(["pseudo_bulk_ID", "samples"])
            
            for key in cluster_files_pb:
                cluster, sb = key
                sb_ID = cluster + "_" + str(sb)
                samples = ",".join(cluster_files_pb[key])
                writer.writerow([sb_ID, samples])
                

def get_files_by_cluster_pb(cluster, pool_ID):
    ext = ".fastq.gz"
    path="FASTQ/"
    return([path + x + ext for x in cluster_files_pb[(cluster, int(pool_ID))]])

if str2bool(config.get("keep_uncollapsed_pseudobulk", False)):

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
            julia = config["julia"],
            output = "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}"
        priority: 10
        shell:
            "{params.julia} {params.bin}/whippet-quant.jl <( cat {input.fastq} ) --force-gz -x {input.index}  -o {params.output}"
            
else:	

    rule quant_pool_pb:
        input:
            fastq = lambda w: get_files_by_cluster_pb(w.cluster, w.pool_ID),
            index = "Whippet/Index/whippet.jls"
        output:
            temp("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.gene.tpm.gz"),
            temp("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.isoform.tpm.gz"),
            temp("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.jnc.gz"),
            temp("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.map.gz"),
            temp("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz")
        params:
            bin = config["whippet_bin_folder"],
            julia = config["julia"],
            output = "Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}"
        priority: 10
        shell:
            "{params.julia} {params.bin}/whippet-quant.jl <( cat {input.fastq} ) --force-gz -x {input.index}  -o {params.output}"
        
#print(expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{cluster}_{pool_ID}.psi.gz", cluster=cluster_files.keys(), pool_ID=list(range(1, n_sb+1  ))))
        
rule get_pseudo_pools:
    input:
        expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.psi.gz", pseudo_pool=sb_IDs)

	
	
rule collapse_pseudo_pools:
  input: 
      gene = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.gene.tpm.tsv.gz",
      isoform = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.isoform.tpm.tsv.gz",
      psi =  "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.psi.tsv.gz",
      table =  "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulk_membership.tsv"


rule merge_quant_gene_sp:
    input:
        files = expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.gene.tpm.gz",  pseudo_pool=sb_IDs),
        jnc =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.jnc.gz", pseudo_pool=sb_IDs),
        mapf =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.map.gz", pseudo_pool=sb_IDs),
        psi =  expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.psi.gz", pseudo_pool=sb_IDs)
    params:
        feature = "Gene"
    output:
        merged = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.gene.tpm.tsv.gz"
    script:
        "../src/merge_quant.py"


rule merge_quant_isoform_sp:
    input:
        files = expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.isoform.tpm.gz", pseudo_pool=sb_IDs)
    params:
        feature = "Isoform"
    output:
        merged = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.isoform.tpm.tsv.gz"
    script:
        "../src/merge_quant.py"
	
rule merge_quant_PSI_sp:
    input:
        files = expand("Whippet/Quant/Single_Cell/Pseudo_bulks/{pseudo_pool}.psi.gz", pseudo_pool=sb_IDs)
    params:
        feature = "PSI"
    output:
        merged = "Whippet/Quant/Single_Cell/Pseudo_bulks/pseudo_bulks.psi.tsv.gz"
    script:
        "../src/merge_quant.py"


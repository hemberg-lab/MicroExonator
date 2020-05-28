###### params #####

#wildcard_constraints:
#    sample="^[A-Za-z0-9_-]*$"

random.seed(123) 

#repeats = 10  #now the repeats are specified on the tsv

###############

import glob, os
import random
import csv
import gzip
from collections import defaultdict


def partition (list_in, n):  # Function to do random pooling
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]


cluster_compare = dict()
cluster_compare_np = dict()
compare_repeats = dict()
## Structure of dictionaries

#cluster_compare = { "Neuronal-vs-Non_Neuronal" :( ["GABA-ergic Neuron",  "Glutamatergic Neuron"], ["Endothelial Cell", "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte Precursor Cell" ] ),
#                    "GABA-ergic_Neuron_vs_Glutamatergic_Neuron" :( ["GABA-ergic Neuron"], ['Glutamatergic Neuron'] )     }


#cluster_compare_np = { "Neuronal-vs-Non_Neuronal" :( 100, 10 ),
#                    "GABA-ergic_Neuron_vs_Glutamatergic_Neuron" :( 50, 50 )     }




with open(config["run_metadata"]) as run:   #Populating the dictionaries

    run_metadata = csv.DictReader(run, delimiter="\t")

    for row in run_metadata:

        A_cluster_names = []
        B_cluster_names = []

        for c in row["A.cluster_names"].split(","):

            A_cluster_names.append(c)

        for c in row["B.cluster_names"].split(","):

            B_cluster_names.append(c)

        cluster_compare[row["Compare_ID"]] = (A_cluster_names, B_cluster_names)
        cluster_compare_np[row["Compare_ID"]] = (int(row["A.number_of_pools"]), int(row["B.number_of_pools"]))
        compare_repeats[row["Compare_ID"]] = int(row["Repeat"])


cluster_files = defaultdict(list)


#cluster_files = {"GABA" : ["fileA", ... ] }

with open(config["cluster_metadata"]) as Single_Cell:

    Single_Cell_clustering = csv.DictReader(Single_Cell, delimiter="\t")

    for row in Single_Cell_clustering:

        cluster_files[row[config["cluster_name"]]].append(row[config["file_basename"]])



target_pool_psi = []
target_pool_delta = []


#cluster_compare = { "Neuronal-vs-Non_Neuronal" :( ["GABA-ergic Neuron",  "Glutamatergic Neuron"], ["Endothelial Cell", "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte Precursor Cell" ] ),
#                    "GABA-ergic_Neuron_vs_Glutamatergic_Neuron" :( ["GABA-ergic Neuron"], ['Glutamatergic Neuron'] )     }


#cluster_compare.items()

#[ [ ["Neuronal-vs-Non_Neuronal"],  ( ["GABA-ergic Neuron",  "Glutamatergic Neuron"], ["Endothelial Cell", "Astrocyte", "Microglia", "Oligodendrocyte", "Oligodendrocyte Precursor Cell" ] ), [ ]  ] ]


#for compare_name, c in cluster_compare.items():  #Getting the target files - key = compare_name ; value = c

    # g1, g2 = c
    #
    # c1_names = []
    # for c1 in g1:
    #
    #     c1_names += cluster_files[c1]
    #
    # c2_names = []
    # for c2 in g2:
    #     c2_names += cluster_files[c2]

#print(cluster_compare)    
    
compare_names = []


for compare_name in cluster_compare.keys():  #Getting the target files - key = compare_name ; value = c

    compare_names.append(compare_name)

    for r in range(compare_repeats[compare_name]):

        delta_name = "Whippet/Delta/Single_Cell/" + compare_name +  "_rep_" +  str(r+1)

        #print(delta_name)
        
        if str2bool(config.get("Only_snakepool", False)):
            
            target_pool_delta.append( delta_name + ".diff.gz")
            
        else:    
       
            target_pool_delta.append( delta_name + ".diff.microexons")


if str2bool(config.get("Only_snakepool", False)):
    
 rule snakepool:   # This rule execute all the nesesary rules to produce the target files
   input:
    target_pool_delta ,#target files
    expand("Whippet/Delta/Single_Cell/Unpooled/{compare_name}.diff.gz", compare_name=compare_names)   
    
else:
      
    rule snakepool:   # This rule execute all the nesesary rules to produce the target files
       input:
        target_pool_delta ,#target files
        expand("Whippet/Delta/Single_Cell/Unpooled/{compare_name}.diff.microexons", compare_name=compare_names)


##### Single cell ####

import random
import csv
from collections import defaultdict

def partition (list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]



cluster_files = defaultdict(list)


#with open("/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator_Final/Whippet/Single_Cell_clustering.txt") as Single_Cell:
with open(config["cluster_metadata"]) as SC:

    Single_Cell_clustering = csv.DictReader(SC, delimiter="\t")

    for row in Single_Cell_clustering:


        cluster_files[row[config["cluster_name"]]].append(row[config["file_basename"]])


delta_unpooled_dict = dict()

for compare_name, c in cluster_compare.items():


    g1, g2 = c

    c1_names = []
    for c1 in g1:

        c1_names += cluster_files[c1]

    c2_names = []
    for c2 in g2:
        c2_names += cluster_files[c2]
        
    delta_unpooled_dict[(compare_name, "A")] = c1_names
    delta_unpooled_dict[(compare_name, "B")] = c2_names


    np_A, np_B = cluster_compare_np[compare_name]


    ## Unpooled analysis

rule delta_unpool:
    input:
        lambda w : expand("Whippet/Quant/{sample}.psi.gz", sample=delta_unpooled_dict[(w.compare_name, "A")]) + expand("Whippet/Quant/{sample}.psi.gz", sample=delta_unpooled_dict[(w.compare_name, "B")]) 
        #expand("Whippet/Quant/{sample}.psi.gz", sample=c1_names) + expand("Whippet/Quant/{sample}.psi.gz", sample=c2_names)
    output:
        "Whippet/Delta/Single_Cell/Unpooled/{compare_name}.run.sh"
    params:
        bin = config["whippet_bin_folder"],
        a = ",".join(expand("Whippet/Quant/{sample}.psi.gz", sample=delta_unpooled_dict[(compare_name, "A")])),
        b = ",".join(expand("Whippet/Quant/{sample}.psi.gz", sample=delta_unpooled_dict[(compare_name, "B")])),
        o = "Whippet/Delta/Single_Cell/Unpooled/{compare_name}"
    shell:
        "echo julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o} > {output}"


rule run_delta_unpool:  #to avoid overload shell comandline
    input:
        "Whippet/Delta/Single_Cell/Unpooled/{compare_name}.run.sh"
    output:
        "Whippet/Delta/Single_Cell/Unpooled/{compare_name}.diff.gz"
    shell:
        "bash {input}"



#pseudo pooling        

pool_dict_quant = dict()
pool_dict_delta = dict()

for compare_name, c in cluster_compare.items():
    
    for r in range(compare_repeats[compare_name]):


        c1_pools = partition(c1_names, np_A)
        c2_pools = partition(c2_names, np_B)

        

        target_pool_psi_A = []
        target_pool_psi_B = []

        delta_name = "Whippet/Delta/Single_Cell/" + compare_name +  "_rep_" +  str(r+1)

        #for pc1, pc2 in zip(c1_pools, c2_pools):

        p = 0

        for pc1 in c1_pools:

            p += 1

            FASTQ_c1 = [ "FASTQ/" + x + ".fastq.gz" for x in  pc1 ]

            PSI_c1 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  pc1 ]

            pool_ID = str(r + 1) + "-"  + str(p)
            pool_dict_quant[(compare_name, pool_ID, "A")] = FASTQ_c1
            #pool_dict_delta[(delta_name, "A")] = PSI_c1

            target_pool_psi_A.append("Whippet/Quant/Single_Cell/" + compare_name + "_A_" + pool_ID + ".psi.gz")

        p = 0
            
        for pc2 in c2_pools:

            p += 1

            FASTQ_c2 = [ "FASTQ/" + x + ".fastq.gz" for x in  pc2 ]
            PSI_c2 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  pc2 ]

            pool_ID = str(r + 1) + "-"  + str(p)
            pool_dict_quant[(compare_name, pool_ID, "B")] = FASTQ_c2
            #pool_dict_delta[(delta_name, "B")] = PSI_c2

            target_pool_psi_B.append("Whippet/Quant/Single_Cell/" + compare_name + "_B_" + pool_ID + ".psi.gz")

        pool_dict_delta[(delta_name, "A")] = target_pool_psi_A
        pool_dict_delta[(delta_name, "B")] = target_pool_psi_B

rule quant_pool:
    input:
        fastq = lambda w: pool_dict_quant[(w.compare_name, w.pool_ID, w.cond)],
        index = "Whippet/Index/whippet.jls"
    output:
        "Whippet/Quant/Single_Cell/{compare_name}_{cond}_{pool_ID}.gene.tpm.gz",
        "Whippet/Quant/Single_Cell/{compare_name}_{cond}_{pool_ID}.isoform.tpm.gz",
        "Whippet/Quant/Single_Cell/{compare_name}_{cond}_{pool_ID}.jnc.gz",
        "Whippet/Quant/Single_Cell/{compare_name}_{cond}_{pool_ID}.map.gz",
        "Whippet/Quant/Single_Cell/{compare_name}_{cond}_{pool_ID}.psi.gz"
    params:
        bin = config["whippet_bin_folder"],
        output = "Whippet/Quant/Single_Cell/{compare_name}_{cond}_{pool_ID}"
    priority: 10
    shell:
        "julia {params.bin}/whippet-quant.jl <( cat {input.fastq} ) --force-gz -x {input.index}  -o {params.output}"
        

        
rule delta_pool:
    input:
        A = lambda w: pool_dict_delta[(w.delta_name, "A")],
        B = lambda w: pool_dict_delta[(w.delta_name, "B")]
    output:
        "{delta_name}.diff.gz"
    params:
        bin = config["whippet_bin_folder"],
        a = lambda w, input: ",".join( input.A ),
        b = lambda w, input: ",".join( input.B ),
        o = "{delta_name}"
    shell:
        "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"
 

        
     #### these rules gereate a single indexed bam per condition which can be used for visualization


        
rule merge_bam:
    input:
        lambda w: expand('Whippet/BAM/{sample}.bam', sample=cluster_files[w.compare_name])
    output:
        temp("Whippet/BAM/Merge/{compare_name}.bam.merge")
    shell:
        "samtools merge  {output} {input}"

rule sort_index_bam:
    input:
        "Whippet/BAM/Merge/{compare_name}.bam.merge"
    output:
        "Whippet/BAM/Merge/{compare_name}.sort.bam"
    shell:
        'samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output}'
        
rule cluster_bams:
    input:
        expand("Whippet/BAM/Merge/{compare_name}.sort.bam", compare_name=cluster_files.keys())  
               
# # # # #               
          
if str2bool(config.get("cluster_sashimi", False)):
    
    gene_nodes = dict()
    node_strand = dict()
    
    compare_sig_nodes = dict()
    sashimis = set([])
    
    psi_file = random.choice(glob.glob('Whippet/Quant/Single_Cell/*.psi.gz'))
    
    
    with gzip.open(psi_file, mode="rt") as f:
        
        reader = csv.DictReader(f, delimiter="\t")
        
        for row in reader:
            gene_nodes[(row["Gene"], int(row["Node"]) ) ] = row["Coord"]
            node_strand[row["Coord"]] = row["Strand"]
            
    for sig_node_file in glob.glob('Whippet/Delta/Single_Cell/Sig_nodes/*'):
        
        compare_name = sig_node_file.split("/")[-1].split(".")[0]
        
        with open(sig_node_file) as file:
            
            reader = csv.DictReader(file, delimiter="\t")
            
            for row in reader:
                compare_sig_nodes[compare_name] = "_".join([row["Gene"], row["Node"], row["Strand"]])
                sashimis.add("Whippet/ggsashimi/" + compare_name + "/" + "_".join([row["Gene"], row["Node"], row["Strand"]]))
                

    rule get_bam_tsv:
        input:
            config["run_metadata"]
        output:
            expand("Whippet/ggsashimi/{compare_name}/{compare_name}.tvs" , compare_name=compare_names)
        script:
            "../src/write_bam_tsv.py"
    
    
    rule get_sig_nodes:
        params:
            path = "Whippet/Delta/Single_Cell/Sig_nodes/"
        output:
            temp(expand("{sashimi}.txt", sashimi=sashimis))
        script:
            "../src/write_sig_node_files.py"
    
    def coord_to_region(gene, node, strand):
        
        node = int(node)
        node_up = 1
        node_down = node
        
        if node > 1:
            node_up = node-1
        if (gene, node+1) in gene_nodes:
            node_down = node+1
        
        if node_up > 1:
            node_up = node_up-1
        if (gene, node_down+1) in gene_nodes:
            node_down = node_down+1   
        
        if (gene, node_up) in gene_nodes:
            node_up_coord = gene_nodes[(gene, node_up)]
        else:
            while ((gene, node_up) in gene_nodes)==False:
                node_up+= -1
                if (gene, node_up) in gene_nodes:
                    node_up_coord = gene_nodes[(gene, node_up)]
                    break
                if node_up < 1:
                    node_up = node
                    node_up_coord = gene_nodes[(gene, node_up)]
                    break

        if (gene, node_down) in gene_nodes:
            node_down_coord = gene_nodes[(gene, node_down)]
        else:
            while ((gene, node_down) in gene_nodes)==False:
                node_down+= 1
                if (gene, node_down) in gene_nodes:
                    node_down_coord = gene_nodes[(gene, node_down)]
                    break
                if node_up > 1000:
                    node_up = node
                    node_down_coord = gene_nodes[(gene, node_down)]
                    break
        
        chrom = node_up_coord.split(":")[0]
        start = node_up_coord.split(":")[1].split("-")[0]
        end = node_down_coord.split(":")[1].split("-")[1]
        
        if strand=="-":
            end = node_up_coord.split(":")[1].split("-")[0]
            start = node_down_coord.split(":")[1].split("-")[1]
            
        return(chrom + ":" + start + "-" + end)

    
    rule ggsashmi:
        input:
            node = "Whippet/ggsashimi/{compare_name}/{gene}_{node}_{strand}.txt",
            tsv = "Whippet/ggsashimi/{compare_name}/{compare_name}.tvs",
            gtf = config["Gene_anontation_GTF"]
        params:
            region = lambda w: coord_to_region(w.gene, w.node, w.strand),
            out = "Whippet/ggsashimi/{compare_name}/{gene}_{node}_{strand}"
        output:
            "Whippet/ggsashimi/{compare_name}/{gene}_{node}_{strand}.pdf"
        shell:
            "python src/sashimi-plot.py -b {input.tsv} -c {params.region} -g {input.gtf} -o {params.out}"
            
    rule get_sashimis:
        input:
            expand("{sashimi}.pdf", sashimi=sashimis)
            
#chr:start-end
    

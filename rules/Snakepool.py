###### params #####

repeats = 10

###############



##### Single cell ####

import random
import csv
from collections import defaultdict

def partition (list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]



cluster_files = defaultdict(list)


#with open("/lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator_Final/Whippet/Tasic_clustering.txt") as Tasic:
with open(config["cluster_metadata"]) as Tasic:

    Tasic_clustering = csv.DictReader(Tasic, delimiter="\t")

    for row in Tasic_clustering:


        cluster_files[row[config["cluster_name"]]].append(row["Run_s"])




for compare_name, c in cluster_compare.items():


    g1, g2 = c

    c1_names = []
    for c1 in g1:

        c1_names += cluster_files[c1]

    c2_names = []
    for c2 in g2:
        c2_names += cluster_files[c2]


    np_A, np_B = cluster_compare_np[compare_name]


    #### these rules gereate a single indexed bam per condition which can be used for visualization


    rule :
        input:
            expand('Whippet/BAM/{sample}.sorted.bam', sample=c1_names)
        output:
            temp("Whippet/BAM/Merge/" + compare_name + ".A.bam")
        shell:
            "samtools merge  {output} {input}"

    rule :
        input:
            "Whippet/BAM/Merge/" + compare_name + ".A.bam"
        output:
            "Whippet/BAM/Merge/" + compare_name + ".A.sort.bam"
        shell:
            'samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output}'


    rule :
        input:
            expand('Whippet/BAM/{sample}.sorted.bam', sample=c2_names)
        output:
            temp("Whippet/BAM/Merge/" + compare_name + ".B.bam")
        shell:
            "samtools merge  {output} {input}"

    rule :
        input:
            "Whippet/BAM/Merge/" + compare_name + ".B.bam"
        output:
            "Whippet/BAM/Merge/" + compare_name + ".B.sort.bam"
        shell:
            'samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output}'

    ## Unpooled analysis



    rule:
        input:
            expand("Whippet/Quant/{sample}.psi.gz", sample=c1_names) + expand("Whippet/Quant/{sample}.psi.gz", sample=c2_names)
        output:
            "Whippet/Delta/Tasic/Unpooled/" + compare_name + ".run.sh"
        params:
            bin = config["whippet_bin_folder"],
            a = ",".join(expand("Whippet/Quant/{sample}.psi.gz", sample=c1_names)),
            b = ",".join(expand("Whippet/Quant/{sample}.psi.gz", sample=c2_names)),
            o = "Whippet/Delta/Tasic/Unpooled/" + compare_name
        shell:
            "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o} > {output}"


    rule:  #to avoid overload shell comandline
        input:
            "Whippet/Delta/Tasic/Unpooled/" + compare_name + ".run.sh"
        output:
            "Whippet/Delta/Tasic/Unpooled/" + compare_name + ".diff.gz"
        shell:
            "bash {input}"




    for r in range(repeats):


        c1_pools = partition(c1_names, np_A)
        c2_pools = partition(c2_names, np_B)

        p = 0

        target_pool_psi_A = []
        target_pool_psi_B = []

        delta_name = "Whippet/Delta/Tasic/" + compare_name +  "_rep_" +  str(r+1)

        #for pc1, pc2 in zip(c1_pools, c2_pools):


        for pc1 in c1_pools:

            p += 1



            FASTQ_c1 = [ config["fastq_path"] + x + ".fastq.gz" for x in  pc1 ]


            PSI_c1 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  pc1 ]

            pool_ID = "pool_" +str(r + 1) + "_"  + str(p)


            target_pool_psi_A.append("Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".psi.gz")




            rule:  #Quantification for A
                input:
                    fastq = FASTQ_c1,
                    index = "Whippet/Index/whippet.jls"
                output:
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".gene.tpm.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".isoform.tpm.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".jnc.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".map.gz",
                    "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID + ".psi.gz"
                params:
                    bin = config["whippet_bin_folder"],
                    output = "Whippet/Quant/Tasic/" + compare_name + "_A_" + pool_ID
                shell:
                    "julia {params.bin}/whippet-quant.jl <( cat {input.fastq} ) -x {input.index} --force-gz -o {params.output}"


        for pc2 in c2_pools:

            p += 1

            FASTQ_c2 = [ config["fastq_path"] + x + ".fastq.gz" for x in  pc2 ]
            PSI_c2 = [ "Whippet/Quant/" + x + ".psi.gz" for x in  pc2 ]

            pool_ID = "pool_" +str(r + 1) + "_"  + str(p)


            target_pool_psi_B.append("Whippet/Quant/Tasic/" + compare_name + "_B_" + pool_ID + ".psi.gz")



            rule:  #Quantification for B
                input:
                    fastq = FASTQ_c2,
                    index = "Whippet/Index/whippet.jls"
                output:
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".gene.tpm.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".isoform.tpm.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".jnc.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".map.gz",
                    "Whippet/Quant/Tasic/" +  compare_name + "_B_" + pool_ID + ".psi.gz"
                params:
                    bin = config["whippet_bin_folder"],
                    output = "Whippet/Quant/Tasic/" + compare_name + "_B_" + pool_ID
                shell:
                    "julia {params.bin}/whippet-quant.jl <( cat {input.fastq} ) -x {input.index} --force-gz -o {params.output}"



        rule:   #Diff
            input:
                target_pool_psi_A + target_pool_psi_B
            output:
                delta_name + ".diff.gz"

            params:
                bin = config["whippet_bin_folder"],
                a = ",".join(target_pool_psi_A),
                b = ",".join(target_pool_psi_B),
                o = delta_name
            shell:
                "julia {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"

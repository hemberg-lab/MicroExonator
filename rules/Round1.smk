
rule bwa_index:
    input:
        "Round1/ME_TAGs.fa"
    output:
        "Round1/ME_TAGs.fa.amb"
    conda:
        "../envs/core.yaml"
    shell:
        "bwa index {input}"

rule Round1_bwa_mem_to_tags:
    input:
        "Round1/ME_TAGs.fa",
        "FASTQ/{sample}.fastq.gz",
        "Round1/ME_TAGs.fa.amb"
    output:
        temp("Round1/{sample}.sam")
    threads: 5
    priority: 100
    params: 
        indel = config["indel_penalty"]
    resources:
        disk = 1
    conda:
        "../envs/core.yaml"
    shell:
        "bwa mem -t {threads} -O {params.indel} -L 25 {input[0]} {input[1]} | awk '$6 ~ /I/' > {output}"


rule Round1_alingment_pre_processing:
    input:
        "Round1/{sample}.sam"
    output:
        temp("Round1/{sample}.sam.pre_processed")
    priority: 100
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/alingment_pre_processing.py {input} F > {output}"

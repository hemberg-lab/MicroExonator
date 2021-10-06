rule Micro_Exon_Tags:
    input:
        "Round1/ME_TAGs.fa",
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric"
    output:
        "Round2/ME_canonical_SJ_tags.de_novo.fa"
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/Micro_exons_tags.py  {input} > {output}"

rule Get_ME_from_annotation:
    input:
        config["Genome_fasta"],
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric",
        config["Gene_anontation_bed12"],
        "data/GT_AG_U2_5.pwm",
        "data/GT_AG_U2_3.pwm",
        config["ME_DB"]
    params:
        bw = config["conservation_bigwig"],
        ME_len = config["ME_len"]
    output:
        "data/ME_canonical_SJ_tags.DB.fa",
        "data/DB.ME_centric"
    conda:
        "../envs/pybedtools.yaml"
    shell:
        "python2 src/Get_annotated_microexons.py  {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {params.bw} {params.ME_len} {input[5]} "


rule merge_tags:
    input:
        "Round2/ME_canonical_SJ_tags.de_novo.fa",
        "data/ME_canonical_SJ_tags.DB.fa"
    output:
        "Round2/ME_canonical_SJ_tags.fa"
    conda:
        "../envs/core.yaml"
    shell:
        "cat {input[0]} {input[1]} > {output}"


rule merge_ME_centric:
    input:
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric",
        "data/DB.ME_centric"
    output:
        "Round2/TOTAL.ME_centric.txt"
    conda:
        "../envs/core.yaml"
    shell:
        "cat {input[0]} {input[1]} > {output}"


rule Round2_bowtie_tags_index:
    input:
        "Round2/ME_canonical_SJ_tags.fa"
    output:
        "Round2/ME_canonical_SJ_tags.fa.1.ebwt"
    conda:
        "../envs/core.yaml"
    shell:
        "bowtie-build {input} {input}"

rule download_fastq2:
    input:
        "download/{sample}.download.sh",
        "Round2/TOTAL.ME_centric.txt"
    params:
        "FASTQ/{sample}.fastq"
    output:
        temp("FASTQ/round2/{sample}.fastq")
    priority: -10
    resources: 
        get_data = 1
    conda:
        "../envs/core.yaml"
    shell:
        #"bash {input[0]}"
        "bash {input[0]} && mv {params} {output}"

def hard_drive_behavior(fastq):
    if config.get("Optimize_hard_drive", False)=="T":
    
        if "validate_fastq_list" in config:
        
            to_validate = set[()]
            
            with open(config["validate_fastq_list"]) as fastq_list:
                reader = csv.reader(fastq_list, delimiter="\t")
                for row in reader:
                    to_validate.add(row[0])
                    
            if fastq in to_validate:
                return("FASTQ/round2/" + fastq + ".fastq.gz.valid")
            else:
                return(  "FASTQ/round2/" + fastq + ".fastq.gz")
                
        else:
            return(  "FASTQ/round2/" + fastq + ".fastq.gz")
    else:

        if "validate_fastq_list" in config:
        
            to_validate = set([])
            
            with open(config["validate_fastq_list"]) as fastq_list:
                reader = csv.reader(fastq_list, delimiter="\t")
                for row in reader:
                    to_validate.add(row[0])
                    
            if fastq in to_validate:
                return("FASTQ/" + fastq + ".fastq.gz.valid")
            else:
                return(  "FASTQ/" + fastq + ".fastq.gz")
        else:

            return("FASTQ/" + fastq + ".fastq.gz")


rule validate_fastq:
    input:
        "FASTQ/{sample}.fastq.gz"
    output:
        "FASTQ/{sample}.fastq.gz.valid"
    shell:
        "python3 src/validate_fastq.py {input}"
    
rule validate_fastq2:
    input:
        "FASTQ/round2/{sample}.fastq.gz"
    output:
        "FASTQ/round2/{sample}.fastq.gz.valid"
    shell:
        "python3 src/validate_fastq.py {input}"

rule Round2_bowtie_to_tags:
    input:
        "Round2/ME_canonical_SJ_tags.fa",
        hard_drive_behavior("{sample}"),
        "Round2/ME_canonical_SJ_tags.fa.1.ebwt"
    output:
        temp("Round2/{sample}.sam")
    threads: 5
    priority: 100
    conda:
        "../envs/core.yaml"
    shell:
        "gzip -dc {input[1]} |  bowtie {input[0]} -p {threads} -q - -S -v 2 --seed 123 | awk '!($6 ~ /I/) && !($6 ~ /D/) && !($6 ~ /S/) && !($6 ~ /*/)' > {output}"


rule Round2_alingment_pre_processing:
    input:
        "Round2/{sample}.sam"
    output:
        temp("Round2/{sample}.sam.pre_processed")
    priority: 100
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/alingment_pre_processing_round2_bowtie.py {input} F > {output}"

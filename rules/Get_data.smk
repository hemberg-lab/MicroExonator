if str2bool(config.get("Keep_fastq_gz", False)):
    rule download_fastq:
        input:
            "download/{sample}.download.sh"
        output:
            "FASTQ/{sample}.fastq.gz"
        resources:
            get_data = 1 
        conda:
            "../envs/core.yaml"
        priority: -10
        shell:
            "bash {input}"

else:
    rule download_fastq:
        input:
            "download/{sample}.download.sh"
        output:
            temp("FASTQ/{sample}.fastq.gz")
        resources:
            get_data = 1 
        conda:
            "../envs/core.yaml"
        priority: -10
        shell:
            "bash {input}"

rule unzip:
   input:
       "FASTQ/{sample}.fastq.gz"
   output:
       temp("FASTQ/{sample}.fastq")
   shell:
       "zcat {input} > {output}"

rule get_fastq:
    input:
        expand("FASTQ/{sample}.fastq.gz", sample=DATA)	        

if "Gene_anontation_bed12" in config:
    pass
else:
    rule generate_bed12:
        input:
            config["Gene_anontation_GTF"]
        output:
            "data/transcriptome.bed12"
        shell:
            "python2 src/GTFtoBED12.py {input} > {output}"
  
    config["Gene_anontation_bed12"] = "data/transcriptome.bed12"
            
        
rule generate_fasta_from_bed12:
    input:
        config["Genome_fasta"],
        config["Gene_anontation_bed12"]
    output:
        "data/transcripts.fa"
    conda:
        "../envs/pybedtools.yaml"
    shell:
        "python2 src/Get_fasta_from_bed12.py {input} > {output}"


rule Splice_Junction_Library:
    input:
        config["Genome_fasta"],
        "data/transcripts.fa",
        config["Gene_anontation_bed12"]
    params:
        ME_len = config["ME_len"]
    output:
        "Round1/ME_TAGs.fa"
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/SJ_tags_generator_for_micro_exons.py {input} {params.ME_len} > {output}"


rule GetPWM:
    input:
        config["Genome_fasta"],
        config["Gene_anontation_bed12"]
    params:
        config["GT_AG_U2_5"],
        config["GT_AG_U2_3"]
    output:
        "data/GT_AG_U2_5.pwm",
        "data/GT_AG_U2_3.pwm"
    conda:
        "../envs/core_py3.yaml"
    shell:
        "python3 src/Get_splicing_PWMs.py {input} {params} {output}"

#if str2bool(config.get("Only_whippet", False))==False:
#    rule gzip_fastq:
#        input:
#            "FASTQ/{sample}.fastq"
#        output:
#            temp("FASTQ/{sample}.fastq.gz")
#        priority: 100
#        shell:
#            "gzip -c {input} > {output}"
          
#else:
#    rule gzip_fastq:
#        input:
#            "FASTQ/{sample}.fastq"
#        output:
#            temp("FASTQ/{sample}.fastq.gz")
#        priority: 100
#        shell:
#            "gzip {input}"            
         


# rule sra_to_fastq:
#     input:
#         config["input_dir"] + "/{sample}.sra"
#     output:
#         temp("data/fastq_paired/{sample}.fastq")
#     shell:
#         "fastq-dump {input} -O data/fastq_paired/"


# rule fastq_gz_to_fastq:
#     input:
#         config["input_dir"] + "/{sample}.fastq.gz"
#     output:
#         temp("data/fastq/{sample}.fastq")
#     shell:
#         "gzip -dc {input} > {output}"
#
# rule fastq_input:
#     input:
#         config["input_dir"] + "/{sample}.fastq"
#     output:
#         "data/fastq/{sample}.fastq"
#     shell:
#         "ln -s {input} {output}"

#rule download_to_fastq:
#    input:
#        "download/{sample}.download.sh"
#    output:
#        "data/fastq/{sample}.fastq"
#    shell:
#        "bash {input}"


# rule split_fastq:
#     input:
#         "data/fastq_paired/{sample}.fastq"
#     output:
#         temp("data/fastq/{sample}.fastq")
#     shell:
#         "python2 src/split_paired_end.py {input} > {output}"

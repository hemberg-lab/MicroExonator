#rule GetPSI:
#    input:
#        "Round2/TOTAL.filter1.ME_SJ_coverage"
#    output:
#        expand("Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage.PSI", sample=DATA)
#    params:
#        path = "Round2/"
#    conda:
#        "../envs/core_py3.yaml"
#    shell:
#      "python3 src/GetPSI.py {input} 5 {params.path}"

#comparison_names = whippet_delta.keys()
      
#rule GetPSI:
#    input:
#        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
#    output:
#        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage.PSI"
#    conda:
#        "../envs/core_py3.yaml"
#    shell:
#      "python3 src/GetPSI.py {input} 5 > {output}"

rule get_GTF:
    input:
        genome_fasta = config["Genome_fasta"],
        annotation_bed12 = config["Gene_anontation_bed12"],
        annotation_GTF = config["Gene_anontation_GTF"],
        out_filtered_ME = "Report/out.high_quality.txt"
    params:
        chrM = False
    output:
        ME_GTF = "Report/out.high_quality.gtf"
    conda:
        "../envs/core.yaml"
    shell:
        "python src/get_isoforms2.py {input} {params}  > {output}"

rule get_GTF_robustly:
    input:
        genome_fasta = config["Genome_fasta"],
        annotation_bed12 = config["Gene_anontation_bed12"],
        annotation_GTF = config["Gene_anontation_GTF"],
        out_filtered_ME = "Report/out.robustly_detected.txt"
    params:
        chrM = False
    output:
        ME_GTF = "Report/out.robustly_detected.gtf"
    conda:
        "../envs/core.yaml"
    shell:
        "python src/get_isoforms2.py {input} {params}  > {output}"

filter_mode = config.get("filter_mode", "original")

def ME_GTF():
    if filter_mode=="unbiased":
        return("Report/out.robustly_detected.gtf")
    else:
        return("Report/out.high_quality.gtf")


if str2bool(config.get("Only_whippet", False)):
      rule whippet_index:
          input:
              Genome = config["Genome_fasta"],
              ME_GTF = config["Gene_anontation_GTF"]
          params:
              bin = config["whippet_bin_folder"],
              julia = config["julia"]
          output:
              "Whippet/Index/whippet.jls",
              "Whippet/Index/whippet.jls.exons.tab.gz"
          log:
              "logs/whippet_index.log"
          shell:
              "{params.julia} {params.bin}/whippet-index.jl --fasta {input.Genome} --gtf {input.ME_GTF} --index {output[0]} 2> {log}"
else:
      rule whippet_index:
          input:
              Genome = config["Genome_fasta"],
              ME_GTF = ME_GTF()
          params:
              bin = config["whippet_bin_folder"],
              julia = config["julia"]
          output:
              "Whippet/Index/whippet.jls",
              "Whippet/Index/whippet.jls.exons.tab.gz"
          log:
              "logs/whippet_index.log"
          shell:
              "{params.julia} {params.bin}/whippet-index.jl --fasta {input.Genome} --gtf {input.ME_GTF} --index {output[0]} 2> {log}" 



#if "Get_Bamfiles" not in config: #To make it optional
#      config["Get_Bamfiles"]="F"
      
if "whippet_flags" not in config:  
      config["whippet_flags"]=""
      

if str2bool(config.get("Get_Bamfiles", False)): #config["Get_Bamfiles"])==True:

      if "whippet_flags" in config:
            quant_flags = config["whippet_flags"]
      else:
            quant_flags = "--sam"

      rule  whippet_quant:
          input:
              #'FASTQ/{sample}.fastq.gz',
              hard_drive_behavior,
              "Whippet/Index/whippet.jls"
          params:
              bin = config["whippet_bin_folder"],
              output = "Whippet/Quant/{sample}",
              flags = "--sam",
              julia = config["julia"]
          output:
              "Whippet/Quant/{sample}.gene.tpm.gz",
              "Whippet/Quant/{sample}.isoform.tpm.gz",
              "Whippet/Quant/{sample}.jnc.gz",
              "Whippet/Quant/{sample}.map.gz",
              "Whippet/Quant/{sample}.psi.gz",
              sam = temp("Whippet/BAM/{sample}.sam")
          priority: 100
          shell:
              "{params.julia} {params.bin}/whippet-quant.jl {input[0]}  -x {input[1]} -o {params.output} {params.flags} > {output.sam}"

else:
      rule  whippet_quant:
          input:
              #'FASTQ/{sample}.fastq.gz',
              hard_drive_behavior,
              "Whippet/Index/whippet.jls"
          params:
              bin = config["whippet_bin_folder"],
              output = "Whippet/Quant/{sample}",
              flags = config["whippet_flags"],
              julia = config["julia"]
          output:
              temp("Whippet/Quant/{sample}.gene.tpm.gz"),
              temp("Whippet/Quant/{sample}.isoform.tpm.gz"),
              temp("Whippet/Quant/{sample}.jnc.gz"),
              temp("Whippet/Quant/{sample}.map.gz"),
              "Whippet/Quant/{sample}.psi.gz"
          priority: 100
          shell:
              "{params.julia} {params.bin}/whippet-quant.jl {input[0]}  -x {input[1]} -o {params.output} {params.flags}"      




rule SamToBam:
    input:
        sam = "Whippet/BAM/{sample}.sam"
    output:
        bam = temp("Whippet/BAM/{sample}.bam")
    shell:
        "samtools sort  -o {output.bam} -O BAM {input}"
        
rule BamIndex:
      input:
        bam = "Whippet/BAM/{sample}.bam"
      output:
        index = "Whippet/BAM/{sample}.bam.bai"
      shell:
        "samtools index {input.bam}"   

            
#expand("Whippet/BAM/{samples}.bam", samples=DATA)
#rule unzip_quant:
#    input:
#        "Whippet/Quant/{sample}.psi.gz"
#    output:
#        temp("Whippet/Quant/{sample}.psi")
#    shell:
#        "zcat {input} > {output}"



def get_downstream_PSI(wildcards):
    if str2bool(config.get("use_corrected_PSI", False)):
        if wildcards.sample in pe_samples:
            return( "Report/quant/corrected/PSI_sparse/bulk/pe/{sample}.corrected.PSI.gz" )
        else: 
            return( "Report/quant/corrected/PSI_sparse/bulk/se/{sample}.corrected.PSI.gz" )
    else:
        return( "Report/quant/{sample}.out_filtered_ME.PSI.uncorrected.gz"  )

rule ME_psi_to_quant:
    input:
        #ME_PSI = "Report/quant/{sample}.out_filtered_ME.PSI.uncorrected.gz",
        ME_PSI = get_downstream_PSI,
        whippet_PSI = "Whippet/Quant/{sample}.psi.gz"
    output:
        ME_psi =  "Whippet/Quant/{sample}.psi.ME.gz"
    shell:
        "python3 src/Replace_PSI_whippet2.py {input} {output} "


#rule gzip_ME_psi_to_quant:
#    input:
#        ME_psi =  "Whippet/Quant/{sample}.psi.ME"
#    output:
#        ME_psi_gz =  "Whippet/Quant/{sample}.psi.ME.gz"
#    shell:
#        "gzip {input}"

#condition1 = config["condition1"].split(",") 
#condition2 = config["condition2"].split(",")  
#comparison_name = config["comparison_name"]

# # # # # # #

def ME_list():
    if filter_mode=="unbiased":
        return("Report/out.robustly_detected.txt")
    else:
        return("Report/out.high_quality.txt")


if str2bool(config.get("Only_snakepool", False))==False:

      rule delta_ME_from_whippet:
          input:
              "Whippet/Index/whippet.jls.exons.tab.gz",
              "Whippet/Delta/{comparison_name}.diff.gz",
              #"Report/out.high_quality.txt"
              ME_list()
          output:
              "Whippet/Delta/{comparison_name}.diff.microexons"
          shell:
              "python src/whippet_delta_to_ME.py {input} > {output}"

      rule delta_ME_from_MicroExonator:
          input:
              "Whippet/Index/whippet.jls.exons.tab.gz",
              "Whippet/Delta/{comparison_name}.ME.diff.gz",
              #"Report/out.high_quality.txt"
              ME_list()
          output:
              "Whippet/Delta/{comparison_name}.diff.ME.microexons"
          shell:
              "python src/whippet_delta_to_ME.py {input} > {output}"

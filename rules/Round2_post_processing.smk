
rule ME_reads:
    input:
        "Round2/{sample}.sam.pre_processed",
        "FASTQ/{sample}.fastq.gz"
    output:
        temp("Round2/{sample}.sam.pre_processed.fastq")
    priority: 100
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/round2_ME_reads_fastq2.py {input}"
        
rule Get_Genome:
    input:
        config["Genome_fasta"]
    output:
        "data/Genome"
    priority: 100
    shell:
        "cp {input} {output}"

rule bowtie_genome_index:
    input:
        "data/Genome"
    output:
        "data/Genome" + ".1.ebwt"
    priority: 100
    threads : 8
    conda:
        "../envs/core.yaml"
    shell:
        "bowtie-build --threads {threads} {input} {input}"
        
if str2bool(config.get("skip_genome_alignment", False)):

    rule bowtie_to_genome:
        input:
            "Round2/{sample}.sam.pre_processed.fastq",
            "data/Genome",
            "data/Genome" + ".1.ebwt"
        output:
            temp("Round2/{sample}.sam.pre_processed.hg19.sam")
        priority: 100
        conda:
            "../envs/core.yaml"
        shell:
            "touch {output}"
else:

    rule bowtie_to_genome:
        input:
            "Round2/{sample}.sam.pre_processed.fastq",
            "data/Genome",
            "data/Genome" + ".1.ebwt"
        output:
            temp("Round2/{sample}.sam.pre_processed.hg19.sam")
        priority: 100
        conda:
            "../envs/core.yaml"
        shell:
            "bowtie {input[1]} -p 1 -q {input[0]} -S -v 2 --seed 123| awk '$2==0 || $2==16'> {output}"


rule Round2_filter:
    input:
        "Round2/{sample}.sam.pre_processed",
        "Round2/{sample}.sam.pre_processed.hg19.sam",
    output:
        temp("Round2/{sample}.sam.pre_processed.filter1")
    priority: 100
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/Filter1_round2.py {input} > {output}"


rule ME_SJ_coverage:
    input:
        "Round2/ME_canonical_SJ_tags.fa",
        "Round2/TOTAL.ME_centric.txt",
        config["Gene_anontation_bed12"],
        "Round2/{sample}.sam.pre_processed.filter1"
    params:
        ME_len = config["ME_len"]
    output:
        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
    priority: 100
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/ME_SJ_coverage.py {input} {params.ME_len} > {output}"

	
rule coverage_to_PSI_report:
    input:
	    "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
    params:
	    config["min_reads_PSI"],
	    config["paired_samples"]    
    output:
	    protected("Report/quant/{sample}.out_filtered_ME.PSI.gz")
    conda:
	    "../envs/core_py3.yaml"
    shell:
	    "python src/counts_to_PSI.py {input} {params} {output}"
	
	
rule Total_sample_exon_counts:
    input:
        expand("Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage", sample=DATA )
    output:
        "Round2/TOTAL.filter1.ME_SJ_coverage"
    conda:
        "../envs/core.yaml"
    shell:
      "cat Round2/*.filter1.ME_SJ_coverage > {output}"

rule write_ME_matches:
    input:
        "Round2/TOTAL.ME_centric.txt"
    output:
        "Round2/TOTAL.ME_centric.ME_matches.txt"
    conda:
        "../envs/core_py3.yaml"
    shell:
        "python3 src/Get_ME_matches.py {input} > {output}"

	
	
	

def get_min_reads():
    if 'min_reads_PSI' in config:
        return(int(config['min_reads_PSI']))
    else:
        return(5)


rule coverage_filter:
    input:
       #"Round2/TOTAL.filter1.ME_SJ_coverage"
       PSI_files = expand("Report/quant/{sample}.out_filtered_ME.PSI.gz", sample=DATA )
    params:
        min_reads_sample = get_min_reads()
    output:
        "Round2/TOTAL.sample_cov_filter.txt"
    script:
        "../src/coverage_sample_filter2.py"

def get_min_conservation():
    if "min_conservation" in config:
        return(int(config["min_conservation"]))
    else:
        return(2) #default value for min_conservation is 2
	
rule Output:
    input:
        ME_table = "Round2/TOTAL.ME_centric.txt",
        ME_coverage = "Round2/TOTAL.sample_cov_filter.txt",
        ME_matches_file = "Round2/TOTAL.ME_centric.ME_matches.txt"
    params:
        wd = config["working_directory"],
        min_number_files_detected = config["min_number_files_detected"],
        skip_mixture = str(str2bool(config.get("skip_mixture_model_filter", False))),
        min_conservation = get_min_conservation()
    output:
        out_filtered_ME = "Report/out_filtered_ME.txt",
        out_low_scored_ME = "Report/out_low_scored_ME.txt",
        out_shorter_than_3_ME = "Report/out_shorter_than_3_ME.txt",
        #"Report/report.html",
        #out_filtered_ME_cov = "Report/out_filtered_ME.cov.txt"
    log:
        "logs/Output.log"
    conda:
        "../envs/R.yaml"
    script:
        "../src/final_filters3.R"        
        
#    shell:
#        '''R -e  'rmarkdown::render("src/final_filters2.Rmd",params = list(ME_table="{params.wd}{input[0]}", ME_coverage="{params.wd}{input[1]}", ME_matches_file="{params.wd}{input[2]}", out_filtered_ME="{params.wd}{output[0]}", out_low_scored_ME="{params.wd}{output[1]}", out_shorter_than_3_ME="{params.wd}{output[2]}", min_number_files_detected={params.min_number_files_detected}, out_filtered_ME_cov="{params.wd}{output[4]}" ), output_file="{params.wd}{output[3]}")' 2> {log} '''


rule high_confident_filters:
    input:
        genome = config["Genome_fasta"],
        transcriptome = config["Gene_anontation_bed12"],
        #"Round2/TOTAL.filter1.ME_SJ_coverage",
        out_filtered = "Report/out_filtered_ME.txt",
        out_low_scored = "Report/out_low_scored_ME.txt",
	PSI_files = expand("Report/quant/{sample}.out_filtered_ME.PSI.gz", sample=DATA )
    output:
        high_qual = "Report/out.high_quality.txt"
    conda:
        "../envs/core_py3.yaml"
    script:
        "../src/high_confident_list.py"
	
#     shell:
#         "python src/high_confident_list.py {input}  > {output}"	
	
# if config.get("split_cov", False):


# 	rule coverage_to_PSI_sample:
# 	    input:
# 		    "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
# 	    params:
# 		    config["min_reads_PSI"],
# 		    config["paired_samples"]    
# 	    output:
# 		    "Round2/{sample}.out_filtered_ME.PSI.txt"
# 	    conda:
# 		    "../envs/core_py3.yaml"
# 	    shell:
# 		    "python src/counts_to_PSI.py {input} {params} > {output}"

			
# 	rule coverage_to_PSI_output:
# 	    input:
# 		    expand("Round2/{sample}.out_filtered_ME.PSI.txt", sample=DATA)
# 	    output:
# 		    "Report/out_filtered_ME.PSI.txt"
# 	    resources:
# 		    split = 1
# 	    priority: 20
# 	    shell:
# 		    "awk '(NR == 1) || (FNR > 1)' {input} > {output}"			

# else:
# 	rule coverage_to_PSI:
# 	    input:
# 		    "Round2/TOTAL.filter1.ME_SJ_coverage"
# 	    params:
# 		    config["min_reads_PSI"],
# 		    config["paired_samples"]    
# 	    output:
# 		    "Report/out_filtered_ME.PSI.txt"
# 	    conda:
# 		    "../envs/core_py3.yaml"
# 	    shell:
# 		    "python src/counts_to_PSI.py {input} {params} > {output}"


rule correct_quant:
    input:
        quant = "Report/quant/{sample}.out_filtered_ME.PSI.gz",
        ME_centric = "Round2/TOTAL.ME_centric.txt",
        spanning_ME_reads =  "Round2/ME_reads/{sample}.ME_spanning_reads.tsv"
    output:
        corrected_quant = "Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz",
        count_spanning_ME_reads = =  "Round2/ME_reads/{sample}.tsv"
    script:
        "../src/correct_quant.py.py"
		

rule get_all_corrected_quant:
    input:
        expand("Report/quant/corrected/{sample}.out_filtered_ME.PSI.gz", sample=DATA)

rule annotation_stats:
    input:
        config["Gene_anontation_bed12"],
        "Report/out.high_quality.txt",
    params:
        30
    output:
        "Report/stats/Microexons.not_consensus",
        "Report/stats/Microexons.annotation.stats"
    conda:
        "../envs/core_py3.yaml"
    shell:
        "python3 src/stats/discovery_stats.py {input} {params}"

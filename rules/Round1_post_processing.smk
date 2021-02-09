

rule row_Micro_Exon_reads:
    input:
        config["Genome_fasta"],
        "Round1/{sample}.sam.pre_processed",
	"FASTQ/{sample}.fastq.gz"
    output:
        temp("Round1/{sample}.sam.row_ME"),
        temp("Round1/{sample}.sam.row_ME.fastq")
    conda:
        "../envs/core.yaml"
    shell:
        "python2 src/row_ME2.py {input} > {output[0]}"


rule hisat2_genome_index:
    input:
        config["Genome_fasta"]
    output:
        "data/Genome.1.ht2"
    threads: 5
    conda:
        "../envs/core.yaml"
    shell:
        "hisat2-build {input} data/Genome"

if str2bool(config.get("skip_genome_alignment", False)):

	rule hisat2_to_Genome:
	    input:
	        "Round1/{sample}.sam.row_ME.fastq",
	        "data/Genome.1.ht2"
	    output:
	        temp("Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam")
	    threads: 1
	    conda:
	        "../envs/core.yaml"
	    shell:
	        "touch {output}"
else:

	rule hisat2_to_Genome:
	    input:
	        "Round1/{sample}.sam.row_ME.fastq",
	        "data/Genome.1.ht2"
	    output:
	        temp("Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam")
	    threads: 1
	    conda:
	        "../envs/core.yaml"
	    shell:
	        "hisat2 -x data/Genome -U {input[0]} > {output}"


rule Round1_filter:
    input:
        config["Genome_fasta"],
        "Round1/{sample}.sam.row_ME",
        "Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam",
        "data/GT_AG_U2_5.pwm",
        "data/GT_AG_U2_3.pwm"
    params:
        bw = config["conservation_bigwig"],
        ME_len = config["ME_len"]
    output:
        protected("Round1/{sample}.sam.row_ME.filter1")
    conda:
        "../envs/pybedtools.yaml"
    shell:
        "python2 src/ME_filter1.py {input} {params.bw} {params.ME_len} > {output}"


rule Micro_Exon_table:
    input:
        expand("Round1/{sample}.sam.row_ME.filter1", sample=DATA )
    output:
        protected("Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric")
    conda:
        "../envs/core.yaml"
    shell:
        "cat Round1/*.sam.row_ME.filter1 | awk 'NF==16' > Round1/TOTAL/TOTAL.sam.row_ME.filter1  &&"
        "python2 src/ME_centric_table.py Round1/TOTAL/TOTAL.sam.row_ME.filter1 > {output}"


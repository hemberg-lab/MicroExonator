
############## Gene Count #######

# rule generate_star_olego:
#     input:
#         "config["Genome_fasta"]",
#
#     shell:
#         "start --runThreadN 5 --runMode genomeGenerate --genomeDir data/ "
#
#
#
# rule generate_star_index:
#     input:
#         "config["Genome_fasta"]",
#
#     shell:
#         "start --runThreadN 5 --runMode genomeGenerate --genomeDir data/ "
#



rule total_hisat2_to_genome:
     input:
         "FASTQ/{sample}.fastq",
         "data/Genome.1.ht2"
     output:
         "Genome_aligments/Hisat2/{sample}.sam"
     threads: 5
     shell:
         "hisat2 -x data/Genome -U {input[0]} -p 5 > {output}"

rule total_olego_to_Genome:
    input:
        "FASTQ/{sample}.fastq"
    output:
        "Genome_aligments/Olego/{sample}.sam"
    threads: 10
    shell:
        "/lustre/scratch117/cellgen/team218/gp7/olego/olego  -t 10 data/Genome_olego {input} > {output}"


rule total_STAR_to_Genome:
    input:
        "FASTQ/{sample}.fastq"
    output:
        "Genome_aligments/STAR/{sample}.samAligned.out.sam"
    threads: 5
    shell:
        "STAR --genomeDir data --readFilesIn {input} --runThreadN 5 --outFileNamePrefix  {output}"

rule mv_STAR:
    input:
        "Genome_aligments/STAR/{sample}.samAligned.out.sam"
    output:
        "Genome_aligments/STAR/{sample}.sam"
    threads: 5
    shell:
        "mv {input} {output}"




rule total_tophat_to_Genome:
    input:
        "FASTQ/{sample}.fastq"
    output:
        dir = "Genome_aligments/Tophat2/{sample}",
        sam = "Genome_aligments/Tophat2/{sample}.sam"
    threads: 5
    shell:
        "tophat2 -p 5 --no-convert-bam --microexon-search -o {output.dir} data/Genome_bowtie2 {input} && mv {output.dir}/accepted_hit.sam {output}"


rule SJ_count:
    input:
        "Genome_aligments/{Software}/{sample}.sam"
    output:
        "Genome_aligments/{Software}/{sample}.sam.SJ_count"
    shell:
        "python2 src/Get_introns_from_sam.py {input} Rd1 40 1000000 8 > {output}"


rule sam_merge:
    input:
        ["Genome_aligments/{Software}/" + x for x in expand("{sample}.sam", sample=DATA ) ]
    output:
        temp("Genome_aligments/{Software}/TOTAL.sam")
    shell:
        "samtools merge {output} {input}"


rule get_exons:
    input:
        "Genome_aligments/{Software}/TOTAL.sam"
    output:
        "Genome_aligments/{Software}/TOTAL.exons.{Software}"
    shell:
        "python2 Get_exons_from_sam.py {input} > {output}"




rule SJ_ground_count:
    input:
        config["fastq_path"] + '{sample}.fastq.gz'
    output:
        "Ground_Truth/{sample}.GT.SJ_count"
    shell:
        "python2 SJ_count_truth.py /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.bed12 simulated_ME_isoforms.bed12 {input}  > {output}"



rule gene_count:
    input:
        "/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.gtf",
        "Genome_aligments/{sample}.sam"
    output:
        "Genome_aligments/{sample}.gene_count.txt"
    threads: 1
    shell:
        "featureCounts -a {input[0]} -o {output} {input[1]}"



rule done_gene_count:
    input:
        expand("Genome_aligments/{sample}.gene_count.txt", sample=DATA )
    output:
        "Round2/done.txt"
    shell:
        "echo done > {output}"
#####

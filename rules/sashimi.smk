import csv

def get_sasshimi_window(ME, w):

	ME_strand, ME_start, ME_end = ME.split("_")[-3:]
	ME_chrom =  "_".join(ME.split("_")[:-3])
					
	start = str(int(ME_start)-w)
	end = str(int(ME_start)+w)
	return(chrom + ":" + start + "-" + end)

if "sashimi_tsv" in config:

    target_ME = set([])
	
    with open("Report/novel_highly_included.tsv") as tsv:
					 
        reader = csv.DictReader()
        for row in reader:
            if float(row["mean_PSI"])>0.9:
                target_ME.add(row["ME"])

        rule ggsashmi_bulk_scripts:
	    #input:
		#node = "ggsashimi/{ME}.txt",
		#gtf = config["Gene_anontation_GTF"],
                #tsv = config["sashimi_tsv"]
		#bams = lambda w: expand("Whippet/BAM/Merge/{cluster}.sort.bam",  cluster=compare_all_clusters[w.compare_name])
            params:
                gtf = config["Gene_anontation_GTF"],
                tsv = config["sashimi_tsv"],
                region = lambda w: get_sasshimi_window(w.ME, 1000),
                out = "ggsashimi/{ME}"
            output:
                "ggsashimi/{ME}.pdf"
            shell:
                "python src/sashimi-plot.py -b {params.tsv} -c {params.region} -g {params.gtf} -o {params.out}"

        rule get_sashimis_bulk:
            input:
                expand("ggsashimi/{ME}.pdf", ME=target_ME)

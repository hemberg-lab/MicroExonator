import csv

def get_sasshimi_window(ME, w):

	ME_strand, ME_start, ME_end = ME.split("_")[-3:]
	chrom =  "_".join(ME.split("_")[:-3])
					
	start = str(int(ME_start)-w)
	end = str(int(ME_start)+w)
	return(chrom + ":" + start + "-" + end)

if "sashimi_tsv" in config:

    target_ME = set([])
	
    with open("Report/novel_highly_included.tsv") as tsv:
					 
        reader = csv.DictReader(tsv, delimiter="\t")
        for row in reader:
            #print(row)
            if float(row["mean_PSI"])>0.9:
                target_ME.add(row["ME"])

        rule ggsashmi_bulk_scripts:
            params:
                gtf = config["Gene_anontation_GTF"],
                tsv = config["sashimi_tsv"],
                region = lambda w: get_sasshimi_window(w.ME, 10000),
                out = "ggsashimi/{ME}",
                pallete = config["sashimi_pallete"]
            output:
                "ggsashimi/{ME}.sh"
                #"ggsashimi/{ME}.pdf"
            shell:
                "echo python src/sashimi-plot.py -b {params.tsv} -c {params.region} -g {params.gtf} -o {params.out} -P {params.pallete} -C 3 -O 3 -A mean > {output}"

        rule run_sashimi:
            input:
                "ggsashimi/{ME}.sh"
            output:
                "ggsashimi/{ME}.pdf"
            shell:
                "bash {input}"

        rule get_sashimis_bulk:
            input:
                expand("ggsashimi/{ME}.pdf", ME=target_ME)

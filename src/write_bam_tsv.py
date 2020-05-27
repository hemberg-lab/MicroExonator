from snakemake.utils import min_version
import csv


with open(snakemake.input[0]) as run:   #Populating the dictionaries

    run_metadata = csv.DictReader(run, delimiter="\t")

    for row in run_metadata:
        
        with open("Whippet/ggsashimi/" + row["Compare_ID"] + "/" + row["Compare_ID"] + ".tvs", "w") as out:

            A_cluster_names = []
            B_cluster_names = []
            
            for c in row["A.cluster_names"].split(","):
                
                out.write(c + "\t" + "Whippet/BAM/Merge/" + c + ".sort.bam" + "\n")

            for c in row["B.cluster_names"].split(","):

                out.write(c + "\t" + "Whippet/BAM/Merge/" + c + ".sort.bam" + "\n")

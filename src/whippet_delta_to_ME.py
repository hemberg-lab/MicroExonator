import csv
from collections import defaultdict

def main(jls_exons_tab, delta ):

    node_exons = dict()


    with open(jls_exons_tab) as F:

        reader = csv.DictReader(F, delimiter="\t")
        
        for row in reader:
            chrom, locus, strand = row["Potential_Exon"].split(":")
            estart, eend = locus.split("-")

            for node in row["Whippet_Nodes"].split(","):

                node_exons[(row["Gene"], node)] = [row["Potential_Exon"], row["Is_Annotated"]]


    with open(delta) as F: 

        reader = csv.DictReader(F, delimiter="\t")

        for row in reader:

            if (row["Gene"], row["Node"] ) in node_exons:

                out =  [row[x] for x in header] + node_exons[(row["Gene"], row["Node"] )]

                TOTAL_DIFF.writerow(out)

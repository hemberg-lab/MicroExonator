import csv, sys
from collections import defaultdict
import gzip



def main(jls_exons_tab, delta, high_qual_ME ):

    node_exons = dict()
    MEs = set([])
    ME_info = dict()
    
    header = ["Gene", "Node", "Coord", "Strand", "Type", "Psi_A.mean", "Psi_B.mean", "DeltaPsi.mean", "DeltaPsi.sd", "Probability.mean", "Probability.sd", "Probability.var", "N.detected.reps", "cdf.beta", "is.diff", "microexon_ID"]


    print("\t".join(["exon_ID"] + header))


    with open(high_qual_ME) as F:

        reader = csv.DictReader(F, delimiter="\t")

        for row in reader:
            MEs.add(row["ME"])
    
    
    with gzip.open(jls_exons_tab, mode="rt") as F:

        reader = csv.DictReader(F, delimiter="\t")
        
        for row in reader:
            chrom, locus, strand = row["Potential_Exon"].split(":")
            estart, eend = locus.split("-")

            for node in row["Whippet_Nodes"].split(","):
                node_exons[(row["Gene"], node)] = [row["Potential_Exon"], row["Is_Annotated"]]


    with open(delta) as F: 

        reader = csv.DictReader(F, delimiter="\t")

        for row in reader:
            
            chrom, pos = row["Coord"].split(":")
            estart, eend = pos.split("-")
            estart = str(int(estart)-1)
            exon_ID = "_".join([chrom, row["Strand"], estart, eend])
            
            #if exon_ID == "chr10_+_127272438_127272444":
            #if "12727243" in exon_ID: 
            #    print(row, exon_ID)
            


            if (row["Gene"], row["Node"] ) in node_exons:
                

                
                Potential_Exon, Is_Annotated = node_exons[(row["Gene"], row["Node"] )]
                
                out =  [row[x] for x in header] + node_exons[(row["Gene"], row["Node"] )]
 

                
                if row["Type"]=="AD":
            
                    nchrom, nstrand, nstart, nend = exon_ID.split("_")


                    echrom, eloci, estrand =  Potential_Exon.split(":")

                    estart, eend =  eloci.split("-")

                    if estrand == "+" and eend == nend:

                        new_exon_ID = "_".join([echrom, estrand, str(int(estart)-1), eend ])

                        exon_ID = new_exon_ID

                    if estrand == "-" and str(int(estart)-1) == nstart:

                        new_exon_ID = "_".join([echrom, estrand, str(int(estart)-1), eend ])

                        exon_ID = new_exon_ID


                elif row["Type"]=="AA":

                    nchrom, nstrand, nstart, nend = exon_ID.split("_")


                    echrom, eloci, estrand =  Potential_Exon.split(":")

                    estart, eend =  eloci.split("-")

                    if estrand == "-" and eend == nend:

                        new_exon_ID = "_".join([echrom, estrand, str(int(estart)-1), eend ])


                        exon_ID = new_exon_ID


                    if estrand == "+" and str(int(estart)-1) == nstart:

                        new_exon_ID = "_".join([echrom, estrand, str(int(estart)-1), eend ])

                        exon_ID = new_exon_ID
                        
                        #if "12727243" in exon_ID: 
                        #    print(row, exon_ID, Potential_Exon)
           

                if exon_ID in MEs:
                   print("\t".join( [row[x] for x in header] + [exon_ID] ))
          
                else:
                   print("\t".join( [row[x] for x in header] + [NA] ))
            
            
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])

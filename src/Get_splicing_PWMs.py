import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from shutil import copyfile

Genome = {}

def Genomictabulator(fasta):

    f = open(fasta)

    for chrfa in SeqIO.parse(f, "fasta"):
        Genome[chrfa.id] = chrfa.seq


    f.close()

def main(bed12, in_GT_AG_U2_5, in_GT_AG_U2_3, out_GT_AG_U2_5, out_GT_AG_U2_3):

    if in_GT_AG_U2_5=="NA" and in_GT_AG_U2_3=="NA":


        GT_AG_U2_5 = defaultdict(int)   #intro-centric
        GT_AG_U2_3 = defaultdict(int)

        for row in csv.reader(open(bed12), delimiter = '\t'):

                csv.field_size_limit(1000000000)

                qstarts = list(map (int, row[11].strip(",").split(",")))[1:-1]
                blocksizes = list(map(int, row[10].strip(",").split(",")))[1:-1]

                start = int(row[1])
                strand = row[5]
                bn = int(row[9])
                chrom = row[0]
                
                if chrom in Genome:

                    for q1, b in zip(qstarts, blocksizes):
                        estart = start + q1
                        eend = start + q1 + b
                        elenght = eend - estart


                        ME5 = str(Genome[chrom][estart-14:estart+3]).upper()  #exon-centric
                        ME3 = str(Genome[chrom][eend-3:eend+10]).upper()


                        if strand == "-":

                            ME5 = str(Genome[chrom][eend-3:eend+14].reverse_complement()).upper()
                            ME3 = str(Genome[chrom][estart-10:estart+3].reverse_complement()).upper()


                        dn = ME3[3:5] + ME5[-5:-3]


                        if dn=="GTAG":

                            for pos, nt in enumerate(ME3):

                                GT_AG_U2_5[(pos, nt)] += 1

                            for pos, nt in enumerate(ME5):

                                GT_AG_U2_3[(pos, nt)] += 1



        with open(out_GT_AG_U2_5, "w") as GT_AG_U2_5_out:


            GT_AG_U2_5_out.write( "\t".join(["A", "C", "G", "T"]) +"\n")

            for i in range(len(GT_AG_U2_5)): #This range is about 4 times biger
                A = GT_AG_U2_5[(i, "A")]
                G = GT_AG_U2_5[(i, "G")]
                C = GT_AG_U2_5[(i, "C")]
                T = GT_AG_U2_5[(i, "T")]

                TOTAL = A + G + C + T

                if TOTAL >0:
                    GT_AG_U2_5_out.write("\t".join(map(str, [ x/TOTAL for x in [A, C, G, T]])) + "\n" )



        with open(out_GT_AG_U2_3, "w") as GT_AG_U2_3_out:


            GT_AG_U2_3_out.write( "\t".join(["A", "C", "G", "T"]) +"\n")

            for i in range(len(GT_AG_U2_3)): #This range is about 4 times biger
                A = GT_AG_U2_3[(i, "A")]
                G = GT_AG_U2_3[(i, "G")]
                C = GT_AG_U2_3[(i, "C")]
                T = GT_AG_U2_3[(i, "T")]

                TOTAL = A + G + C + T

                if TOTAL >0:
                    GT_AG_U2_3_out.write("\t".join(map(str, [ x/TOTAL for x in [A, C, G, T]])) + "\n" )


    else:

        copyfile(in_GT_AG_U2_5, out_GT_AG_U2_5)
        copyfile(in_GT_AG_U2_3, out_GT_AG_U2_3)




if __name__ == '__main__':
    Genomictabulator(sys.argv[1])
    main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

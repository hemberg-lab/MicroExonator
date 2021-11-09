import sys
import csv
import gzip
csv.field_size_limit(300000000000000)


def main(ME_PSI, whippet_PSI, output):

    Coord_info = dict()

    with gzip.open(ME_PSI, mode="rt") as file:

        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:
            
            if "ME_coords" in row:
                ME = row["ME_coords"]
            else:
                ME = row["ME"]
            
            ME_strand, ME_start, ME_end = ME.split("_")[-3:]  #to avoid errors with chromosomes like chr1_GL456210_random
            ME_chrom = "_".join(ME.split("_")[:-3])

            Coord = ME_chrom + ":" + str(int(ME_start)+1) + "-" + ME_end

            info = ME, Coord, row['PSI'], row['CI_Lo'], row['CI_Hi']

            Coord_info[Coord] = info

    header = [ 'Gene', 'Node', 'Coord', 'Strand', 'Type', 'Psi', 'CI_Width', 'CI_Lo,Hi', 'Total_Reads', 'Complexity', 'Entropy', 'Inc_Paths', 'Exc_Paths', 'Edges'  ]

    with gzip.open(whippet_PSI, mode="rt") as whippet_PSI, gzip.open(output, mode="wt") as out:

        out.write("\t".join(header) + "\n" )

        reader = csv.DictReader(whippet_PSI, delimiter="\t")

        for row in reader:

            whippet_row = [row[x] for x in header ]

            if row['Coord'] in Coord_info:

                ME, Coord, ME_PSI, CI_Lo, CI_Hi = Coord_info[row['Coord']]

                if ME_PSI == "NA":

                    out.write("\t".join(whippet_row) + "\n" )

                else:

                    ME_CI_Lo_Hi = ",".join([ CI_Lo, CI_Hi])

                    ME_CI_Width = str(float(CI_Hi) - float(CI_Lo))

                    ME_out = [row['Gene'], row['Node'], Coord, row['Strand'], row['Type'], ME_PSI, ME_CI_Width, ME_CI_Lo_Hi, row['Total_Reads'], row['Complexity'], row['Entropy'], row['Inc_Paths'], row['Exc_Paths'], row['Edges']]

                    out.write("\t".join(ME_out) + "\n")

            else:

                out.write("\t".join(whippet_row) + "\n")


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])

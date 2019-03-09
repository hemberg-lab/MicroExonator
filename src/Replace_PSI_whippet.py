import sys
import csv
csv.field_size_limit(300000000000000)

def main(ME_PSI, whippet_PSI):

    Coord_info = dict()

    reader1 = csv.reader(open(ME_PSI), delimiter="\t")

    header1 = next(reader1)

    for row in reader1:


        ME, Coord, PSI, CI_Lo, CI_Hi, Class = row

        Coord_info[Coord] = row

    reader2 = csv.reader(open(whippet_PSI), delimiter="\t")
    header2 = next(reader2)

    print( "\t".join(header2) )

    for row in reader2:

        Gene, Node, Coord, Strand, Type, Psi, CI_Width, CI_Lo_Hi, Total_Reads, Complexity, Entropy, Inc_Paths, Exc_Paths, Edges = row

        if Coord in Coord_info:

            ME, Coord, ME_PSI, CI_Lo, CI_Hi, Class = Coord_info[Coord]

            if ME_PSI == "NA":

                print("\t".join(row))

            else:

                ME_CI_Lo_Hi = ",".join([ CI_Lo, CI_Hi])

                ME_CI_Width = str(float(CI_Hi) - float(CI_Lo))

                out = [Gene, Node, Coord, Strand, Type, ME_PSI, ME_CI_Width, ME_CI_Lo_Hi, Total_Reads, Complexity, Entropy, Inc_Paths, Exc_Paths, Edges]

                print("\t".join(out))

        else:

            print("\t".join(row))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2]  )

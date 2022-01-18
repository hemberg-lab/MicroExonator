import sys
import csv
import gzip
from snakemake.utils import min_version

csv.field_size_limit(100000000)
csv.field_size_limit()



def main(mode, out_file, file_list  ):
    with gzip.open(out_file, 'wt') as out:

        for file in file_list:

            with gzip.open(file, mode="rt") as f:

                if mode=="Isoform" or mode=="Gene":
                    header = ["Sample", mode, "TpM", "Read_Counts"]
                elif mode=="PSI":
                    header = ['Sample', 'Gene', 'Node', 'Coord', 'Strand', 'Type', 'Psi', 'CI_Width', 'CI_Lo,Hi', 'Total_Reads', 'Complexity', 'Entropy', 'Inc_Paths', 'Exc_Paths', 'Edges']

                writer = csv.DictWriter(out, fieldnames=header, extrasaction='ignore', delimiter="\t")
                writer.writeheader()

                sample = file.strip(snakemake.params["trim"])
                #sample = file.strip(file)

                #sample = "".join(file.split("/")[-1].split(".")[:-2])  #files needs to finish with *.fastq.gz
                reader = csv.DictReader(f, delimiter="\t")

                for row in reader:

                    row["Sample"] = sample
                    writer.writerow(row)

#print(snakemake.input)

if __name__ == '__main__':
    main(snakemake.params["feature"], snakemake.output["merged"],  snakemake.input["files"])

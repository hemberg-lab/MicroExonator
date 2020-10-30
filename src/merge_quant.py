import sys
import csv
import gzip


csv.field_size_limit(100000000)
csv.field_size_limit()


def main(mode, out_file, file_list  ):
    with open(out_file, 'w') as out:

        for file in file_list:

            with gzip.open(file, mode="rt") as f:

                header = ["Sample", mode, "TpM", "Read_Counts"]

                writer = csv.DictWriter(out, fieldnames=header, extrasaction='ignore', delimiter="\t")
                writer.writeheader()

                sample = file.split("/")[-1].split(".")[0]
                reader = csv.DictReader(f, delimiter="\t")

                for row in reader:

                    row["Sample"] = sample
                    writer.writerow(row)
                
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3:])

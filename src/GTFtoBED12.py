import sys
import csv
from collections import defaultdict

def main(gtf_file):
    
    with open(gtf_file) as gtf:
        
        transcript_coords = dict()
        transcript_qstarts_blocksize = defaultdict(list)
        reader = csv.reader(gtf, delimiter="\t")
        
        for row in reader:
            
            if row[0][0]!="#":
                
                chrom = row[0]
                group = row[1]
                blocktype = row[2]
                block_start = int(row[3]) - 1
                block_end = int(row[4])
                block_size = block_end - block_start
                strand = row[6]

                tags = row[8].strip(" ").split(";")

                for t in tags:
                    pair =  t.strip(" ").split(" ")
                    if pair!=['']:
                        ID_type, ID  = pair
                        if ID_type == "transcript_id":
                            transcript = ID.strip('"')

                if blocktype == 'transcript':

                    transcript_coords[transcript] = (chrom, block_start, block_end, strand)

                if blocktype == 'exon':

                    exon_size = block_end - block_start


                    transcript_qstarts_blocksize[transcript].append((block_start, exon_size))


        for trancript in transcript_coords:


            chrom, start, end, strand = transcript_coords[trancript]

            n_blocks = len(transcript_blocksizes[trancript])


            q_b_tuples = sorted(transcript_qstarts_blocksize[trancript] , key=lambda x: x[0])

            qstarts_list = [x[0] for x in q_b_tuples ]
            blocksizes_list = [x[1] for x in q_b_tuples ]

            qstarts = ",".join(map(str, [x - start for x in qstarts_list] ))
            blocksizes = ",".join(map(str, blocksizes_list))


            bed12 = [chrom, start, end, transcript, "0", strand, start, end, "0", n_blocks, blocksizes, qstarts]

            print( "\t".join(map(str, bed12)))

if __name__ == '__main__':
	main(sys.argv[1])

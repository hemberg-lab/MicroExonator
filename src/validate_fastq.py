from Bio.SeqRecord import SeqRecord
import gzip
from Bio import bgzf

def main(row_fastq):
	
    with gzip.open(row_fastq, "rt") as f, gzip.open(row_fastq + ".gz", 'wt') as out:

        for read in SeqIO.parse(f, "fastq"):
 
            fastq_out = SeqRecord( read.seq, read.id, description = "" )
            fastq_out.letter_annotations["phred_quality"] = read.letter_annotations["phred_quality"][ :len(read.seq)]
            
            if len(read.seq)==len(read.letter_annotations["phred_quality"]):
                out.write(fastq_out.format("fastq"))
                
                
if __name__ == '__main__':
    main(sys.argv[1])

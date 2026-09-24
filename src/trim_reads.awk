# Trim FASTQ reads (and qualities) to L nt and record a read-length histogram.
# Usage: awk -v L=150 -v out=lengths.tsv -f trim_reads.awk reads.fastq
NR % 4 == 2 { hist[length($0)]++ }
NR % 2 == 0 && length($0) > L { $0 = substr($0, 1, L) }
{ print }
END {
    print "length\treads" > out
    for (l in hist) print l "\t" hist[l] > out
}

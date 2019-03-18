import sys,csv,pdb

##### Use as #####
    #python bed12 MicroExonator_high_qual ME_len
    #Compatible with python2 and python3
######

#def main(bed12, ME_high_cual , ME_len):

def bed12_microexon_extract(bed12_path,ME_len):
	ME_annotated = set()
	for row in csv.reader(open(bed12_path), delimiter = '\t'):
		blocksizes = list(map(int, row[10].strip(",").split(",")))
		qstarts = list(map(int, row[11].strip(",").split(",")))

		start = int(row[1])
		end = int(row[2])
		strand = row[5]
		bn = int(row[9])
		chrom = row[0]

		for q1, q2, q3, b1, b2, b3 in zip(qstarts, qstarts[1:] , qstarts[2:], blocksizes, blocksizes[1:], blocksizes[2:]):

			estart = start + q2
			eend = start + q2 + b2
			elength = eend - estart
			exon = (chrom, strand, str(estart), str(eend))
			transcript = row[3]
			SJ_start = start + q1 + b1
			SJ_end = start + q3
			ME_intron = " ".join([chrom, str(SJ_start), str(SJ_end), "SJ", "0", strand])

			if elength <= ME_len:
               			ME_annotated.add(exon)

	return ME_annotated

def microexonator_output_reader(high_qual_filters):
	ME_MicroExonator = set()
	reader = csv.reader(open(high_qual_filters), delimiter = '\t')
	next(reader)
	for row in reader:
		ME = row[0]
		chrom = "_".join(ME.split("_")[:-3])  # To avoid errors with chrom_
		strand, start, end = ME.split("_")[-3:]
		ME = (chrom, strand, start, end)
		ME_MicroExonator.add(ME)
	return ME_MicroExonator

def compare(set_Microexonator,set_consensus):
	set_novel = set_Microexonator-set_consensus
	datafile=open("Microexons.not_consensus","w")
	for line in set_novel:
		datafile.write('\t'.join([str(x) for x in line])+'\n')
	datafile.close()

	set_total = set_Microexonator+set_consensus
	datafile=open("Microexons.annotation.stats","w")
	for line in set_total:
		if line in set_novel:
			datafile.write('\t'.join([str(x) for x in line])+'\t'+str(1)+'\n')
		else:
			datafile.write('\t'.join([str(x) for x in line])+'\t'+str(0)+'\n')
	datafile.close()
	return


if __name__ == '__main__':
	consensus_file = sys.argv[1]
	MicroExonator_file = sys.argv[2]
	length_limit =int(sys.argv[3]) 
	Consensus_microexons = bed12_microexon_extract(consensus_file,length_limit)
	MicroExonator_microexons=microexonator_output_reader(MicroExonator_file)
	compare(MicroExonator_microexons,Consensus_microexons)

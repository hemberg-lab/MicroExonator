import sys
import csv
from collections import defaultdict



keys = set([])
exon_start = []
exon_size = []
transcript_start_end = {}
transcript_id_info = {}

def converter (GTF):
	#reader = csv.reader(open(GTF), delimiter = '\t')

	f = open(GTF)
	for line in f:

	    if line[0]!="#":
                    
                row = line.replace(";","\t").replace(" ","\t").replace("\t\t","\t").split("\t")
                chr = row[0]
                group = row[1]
                blocktype = row[2]
                block_start = int(row[3]) - 1
                block_end = int(row[4])
                block_size = block_end - block_start
                strand = row[6]
                        
                if blocktype == 'transcript':

                    #transcript_start_end[transcript_id] = (block_start, block_end)
                    try:
                        gene_id = row[9]
                        transcript_id = row[11].strip('"')
                        transcript_start_end[transcript_id] = (block_start, block_end)
                    except IndexError:
                        pass
                    except ValueError:
                        pass

                if blocktype == 'exon':
                    keys.add(transcript_id)
                    transcript_id_info[transcript_id] = (chr, group, blocktype, strand, gene_id)
                    exon_start.append((transcript_id, block_start))
                    exon_size.append((transcript_id, block_size))


	#dict_id_info = dict(transcript_id_info)
	#dict_transcript_start_end = dict(transcript_start_end)

	dict_exon_start = defaultdict(list)
	for transcript_id, block_start in exon_start:
		dict_exon_start[transcript_id].append(block_start)

	dict_exon_size = defaultdict(list)
	for transcript_id, block_size in exon_size:
		dict_exon_size[transcript_id].append(block_size)

	#print(keys)

	for key in transcript_id_info:

		chr = transcript_id_info[key][0]
		name = key
		strand = transcript_id_info[key][3]
		blocknum = len(dict_exon_size[key])
		blocksizes = dict_exon_size[key]
		exon_starts = dict_exon_start[key]
		if strand == '-':
			blocksizes = dict_exon_size[key][::-1]
			exon_starts = dict_exon_start[key][::-1]
		start = transcript_start_end[key][0]
		end = transcript_start_end[key][1]
		tstarts = [ x-start for x in exon_starts]

		group = transcript_id_info[key][1]
		blocktype = transcript_id_info[key][2]
		gene_id = transcript_id_info[key][4]
#               gene_type = transcript_id_info[key][5]
#		gene_status = transcript_id_info[key][6]
#		gene_name = transcript_id_info[key][7]
#		transcript_type = transcript_id_info[key][8]
#		transcript_status = transcript_id_info[key][9]
#		transcript_name = transcript_id_info[key][10]
#		level = transcript_id_info[key][11]


		#BED12 = [chr, str(start), str(end), name.strip('"'), gene_name.strip('"'), strand, transcript_name.strip('"'), transcript_type.strip('"'), str(level), str(blocknum), ",".join(map(str, blocksizes)) + "," , ",".join(map(str, tstarts)) + ",", group, gene_id.strip('"'), gene_status.strip('"'), transcript_status.strip('"'), gene_type.strip('"') ]

		BED12 = [chr, str(start), str(end), name.strip('"'), "0", strand, str(start), str(end), "0", str(blocknum), ",".join(map(str, blocksizes)) ]

		print "\t".join(BED12)

		#chr, start, end, name, gene_name, strand, transcript_name, transcript_type, level, blocknum, blocksizes, tstarts, group, gene_id, gene_status, transcript_status, gene_type




#		print key, dict_exon_size[key], dict_exon_start[key], dict_id_info[key]


	f.close()








if __name__ == '__main__':
	converter (sys.argv[1])

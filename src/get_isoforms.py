import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

annotated_ME  = set([])

chrom_sizes = {}



def make_chrom_sizes(fasta):

    f = open(fasta)

    for chrfa in SeqIO.parse(f, "fasta"):
        chrom_sizes[chrfa.id] = len(chrfa.seq)
    f.close()


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def main(annotation_bed12, annotation_gtf, out_filtered_ME, chrM):
	
    #chrM generates some indixing errors in mm10, that is why it should be excluded

    with open(annotation_bed12) as F:

        reader = csv.reader(F, delimiter="\t")

        for row in reader:
            chrom, start, end, transcript, score, strand, trickStart, trickEnd, score2, blocknumber, blocksizes, qstarts = row[:12]

            start = int(start)

            qstarts = list(map (int, qstarts.strip(",").split(",")))
            blocksizes = list(map(int, blocksizes.strip(",").split(",")))

            for q1, b in zip(qstarts, blocksizes):

                estart = start + q1
                eend = start + q1 + b
                elenght = eend - estart

                exon = "_".join(map(str, [chrom, strand, estart, eend]))

                annotated_ME.add(exon)

    ME_starts = dict()
    ME_ends =  dict()

    with open(out_filtered_ME) as F:

        reader = csv.DictReader(F, delimiter="\t")

        for row in reader:

            chrom = "_".join(row["ME"].split("_")[:-3])
            starnd, start, end = row["ME"].split("_")[-3:]

            ME_start = [chrom, start, strand]
            ME_end = [chrom, end, strand]


            for SJ in row["total_SJs"].split(","):

                try:

                    SJ_start, SJ_end = SJ.split(":")[1].split("-")

                except ValueError:

                    SJ_start, SJ_end = SJ.split(":")[1].split("+")


                if row["ME"] not in annotated_ME:
                    
                    if row["ME"]=="chr1_-_162075774_162075786":
                        print(row)

                    ME_starts[(chrom, int(SJ_start), strand)] = row["ME"]
                    ME_ends[(chrom, int(SJ_end), strand)] = row["ME"]

    transcript_to_gene = dict()
    gene_coordinates = dict()
    transcript_coordinates  = dict()

    for row in csv.reader(open(annotation_gtf), delimiter = '\t'):

        if row[0][0]!="#":

            chrom = row[0]
            strand =  row[6]
            feature = row[2]
            start = row[3]
            end = row[4]

            tags = row[8].strip(";").split("; ")

            tag_dict = dict()

            for t in tags:
                tag = t.strip(" ")
                if len(tag.split(' '))==2:
                    field, value = tag.split(' ')
                    value = value.strip('"')
                    tag_dict[field] = value


            gene_id = tag_dict["gene_id"]


            if "transcript_id" in tag_dict:

                transcript_id = tag_dict["transcript_id"]


                transcript_to_gene[transcript_id] = gene_id


            if feature=="gene":

                gene_coordinates[gene_id] = (chrom, start, end, strand)

            if feature=="transcript":

                transcript_coordinates[transcript_id] = (chrom, start, end, strand)

    non_ME_transcripts = defaultdict(list)
    ME_transcripts = defaultdict(list)

    gene_ME_intron = set([])


    with open(annotation_bed12) as F:

        reader = csv.reader(F, delimiter="\t")

        for row in reader:
            chrom, start, end, transcript_id, score, strand, trickStart, trickEnd, score2, blocknumber, blocksizes, qstarts = row[:12]

            qstarts = list(map (int, qstarts.strip(",").split(",")))
            blocksizes = list(map(int, blocksizes.strip(",").split(",")))

            start = int(start)
            end = int(end)

            blocknumber = int(blocknumber)


            if transcript_id in transcript_to_gene:


                gene_id = transcript_to_gene[transcript_id]
                g_chrom, g_start, g_end, g_strand = gene_coordinates[gene_id]
                t_chrom, t_start, t_end, t_strand = transcript_coordinates[transcript_id]


                past_intron = ""
                
                if transcript_id=="ENSMUST00000159763":
                    print("ENSMUST00000159763")


                    for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
                        estart = start + q1  + 1 #GTF is 1-based
                        eend = start + q1 + b
                        elenght = eend - estart

                        exon = [chrom, strand, estart, eend]

                        istart = start + q1 + b
                        iend = start + q2
                        ilen = iend - istart

                        intron  =  (chrom, strand, istart, iend)
                        non_ME_transcripts[transcript_id].append(exon)

                        if (tuple(exon) in set( map( tuple, ME_transcripts[transcript_id]) ) )==False:

                            ME_transcripts[transcript_id].append(exon)

                        if (chrom, eend, strand) in ME_starts:


                            ME = ME_starts[(chrom, eend, strand)]

                            ME = ME.split("_")
                            ME_chrom = "_".join(ME[:-3])
                            ME_strand, ME_start, ME_end  = ME[-3:]
                            ME_start = int(ME_start)
                            ME_end = int(ME_end)

                            ME_start += 1 ## GTF 1-based

                            ME = [ME_chrom, ME_strand, ME_start, ME_end ]
                            print(ME)


                            if ((gene_id, intron) in  gene_ME_intron) == False and int(t_start)<ME_start and int(t_end) > ME_end:

                                if (tuple(ME) in set( map( tuple, ME_transcripts[transcript_id]) ) )==False:

                                    ME_transcripts[transcript_id].append(ME)

                                    gene_ME_intron.add((gene_id, intron))


                        if (chrom, estart-1, strand) in ME_ends:


                            ME = ME_ends[(chrom, estart-1, strand)].split("_")
                            ME_chrom = "_".join(ME[:-3])
                            ME_strand, ME_start, ME_end  = ME[-3:]
                            ME_start = int(ME_start)
                            ME_end = int(ME_end)

                            ME_start += 1 ## GTF 1-based

                            ME = [ME_chrom, ME_strand, ME_start, ME_end ]
                            print(ME)

                            if ((gene_id, past_intron) in  gene_ME_intron) == False and int(t_start)<ME_start and int(t_end) > ME_end:

                                if (tuple(ME) in set( map( tuple, ME_transcripts[transcript_id]) ) )==False:

                                    ME_transcripts[transcript_id].insert(-2, ME)

                                    gene_ME_intron.add((gene_id, past_intron))


                        past_intron = intron


                    #last exon


                    q1, b = qstarts[-1], blocksizes[-1]
                    estart = start + q1  + 1 #GTF is 1-based
                    eend = start + q1 + b
                    elenght = eend - estart

                    exon = [chrom, strand, estart, eend]

                    non_ME_transcripts[transcript_id].append(exon)

                    ME_transcripts[transcript_id].append(exon)
				

            if transcript_id.split(".")[0] in transcript_to_gene:  ### zebrafish fix
		
		
                transcript_id = transcript_id.split(".")[0]

                gene_id = transcript_to_gene[transcript_id]
                g_chrom, g_start, g_end, g_strand = gene_coordinates[gene_id]
                t_chrom, t_start, t_end, t_strand = transcript_coordinates[transcript_id]


                past_intron = ""


                for q1, q2, b in zip(qstarts, qstarts[1:], blocksizes):
                    estart = start + q1  + 1 #GTF is 1-based
                    eend = start + q1 + b
                    elenght = eend - estart

                    exon = [chrom, strand, estart, eend]

                    istart = start + q1 + b
                    iend = start + q2
                    ilen = iend - istart

                    intron  =  (chrom, strand, istart, iend)
                    non_ME_transcripts[transcript_id].append(exon)
                    

                    

                    if (tuple(exon) in set( map( tuple, ME_transcripts[transcript_id]) ) )==False:

                        ME_transcripts[transcript_id].append(exon)

                    if (chrom, eend, strand) in ME_starts:

                        ME = ME_starts[(chrom, eend, strand)]

                        ME = ME.split("_")
                        ME_chrom = "_".join(ME[:-3])
                        ME_strand, ME_start, ME_end  = ME[-3:]
			
                        ME_start = int(ME_start)
                        ME_end = int(ME_end)

                        ME_start += 1 ## GTF 1-based			
			
                        ME = [ME_chrom, ME_strand, ME_start, ME_end ]


                        if ((gene_id, intron) in  gene_ME_intron) == False and int(t_start)<ME_start and int(t_end) > ME_end:

                            if (tuple(ME) in set( map( tuple, ME_transcripts[transcript_id]) ) )==False:

                                ME_transcripts[transcript_id].append(ME)

                                gene_ME_intron.add((gene_id, intron))


                    if (chrom, estart, strand) in ME_ends:

                        ME = ME_ends[(chrom, estart-1, strand)].split("_")
                        ME_chrom = "_".join(ME[:-3])
                        ME_strand, ME_start, ME_end  = ME[-3:]
                        ME_start = int(ME_start)
                        ME_end = int(ME_end)
			
                        ME_start += 1 ## GTF 1-based				
			
                        ME = [ME_chrom, ME_strand, ME_start, ME_end ]

                        if ((gene_id, past_intron) in  gene_ME_intron) == False and int(t_start)<ME_start and int(t_end) > ME_end:

                            if (tuple(ME) in set( map( tuple, ME_transcripts[transcript_id]) ) )==False:

                                ME_transcripts[transcript_id].insert(-2, ME)

                                gene_ME_intron.add((gene_id, past_intron))


                    past_intron = intron


                #last exon


                q1, b = qstarts[-1], blocksizes[-1]
                estart = start + q1  + 1 #GTF is 1-based
                eend = start + q1 + b
                elenght = eend - estart

                exon = [chrom, strand, estart, eend]

                non_ME_transcripts[transcript_id].append(exon)

                ME_transcripts[transcript_id].append(exon)	# fix zebrafish	
		
		
    printed_genes = set()

#     for transcript_id in non_ME_transcripts:

#         new_MEs = set([])
#         new_ME_starts = set([])
#         new_ME_ends = set([])

#         gene_id = transcript_to_gene[transcript_id]
#         g_chrom, g_start, g_end, g_strand = gene_coordinates[gene_id]
#         t_chrom, t_start, t_end, t_strand = transcript_coordinates[transcript_id]
        

        
#         if g_chrom in chrom_sizes:  #This is to avoid indexing problems with chromosomes that are in the anotation but not at the genome
            
#             if chrom_sizes[g_chrom] > int(g_end): # this should solve the mouse mithocodrial genes issues
            
#                 if  gene_id not in printed_genes:

#                     if chrM==False:
#                         if g_chrom!="chrM":
#                             print("\t".join(map(str, [ g_chrom, "MicroExonator", "gene", g_start, g_end, ".", g_strand, ".", "gene_id " +'"'+ gene_id +'"'+ ";" ])))
#                             printed_genes.add(gene_id)		
#                     else:
#                         print("\t".join(map(str, [ g_chrom, "MicroExonator", "gene", g_start, g_end, ".", g_strand, ".", "gene_id " +'"'+ gene_id +'"'+ ";" ])))
#                         printed_genes.add(gene_id)


#                 if chrM==False:
#                     if t_chrom!="chrM":
#                         print("\t".join(map(str, [  t_chrom, "MicroExonator", "transcript", t_start, t_end, ".", t_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id +'"'+ ";" ])))
#                 else:
#                     print("\t".join(map(str, [  t_chrom, "MicroExonator", "transcript", t_start, t_end, ".", t_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id +'"'+ ";" ])))


#                 for e in non_ME_transcripts[transcript_id]:

#                     e_chrom, e_strand, e_start, e_end = e


#                     if chrM==False:
#                         if t_chrom!="chrM":
#                             print("\t".join(map(str, [e_chrom, "MicroExonator", "exon", e_start, e_end, ".", e_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id +'"'+ ";"  ])))
#                     else:
#                         print("\t".join(map(str, [e_chrom, "MicroExonator", "exon", e_start, e_end, ".", e_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id +'"'+ ";"  ])))

#                 if len(ME_transcripts[transcript_id]) - len(non_ME_transcripts[transcript_id]) > 0:
#                 #if transcript_id in ME_transcripts:

#                     for e in ME_transcripts[transcript_id]:

#                         if tuple(e) not in set(map(tuple, non_ME_transcripts[transcript_id])):

#                             ME_chrom, ME_strand, ME_start, ME_end = e

#                             new_MEs.add(tuple(e))
#                             new_ME_starts.add(ME_start)
#                             new_ME_ends.add(ME_end)


#                     if len(new_ME_starts)!=0:
#                         transcript_id_ME =  "_".join( [  transcript_id, ".".join(map (str, new_ME_starts ))  ])
#                     else:
#                         transcript_id_ME = transcript_id
                    

#                     if chrM==False:
#                         if t_chrom!="chrM":
#                             print("\t".join(map(str, [  t_chrom, "MicroExonator", "transcript", t_start, t_end, ".", t_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id_ME +'"'+ ";" ])))
#                     else:
#                         print("\t".join(map(str, [  t_chrom, "MicroExonator", "transcript", t_start, t_end, ".", t_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id_ME +'"'+ ";" ])))



#                     for e in ME_transcripts[transcript_id]:

#                         e_chrom, e_strand, e_start, e_end = e


#                         if tuple(e) in new_MEs:

#                             if chrM==False:
#                                 if e_chrom!="chrM":
#                                     print("\t".join(map(str, [e_chrom, "MicroExonator", "exon", e_start, e_end, ".", e_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id_ME +'"'+ ";"  ])))
#                             else:
#                                 print("\t".join(map(str, [e_chrom, "MicroExonator", "exon", e_start, e_end, ".", e_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id_ME +'"'+ ";"  ])))

#                         elif e_start not in new_ME_starts and e_end not in new_ME_ends:

#                             if chrM==False:
#                                 if e_chrom!="chrM":
#                                     print("\t".join(map(str, [e_chrom, "MicroExonator", "exon", e_start, e_end, ".", e_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id_ME +'"'+ ";"  ])))
#                             else:
#                                 print("\t".join(map(str, [e_chrom, "MicroExonator", "exon", e_start, e_end, ".", e_strand, ".", "gene_id " +'"'+ gene_id +'"'+ "; " + "transcript_id " +'"'+ transcript_id_ME +'"'+ ";"  ])))



if __name__ == '__main__':
	make_chrom_sizes(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4], str2bool(sys.argv[5]))

import sys
import csv
import argparse




samples = []

def urls_to_download(urls, split=False):

    splits = []


    if split:
        for row in csv.reader(open(split)):
            splits.append(row[0])

    for row in csv.reader(open(urls)):

        D_folder = "download/"
        O_folder = "data/fastq/"
        url = row[0]
        FILE =  url.split("/")[-1]
        basename = FILE.split(".")[0]
        ext = ".".join(FILE.split(".")[1:])

        samples.append(FILE)
        download_script = open(D_folder + basename + ".download.sh", "w")
        wget = ["wget -r", url, "-O", D_folder + FILE]
        to_fastq = []
        split_pairs = []
        move = []
        rm = []


        if  ext=="sra":

            to_fastq = ["fastq-dump -O download", D_folder + FILE]


            if FILE in splits:
                split_pairs = ["python2 src/split_paired_end.py", D_folder +  basename + ".fastq > ", D_folder + basename + ".fastq.split"]
                move = ["mv", D_folder + basename + ".fastq.split", O_folder + basename + ".fastq"]
                rm = ["rm", D_folder + FILE, D_folder + basename + ".fastq"]

            else:
                move = ["mv", D_folder + basename + ".fastq", O_folder + basename + ".fastq"]
                rm = ["rm", D_folder + FILE]



        if  ext=="fastq.gz":

            to_fastq = ["gzip -d", D_folder + FILE]


            if FILE in splits:
                split_pairs = ["python2 src/split_paired_end.py", D_folder +  basename + ".fastq > ", D_folder + basename + ".fastq.split"]
                move = ["mv", D_folder + basename + ".fastq.split", O_folder + basename + ".fastq"]
                rm = ["rm", D_folder + basename + ".fastq"]

            else:
                move = ["mv", D_folder + basename + ".fastq", O_folder + basename + ".fastq"]



        download_script.write("#!/bin/bash" + "\n")
        download_script.write("#" + ext + " detected" + "\n")
        download_script.write(" ".join(wget) + "\n")
        download_script.write(" ".join(to_fastq)+ "\n")
        download_script.write(" ".join(split_pairs) + "\n")
        download_script.write(" ".join(move)+ "\n")
        download_script.write(" ".join(rm)+ "\n")

        download_script.close()

def accession_to_download(accession, split=False):

    splits = []


    if split:
        for row in csv.reader(open(split)):
            splits.append(row[0].split(".")[0]) #Only get the basename

    for row in csv.reader(open(accession)):

        D_folder = "download/"
        O_folder = "data/fastq/"

        SRA = row[0]

        samples.append(SRA)
        download_script = open(D_folder + SRA + ".download.sh", "w")

        to_fastq = []
        merge = []
        move = ["mv", D_folder + SRA + ".fastq", O_folder + SRA + ".fastq"]
        rm = []

        if SRA in splits:

            to_fastq = ["fastq-dump --split-files -O", D_folder, "--accession", SRA]
            merge = ["python src/merge_pairs.py", D_folder + SRA + "_1.fastq", D_folder + SRA + "_2.fastq", D_folder + SRA + ".fastq"]
            rm = ["rm", D_folder + SRA + "_1.fastq", D_folder + SRA + "_2.fastq"]

        else:

            to_fastq = ["fastq-dump -O", D_folder, "--accession", SRA]



        download_script.write("#!/bin/bash" + "\n")
        download_script.write(" ".join(to_fastq)+ "\n")
        download_script.write(" ".join(merge) + "\n")
        download_script.write(" ".join(move)+ "\n")
        download_script.write(" ".join(rm)+ "\n")

        download_script.close()




def main():

    with open("config.yaml", "w") as out:

        parser = argparse.ArgumentParser()
        parser.add_argument("-g", "--genome",  help="Multi-fastq that contains the chromosomes", required=True)
        parser.add_argument("-gab", "--bed12",  help="Gene anotation formated as bed12", required=True)
        parser.add_argument("-gaf", "--transcripts",  help="Multi-fasta containing the transcript sequences", required=True)
        parser.add_argument("-U2_5", "--GT_AG_U2_5",  help="GT-AG U5 5' PWM", required=True)
        parser.add_argument("-U2_3", "--GT_AG_U2_3",  help="GT-AG U5 3' PWM", required=True)
        parser.add_argument("-p", "--phylop",  help="pylop", required=True)
        parser.add_argument("-l", "--ME_len", help="Maximun micro-exon lenght", required=True)
        parser.add_argument("-wd", "--working_directory", help="path to the working directory folder", required=True)
        #parser.add_argument("-p2", "--phylop2",  help="pylop2", required=True) #

        parser.add_argument("-db", "--data_base", help="bed12 or bed6 file with microexon anotation", required=False)
        parser.add_argument("-i", "--input_dir",  help="Path to folder containing the RNA-samples", required=False)
        parser.add_argument('samples', nargs='*', help="RNA-seq samples formated as .fastq, .SRA, or .fastq.gz. They need to be named with the right extension")



        parser.add_argument("-u", "--url",  help="list of URLs that link to fastq.gz or SRA", required=False)
        parser.add_argument("-s", "--split_paired_end",  help="list of fastq.gz or SRA files which correspond to collapsed paired-end files", required=False)


        parser.add_argument("-AI", "--accession", help="SRA accession list", required=False)


        args = parser.parse_args()

        out.write("Genome_fasta" + " : " +  args.genome + "\n")
        out.write("Gene_anontation_bed12" + " : " +  args.bed12 + "\n")
        out.write("Gene_anontation_fasta" + " : " +  args.transcripts + "\n")

        out.write("GT_AG_U2_5" + " : " +  args.GT_AG_U2_5 + "\n")
        out.write("GT_AG_U2_3" + " : " +  args.GT_AG_U2_3 + "\n")
        out.write("vertebrates_phylop" + " : " +  args.phylop + "\n")
        out.write("ME_len" + " : " +  args.ME_len + "\n")
        out.write("working_directory" + " : " +  args.working_directory.rstrip("/") + "/" + "\n")


        if args.data_base is not None:
            out.write("ME_DB" + " : " + args.data_base + "\n")


        if args.url is None and args.accession is None:

            out.write("input_dir" + " : " + args.input_dir.rstrip("/") + "\n")
            out.write("samples" + " : " +  str(args.samples) + "\n")

        if args.url is not None:


            if args.split_paired_end is None:

                urls_to_download(args.url)

            else:

                urls_to_download(args.url, args.split_paired_end)


        if args.accession is not None:


            if args.split_paired_end is None:

                accession_to_download(args.accession)

            else:

                accession_to_download(args.accession, args.split_paired_end)




        if args.url is not None or args.accession is not None:

            out.write("input_dir" + " : download" + "\n")
            out.write("samples" + " : " +  str(samples) + "\n")


main()


#python config.py -g /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/mm10.fa -gab /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.v19.chr_patch_hapl_scaff.annotation.bed12 -gaf /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.transcripts.fa -U2_5 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix -U2_3 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix -p /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw -i ../data ../data/*.fastq.gz



#python config.py -g /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/mm10.fa -gab /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.bed12 -gaf /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.transcripts.fa -U2_5 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix -U2_3 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix -p1 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw -p2 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60wayPlacental.bw -i ../../Single_cell/DRG_neurons/ ../../Single_cell/DRG_neurons/SRR1660*

#python config.py -g /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/mm10.fa -gab /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.bed12 -gaf /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Gene_annotation/gencode.vM11.transcripts.fa -U2_5 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix -U2_3 /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix -p /lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw -u url/URL_list/test.url.txt -AI url/URL_list/test.sra.txt -s url/URL_list/split.txt -l 30 -wd /lustre/scratch117/cellgen/team218/gp7/Micro-exons/Software/Micro-Exonator

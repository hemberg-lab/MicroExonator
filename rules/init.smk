
import glob, os
import random
import csv
from collections import defaultdict




try:
    os.mkdir("FASTQ")
except FileExistsError:
    pass

try:
    os.mkdir("download")
except FileExistsError:
    pass
    
    
try:
    os.mkdir("logs")
except FileExistsError:
    pass


#DATA = set([]) #Defined now at the beining of MicroExonator.smk

if os.path.isfile('./NCBI_accession_list.txt'):


    with open("NCBI_accession_list.txt") as file :

        reader = csv.reader(file, delimiter="\t")

        for row in reader:

            RUN = row[0]
            DATA.add(RUN)

            file_name = "download/" + RUN + ".download.sh"
            command = "fastq-dump.2.11.0 --split-files -O FASTQ --gzip --defline-qual '+'"

            if len(glob.glob(file_name))==0: #Check if the file is there, as if this file is overwriten everything will start from scratch

                download_file =  open(file_name, "w")

                download_file.write("#!/bin/bash" + "\n")
                download_file.write('srr="' + RUN + '"' + "\n" )
                download_file.write(command + " " + RUN + "\n")
                download_file.write( "numLines=$(fastq-dump.2.11.0 -X 1 -Z --split-spot $srr | wc -l)" + "\n")
                download_file.write( "if [ $numLines -eq 8 ]; then cat FASTQ/${srr}_1.fastq.gz FASTQ/${srr}_2.fastq.gz > FASTQ/$srr.fastq.gz && rm FASTQ/${srr}_1.fastq.gz FASTQ/${srr}_2.fastq.gz; fi"  + "\n")
                download_file.write( "if [ -f FASTQ/${srr}_1.fastq.gz ]; then mv FASTQ/${srr}_1.fastq.gz FASTQ/${srr}.fastq.gz ; elif [ -f FASTQ/${srr}_2.fastq.gz ]; then mv FASTQ/${srr}_2.fastq.gz FASTQ/${srr}.fastq.gz; fi"  + "\n")

#                download_file.write("gzip FASTQ/${srr}.fastq" + "\n")
                

if os.path.isfile("./local_samples.tsv"):


    with open("./local_samples.tsv") as file :

        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:

            DATA.add(row["sample"])

            file_name = "download/" + row["sample"]  + ".download.sh"


            if len(glob.glob(file_name))==0: #Check if the file is there, as if this file is overwriten everything will start from scratch

                if row["path"].split(".")[-1] == "gz":

                    download_file =  open(file_name, "w")

                    download_file.write("#!/bin/bash" + "\n")
                    download_file.write("cp " + row["path"] +  " FASTQ/" + row["sample"]  + ".fastq.gz" + "\n")
                
                
                elif row["path"].split(".")[-1] == "bz2":
                    download_file =  open(file_name, "w")

                    download_file.write("#!/bin/bash" + "\n")
                    download_file.write("bzcat " + row["path"] +  " | gzip > FASTQ/" + row["sample"]  + ".fastq.gz" + "\n")

                elif row["path"].split(".")[-1] == "fastq":
                
                    download_file =  open(file_name, "w")

                    download_file.write("#!/bin/bash" + "\n")
                    download_file.write("cat " + row["path"] +  " | gzip > FASTQ/" + row["sample"]  + ".fastq.gz" + "\n")

                else:
                    print("Error: only fastq, fastq.gz or fastq.bz2 are suported!")

if os.path.isfile("./sample_url.tsv"):
    
    with open("./sample_url.tsv") as file :

        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:

            DATA.add(row["sample"])

            file_name = "download/" + row["sample"]  + ".download.sh"


            if len(glob.glob(file_name))==0: #Check if the file is there, as if this file is overwriten everything will start from scratch

                download_file =  open(file_name, "w")

                download_file.write("#!/bin/bash" + "\n")
                download_file.write("wget -r " + row["url"] +  " -O FASTQ/" + row["sample"]  + ".fastq.gz" + "\n")
                #download_file.write("gzip -d FASTQ/" + row["sample"]  + ".fastq.gz" + "\n")




if os.path.isfile("./local_samples_cram.tsv"):

    


    with open("./local_samples_cram.tsv") as file :

        reader = csv.DictReader(file, delimiter="\t")
        
        files_cramtofastq = defaultdict(list)
        files_subcat  = defaultdict(list)

        for row in reader:
            
            basename = row["sample"]
            
            if "multiplex" in row:
                
                basename = "_".join([row["sample"],  row["multiplex"]])
            

            DATA.add(row["sample"])

            file_name = "download/" + row["sample"]  + ".download.sh"



            if len(glob.glob(file_name))==0: #Check if the file is there, as if this file is overwriten everything will start from scratch

                

                if "multiplex" in row: 
                
                    out_name =  row["sample"] 
                
                    #if row["multiplex"] == config["multiplex_number"]:
                    
                    #download_file.write("cat " + subcat + "> FASTQ/" + row["sample"]  + ".fastq " + "\n")

                    cramtofastq = "samtools view -b " + row["path"] +  " | samtools bam2fq -1 FASTQ/" + basename  + ".rd1.fastq -2 FASTQ/" +  basename + ".rd2.fastq -0 /dev/null -s /dev/null -n -F 0x900 - &&" + "\n"
                    subcat = "FASTQ/" + basename  + ".rd1.fastq FASTQ/" +  basename + ".rd2.fastq "

                    files_cramtofastq[out_name].append(cramtofastq)
                    files_subcat[out_name].append(subcat)
                else:
                
                    download_file =  open(file_name, "w")

                    download_file.write("#!/bin/bash" + "\n")
                    download_file.write("samtools view -b " + row["path"] +  " | samtools bam2fq -1 FASTQ/" + basename  + ".rd1.fastq -2 FASTQ/" +  basename + ".rd2.fastq -0 /dev/null -s /dev/null -n -F 0x900 - &&" + "\n")
                
                    download_file.write("cat  FASTQ/" + basename  + ".rd1.fastq FASTQ/" +  basename + ".rd2.fastq > " + " FASTQ/" + row["sample"]  + ".fastq " + "\n") 

        
        for out_name in files_cramtofastq:

            file_name = "download/" + out_name  + ".download.sh"


            download_file =  open(file_name, "w")
            download_file.write("#!/bin/bash" + "\n")

            cramtofastqs = files_cramtofastq[out_name]
            subcats = files_subcat[out_name]

            for cramtofastq in cramtofastqs:

                download_file.write(cramtofastq)

            cat = "cat " + "".join(subcats) +  " > FASTQ/" + out_name + ".fastq &&" + "\n"
            rm  = "rm " + "".join(subcats)  + "\n"

            download_file.write(cat)
            download_file.write(rm)


#print(DATA)

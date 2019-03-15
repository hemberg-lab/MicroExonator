import csv
import sys
from collections import defaultdict




def main(gene_model_bed12, out_filtered_ME_cov, out_filtered_ME, out_low_scored_ME):

    estart_exons = defaultdict(set)
    eend_exons = defaultdict(set)

    total_ME_up = defaultdict(int)
    total_ME_down = defaultdict(int)

    high_confident_ME = []

    with open(gene_model_bed12) as bedfile, \
        open(out_filtered_ME_cov) as ME_out_cov, \
        open(out_filtered_ME) as ME_out1, \
        open(out_low_scored_ME) as ME_low1, \
        open(out_filtered_ME) as ME_out, \
        open(out_low_scored_ME) as ME_low:

        reader = csv.reader(bedfile, delimiter="\t")

        for row in reader:

            csv.field_size_limit(1000000000)

            qstarts = list(map (int, row[11].strip(",").split(",")))[1:-1]
            blocksizes = list(map(int, row[10].strip(",").split(",")))[1:-1]

            start = int(row[1])
            strand = row[5]
            bn = int(row[9])
            chrom = row[0]


            for q1, b in zip(qstarts, blocksizes):
                estart = start + q1
                eend = start + q1 + b
                elenght = eend - estart
                exon = (chrom, strand, estart, eend)

                estart_exons[(chrom, strand, estart)].add(exon)
                eend_exons[(chrom, strand, eend)].add(exon)


        #### Counting exons that have the same start/end


        reader = csv.DictReader(ME_out1, delimiter="\t")


        for row in reader:

            chrom = "_".join(row["ME"].split("_")[:-3]) 
            strand, estart, eend = row["ME"].split("_")[-3:]
            exon = (chrom, strand, estart, eend)

            estart_exons[(chrom, strand, estart)].add(exon)
            eend_exons[(chrom, strand, eend)].add(exon)



        reader = csv.DictReader(ME_low1, delimiter="\t")


        for row in reader:

            if row["ME_type"]=="RESCUED":
                
                chrom = "_".join(row["ME"].split("_")[:-3]) 
                strand, estart, eend = row["ME"].split("_")[-3:]
                exon = (chrom, strand, estart, eend)

                estart_exons[(chrom, strand, estart)].add(exon)
                eend_exons[(chrom, strand, eend)].add(exon)


        ## Summing exon coverage


        reader = csv.DictReader(ME_out_cov, delimiter="\t")

        for row in reader:

            chrom = "_".join(row["ME"].split("_")[:-3]) 
            strand, estart, eend = row["ME"].split("_")[-3:]
            exon = (chrom, strand, estart, eend)

            sum_ME_SJ_coverage_up = int(row["sum_ME_SJ_coverage_up"])
            sum_ME_SJ_coverage_down = int(row["sum_ME_SJ_coverage_down"])


            total_ME_up[row["ME"]] += sum_ME_SJ_coverage_up
            total_ME_down[row["ME"]] += sum_ME_SJ_coverage_down




        print("ME", "transcript", "sum_total_coverage", "total_SJs", "total_coverages", "len_micro_exon_seq_found", "micro_exon_seq_found", "total_number_of_micro_exons_matches", "U2_scores",  "mean_conservations_vertebrates", "P_MEs", "total_ME",   "ME_P_value", "ME_type", sep="\t")

        reader = csv.DictReader(ME_out, delimiter="\t")

        for row in reader:

            chrom = "_".join(row["ME"].split("_")[:-3]) 
            strand, estart, eend = row["ME"].split("_")[-3:]
            exon = (chrom, strand, estart, eend)

            sum_ME_SJ_coverage_up = total_ME_up[row["ME"]]
            sum_ME_SJ_coverage_down =  total_ME_down[row["ME"]]

            abs_up_down_diff = "NA"


            # if row["ME"]=="chr11_-_41913982_41913993":
            #
            #     print(len(estart_exons[(chrom, strand, estart)]), len(eend_exons[(chrom, strand, eend)]),  abs(sum_ME_SJ_coverage_up-sum_ME_SJ_coverage_down)/(sum_ME_SJ_coverage_up+sum_ME_SJ_coverage_down) )
            #     print( row["ME"], row["transcript"], row["sum_total_coverage"], row["total_SJs"], row["total_coverages"], row["len_micro_exon_seq_found"], row["micro_exon_seq_found"], row["total_number_of_micro_exons_matches"], row["U2_scores"], row["mean_conservations_vertebrates"], row["P_MEs"], row["total_ME"], row["ME_P_value"], row["ME_type"], sep="\t")
            #



            if sum_ME_SJ_coverage_up+sum_ME_SJ_coverage_down>0:


                abs_up_down_diff = abs(sum_ME_SJ_coverage_up-sum_ME_SJ_coverage_down)/(sum_ME_SJ_coverage_up+sum_ME_SJ_coverage_down)


                if (len(estart_exons[(chrom, strand, estart)]) + len(eend_exons[(chrom, strand, eend)]) > 2) and abs_up_down_diff > 0.95:

                    #ME_black_list.write(row["ME"]+"\n")
                    pass

                else:
                    print(row["ME"], row["transcript"], row["sum_total_coverage"], row["total_SJs"], row["total_coverages"], row["len_micro_exon_seq_found"], row["micro_exon_seq_found"], row["total_number_of_micro_exons_matches"], row["U2_scores"], row["mean_conservations_vertebrates"], row["P_MEs"], row["total_ME"], row["ME_P_value"], row["ME_type"], sep="\t")



        reader = csv.DictReader(ME_low, delimiter="\t")

        for row in reader:

            chrom = "_".join(row["ME"].split("_")[:-3]) 
            strand, estart, eend = row["ME"].split("_")[-3:]
            exon = (chrom, strand, estart, eend)

            sum_ME_SJ_coverage_up = total_ME_up[row["ME"]]
            sum_ME_SJ_coverage_down =  total_ME_down[row["ME"]]

            abs_up_down_diff = "NA"


            if sum_ME_SJ_coverage_up+sum_ME_SJ_coverage_down>0 and row["ME_type"]=="RESCUED":

                abs_up_down_diff = abs(sum_ME_SJ_coverage_up-sum_ME_SJ_coverage_down)/(sum_ME_SJ_coverage_up+sum_ME_SJ_coverage_down)

                if (len(estart_exons[(chrom, strand, estart)]) + len(eend_exons[(chrom, strand, eend)]) > 2) and abs_up_down_diff > 0.95:

                    #ME_black_list.write(row["ME"]+"\n")
                    pass

                else:
                    print(row["ME"], row["transcript"], row["sum_total_coverage"], row["total_SJs"], row["total_coverages"], row["len_micro_exon_seq_found"], row["micro_exon_seq_found"], row["total_number_of_micro_exons_matches"], row["U2_scores"], row["mean_conservations_vertebrates"], row["P_MEs"], row["total_ME"], row["ME_P_value"], row["ME_type"], sep="\t")



main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

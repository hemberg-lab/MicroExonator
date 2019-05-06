import sys
import csv

def binP(N, p, x1, x2):
    p = float(p)
    q = p/(1-p)
    k = 0.0
    v = 1.0
    s = 0.0
    tot = 0.0

    while(k<=N):
            tot += v
            if(k >= x1 and k <= x2):
                    s += v
            if(tot > 10**30):
                    s = s/10**30
                    tot = tot/10**30
                    v = v/10**30
            k += 1
            v = v*q*(N+1-k)/k
    return s/tot

def calcBin(vx, vN, vCL = 95):
    '''
    Calculate the exact confidence interval for a binomial proportion

    Usage:
    >>> calcBin(13,100)
    (0.07107391357421874, 0.21204372406005856)
    >>> calcBin(4,7)
    (0.18405151367187494, 0.9010086059570312)
    '''
    vx = float(vx)
    vN = float(vN)
    #Set the confidence bounds
    vTU = (100 - float(vCL))/2
    vTL = vTU

    vP = vx/vN
    if(vx==0):
            dl = 0.0
    else:
            v = vP/2
            vsL = 0
            vsH = vP
            p = vTL/100

            while((vsH-vsL) > 10**-5):
                    if(binP(vN, v, vx, vN) > p):
                            vsH = v
                            v = (vsL+v)/2
                    else:
                            vsL = v
                            v = (v+vsH)/2
            dl = v

    if(vx==vN):
            ul = 1.0
    else:
            v = (1+vP)/2
            vsL =vP
            vsH = 1
            p = vTU/100
            while((vsH-vsL) > 10**-5):
                    if(binP(vN, v, 0, vx) < p):
                            vsH = v
                            v = (vsL+v)/2
                    else:
                            vsL = v
                            v = (v+vsH)/2
            ul = v
    return (dl, ul)


def main(ME_SJ_coverage, min_sum_PSI):  #, path):

    with open(ME_SJ_coverage) as F:

        reader = csv.reader(F, delimiter="\t")

        #print( "\t".join( [ "ME", "Coord", "PSI", "CI_Lo", "CI_Hi", "is_alternative_5", "alternatives_5",  "is_alternative_3", "alternatives_3" ]) )

        for row in reader:

            FILE, ME, total_SJs, ME_SJ_coverages, sum_ME_coverage, sum_ME_SJ_coverage_up_down_uniq, sum_ME_SJ_coverage_up, sum_ME_SJ_coverage_down, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3,  total_cov_alternatives_3 = row


            SUM_PSI = float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5)
            if SUM_PSI>=min_sum_PSI:

                PSI= float(sum_ME_coverage)/(float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5))

                CI_Lo, CI_Hi = calcBin(float(sum_ME_coverage),  SUM_PSI)

            else:

                PSI = "NA"
                CI_Lo, CI_Hi = ["NA", "NA"]


            ME_strand, ME_start, ME_end = ME.split("_")[-3:]  #to avoid errors with chromosomes like chr1_GL456210_random
            ME_chrom = "_".join(ME.split("_")[:-3])


            Coord = ME_chrom + ":" + str(int(ME_start)+1) + "-" + ME_end

	
#             if path[-1]!="/":
#                 path += "/"		
#             with open(path + FILE + ".sam.pre_processed.filter1.ME_SJ_coverage.PSI", "a") as out:	
	
            if is_alternative_5=="True":


                for alt5, alt5_cov in zip(alternatives_5.split(","), cov_alternatives_5.split(",")):

                    alt5_crom, alt5_loci = alt5.split(":")

                    if "+" in alt5_loci:


                        alt5_start, alt5_end = alt5_loci.split("+")
                        alt5_strand = "+"

                        alt5_Coord_start = str(int(alt5_start) + 1)
                        alt5_Coord_end = str(int(ME_start) )

                        alt5_Coord = alt5_crom + ":" + alt5_Coord_start + "-" + alt5_Coord_end


                    elif "-" in alt5_loci:

                        alt5_start, alt5_end = alt5_loci.split("-")
                        alt5_strand = "-"

                        alt5_Coord_start = str(int(ME_end) + 1)
                        alt5_Coord_end = str(int(alt5_end))

                        alt5_Coord = alt5_crom + ":" + alt5_Coord_start + "-" + alt5_Coord_end



                    alt5_SUM_PSI = float(alt5_cov) + (float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+(float(total_cov_alternatives_5)- float(alt5_cov) ) + float(sum_ME_coverage)  )

                    if alt5_SUM_PSI>=min_sum_PSI:

                        alt5_PSI= float(alt5_cov)/(float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+ float(total_cov_alternatives_5) + float(sum_ME_coverage)  )

                        alt5_CI_Lo, alt5_CI_Hi = calcBin(float(alt5_cov),  SUM_PSI)

                    else:

                        alt5_PSI = "NA"
                        alt5_CI_Lo, alt5_CI_Hi = ["NA", "NA"]


                    #out.write( "\t".join( map(str, [ alt5, alt5_Coord, alt5_PSI, alt5_CI_Lo, alt5_CI_Hi, "alt5" ])) )
                    print("\t".join( map(str, [ alt5, alt5_Coord, alt5_PSI, alt5_CI_Lo, alt5_CI_Hi, "alt5" ])))


            if is_alternative_3=="True":

                for alt3, alt3_cov in zip(alternatives_3.split(","), cov_alternatives_3.split(",")):

                    alt3_crom, alt3_loci = alt3.split(":")

                    if "-" in alt3_loci:


                        alt3_start, alt3_end = alt3_loci.split("-")
                        alt3_strand = "-"

                        alt3_Coord_start = str(int(alt3_start) + 1)
                        alt3_Coord_end = str(int(ME_start))

                        alt3_Coord = alt3_crom + ":" + alt3_Coord_start + "-" + alt3_Coord_end


                    elif "+" in alt3_loci:

                        alt3_start, alt3_end = alt3_loci.split("+")
                        alt3_strand = "+"

                        alt3_Coord_start = str(int(ME_end) + 1)
                        alt3_Coord_end = str(int(alt3_end))

                        alt3_Coord = alt3_crom + ":" + alt3_Coord_start + "-" + alt3_Coord_end



                    alt3_SUM_PSI = float(alt3_cov) + (float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+(float(total_cov_alternatives_5)- float(alt3_cov) ) + float(sum_ME_coverage)  )

                    if alt3_SUM_PSI>=min_sum_PSI:

                        alt3_PSI= float(alt3_cov)/(float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+ float(total_cov_alternatives_5) + float(sum_ME_coverage)  )

                        alt3_CI_Lo, alt3_CI_Hi = calcBin(float(alt3_cov),  SUM_PSI)

                    else:

                        alt3_PSI = "NA"
                        alt3_CI_Lo, alt3_CI_Hi = ["NA", "NA"]


                    #out.write( "\t".join( map(str, [ alt3, alt3_Coord, alt3_PSI, alt3_CI_Lo, alt3_CI_Hi, "alt3" ])) + "\n" )
                    print( "\t".join( map(str, [ alt3, alt3_Coord, alt3_PSI, alt3_CI_Lo, alt3_CI_Hi, "alt3" ])))




            #out.write( "\t".join( map(str, [ ME, Coord, PSI, CI_Lo, CI_Hi, "ME" ])) + "\n" )
            print( "\t".join( map(str, [ ME, Coord, PSI, CI_Lo, CI_Hi, "ME" ])))


if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2])) #, sys.argv[3]  )

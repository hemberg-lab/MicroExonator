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


def main(total_cov, min_sum_PSI, paired):
	
  print("File", "Unique_ME_reads", "ME_coords", "SJ_coords", "ME_coverages", "SJ_coverages", "PSI", "CI_Lo", "CI_Hi", "Alt5", "Alt3", "Alt5_coverages", "Alt3_coverages", sep="\t")

  paired_files = set([])
  pair12 = dict()
	
  if paired!="F":
		
    for row in csv.reader(open(paired), delimiter="\t" ):
			
      pair1, pair2 = row
			
      paired_files.add(pair1)
      paired_files.add(pair2)
      pair12[pair1] = pair2
		
  with open(total_cov) as file :
    
    paired_info = dict()
        
    #reader = csv.DictReader(file, delimiter="\t")  #This code used to work for input this header, and now is being re-adapted
    reader = csv.reader(file, delimiter="\t")
    
    for line in reader:
      
      #FILE_NAME, ME, total_SJs, ME_SJ_coverages, sum_ME_coverage, sum_ME_SJ_coverage_up_down_uniq, sum_ME_SJ_coverage_up, sum_ME_SJ_coverage_down, SJ_coverages, sum_SJ_coverage, is_alternative_5, is_alternative_3, alternatives_5, cov_alternatives_5, total_cov_alternatives_5, alternatives_3, cov_alternatives_3, total_cov_alternatives_3 = row
      header = ['FILE_NAME', 'ME', 'total_SJs', 'ME_SJ_coverages', 'sum_ME_coverage', 'sum_ME_SJ_coverage_up_down_uniq', 'sum_ME_SJ_coverage_up', 'sum_ME_SJ_coverage_down', 'SJ_coverages', 'sum_SJ_coverage', 'is_alternative_5', 'is_alternative_3', 'alternatives_5', 'cov_alternatives_5', 'total_cov_alternatives_5', 'alternatives_3', 'cov_alternatives_3', 'total_cov_alternatives_3'] 
      row = dict()
      for col, col_name in zip(line, header):
          row[col_name] = col
      

      sum_ME_coverage = row["sum_ME_coverage"]
      sum_SJ_coverage = row["sum_SJ_coverage"]
      total_cov_alternatives_3 = row["total_cov_alternatives_3"]
      total_cov_alternatives_5 = row["total_cov_alternatives_5"]
	
      Unique_ME_reads = row["sum_ME_SJ_coverage_up_down_uniq"]
      
      SUM_PSI = float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5)
      if SUM_PSI>=min_sum_PSI:

          PSI= float(sum_ME_coverage)/(float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5))

          CI_Lo, CI_Hi = calcBin(float(sum_ME_coverage),  SUM_PSI)

      else:

          PSI = "NA"
          CI_Lo, CI_Hi = ["NA", "NA"]
          
      if row["FILE_NAME"] in paired_files:
        
        info = [Unique_ME_reads, sum_ME_coverage, sum_SJ_coverage, total_cov_alternatives_3, total_cov_alternatives_5, row["FILE_NAME"], row["ME"], row["total_SJs"], row["ME_SJ_coverages"], row["SJ_coverages"], row["alternatives_5"], row["alternatives_3"], row["cov_alternatives_5"], row["cov_alternatives_3"]]
        paired_info[(row["FILE_NAME"], row["ME"])] = info

      else:
      
        print(row["FILE_NAME"], row["sum_ME_SJ_coverage_up_down_uniq"], row["ME"], row["total_SJs"], row["ME_SJ_coverages"], row["SJ_coverages"], PSI, CI_Lo, CI_Hi, row["alternatives_5"], row["alternatives_3"], row["cov_alternatives_5"], row["cov_alternatives_3"], sep='\t')
  
    for p, p1info in paired_info.items():

      FILE_NAME_1, ME_1 = p

      if FILE_NAME_1 in pair12:
        FILE_NAME_2 = pair12[FILE_NAME_1]
        p2info = paired_info[(FILE_NAME_2, ME_1)]
        
        Unique_ME_reads_1, sum_ME_coverage_1, sum_SJ_coverage_1, total_cov_alternatives_3_1, total_cov_alternatives_5_1, FILE_NAME, ME, total_SJs, ME_SJ_coverages, SJ_coverages, alternatives_5, alternatives_3, cov_alternatives_5, cov_alternatives_3 = p1info
        Unique_ME_reads_2, sum_ME_coverage_2, sum_SJ_coverage_2, total_cov_alternatives_3_2, total_cov_alternatives_5_2 = p2info[:5]
        
        Unique_ME_reads = int(Unique_ME_reads_1) + int(Unique_ME_reads_2)
        sum_ME_coverage = int(sum_ME_coverage_1) + int(sum_ME_coverage_2)
        sum_SJ_coverage = int(sum_SJ_coverage_1) + int(sum_SJ_coverage_2)
        total_cov_alternatives_3 = int(total_cov_alternatives_3_1) + int(total_cov_alternatives_3_2)
        total_cov_alternatives_5 = int(total_cov_alternatives_5_1) + int(total_cov_alternatives_5_2)
  
        SUM_PSI = float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5)
        if SUM_PSI>=min_sum_PSI:

            PSI= float(sum_ME_coverage)/(float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5))

            CI_Lo, CI_Hi = calcBin(float(sum_ME_coverage),  SUM_PSI)

        else:

            PSI = "NA"
            CI_Lo, CI_Hi = ["NA", "NA"]
            
        print(FILE_NAME, Unique_ME_reads, ME, total_SJs, ME_SJ_coverages, SJ_coverages, PSI, CI_Lo, CI_Hi, alternatives_5, alternatives_3, cov_alternatives_5, cov_alternatives_3, sep="\t")

  
if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2]), paired=sys.argv[3] )

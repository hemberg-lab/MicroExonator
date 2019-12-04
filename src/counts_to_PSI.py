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


def main(total_cov):
  
  with open(total_cov) as file :
    
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
      
      #FILE_NAME       ME      total_SJs       ME_SJ_coverages sum_ME_coverage sum_ME_SJ_coverage_up_down_uniq sum_ME_SJ_coverage_up   sum_ME_SJ_coverage_down SJ_coverages    sum_SJ_coverage is_alternative_5        is_alternative_3        alternatives_5  cov_alternatives_5      total_cov_alternatives_5        alternatives_3  cov_alternatives_3      total_cov_alternatives_3
      
      sum_ME_coverage = row["sum_ME_coverage"]
      sum_SJ_coverage = row["sum_SJ_coverage"]
      total_cov_alternatives_3 = row["total_cov_alternatives_3"]
      total_cov_alternatives_5 = row["total_cov_alternatives_5"]
      
      SUM_PSI = float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5)
      if SUM_PSI>=min_sum_PSI:

          PSI= float(sum_ME_coverage)/(float(sum_ME_coverage)+float(sum_SJ_coverage)+float(total_cov_alternatives_3)+float(total_cov_alternatives_5))

          CI_Lo, CI_Hi = calcBin(float(sum_ME_coverage),  SUM_PSI)

      else:

          PSI = "NA"
          CI_Lo, CI_Hi = ["NA", "NA"]
          
      print(row["FILE_NAME"], row["ME"], PSI, CI_Lo, CI_Hi)
  
  
  
  
if __name__ == '__main__':
	main(sys.argv[1] )

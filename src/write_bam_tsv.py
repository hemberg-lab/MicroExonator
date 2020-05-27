from snakemake.utils import min_version
import csv


with snakemake.input[0] as f:
  
  reader = csv.DictReader(f, delimiter="\t")
  
  for row in reader:
  

  
  
  

from snakemake.utils import min_version
import glob, os, csv


for path in glob.glob(snakemake.params[0] + "/*.*"):
  
  compare_name = path.split("/")[-1].split(".")[0]
  
  with open(path) as file  :
    
    reader = csv.DictReader(file, delimiter="\t")
    
    for row in reader:
      
      with open("Whippet/ggsashimi/" + compare_name + "/" + "_".join([row["Gene"], row["Node"], row["Strand"], ".txt"] ), "w") as output:
        
        output.write(" ".join([row["Gene"], row["Node"], row["Strand"]] ))
      
        
  
  

## Integreate Whippet and Snakemake ##

First, you need to create a environment that has `snakemake` and the version of `julia` that is compatible with `Whipet v0.11`. For this, we can create a conda environment from `Whippet/julia_0.6.1.yaml`, which is yaml file with the recipe to create a stable environment with julia 0.6.1 and snakemake. 

`conda env create -f Whippet/julia_0.6.1.yaml`

Then, activate the newly created enviroment:

`source activate julia_0.6.1`

Enter julia's interactive mode:

`julia`

Install Whippet:

`Pkg.add("Whippet")`

Exit interactive julia session (`control + d`) and find Whippet's binary folder, that should be inside of your miniconda environment folder. Once you find the path to this folder, add it to `config.yaml` writing the following lines:

    whippet_bin_folder : path/to/miniconda/envs/julia/share/julia/site/v0.6/Whippet/bin/whippet-quant.jl
    Gene_anontation_GTF : path/to/gene_annotation.gtf
    condition1 : sample1,sample2,sample3... 
    condition2 : sample4,sample5,sample6... 
    comparison_name :  cond1_vs_cond2

Where:
* condition1 and condition2 are coma-separated list of the sample names that you want to compare.
* comparison_name is the name of the comparison.


comparison_names = whippet_delta.keys()


if "whippet_delta" in config:

      if str2bool(config.get("Only_whippet", False)):
            rule differential_inclusion:
                input:
                    expand("Whippet/Delta/{comparison_name}.diff.gz", comparison_name=comparison_names)
      else:
            rule differential_inclusion:
                input:
                    expand("Whippet/Delta/{comparison_name}.diff.microexons", comparison_name=comparison_names),
                    expand("Whippet/Delta/{comparison_name}.diff.ME.microexons", comparison_name=comparison_names)


rule whippet_delta:
    input:
        lambda wildcards : expand("Whippet/Quant/{sample}.psi.gz", sample= whippet_delta[wildcards.comparison_name]["A"].split(",")),
        lambda wildcards : expand("Whippet/Quant/{sample}.psi.gz", sample= whippet_delta[wildcards.comparison_name]["B"].split(","))
    output:
        "Whippet/Delta/{comparison_name}.diff.gz"
    params:
        bin = config["whippet_bin_folder"],
        a = lambda wildcards : ",".join(expand("Whippet/Quant/{sample}.psi.gz", sample= whippet_delta[wildcards.comparison_name]["A"].split(","))),
        b = lambda wildcards : ",".join(expand("Whippet/Quant/{sample}.psi.gz", sample= whippet_delta[wildcards.comparison_name]["B"].split(","))),
        o = lambda wildcards : "Whippet/Delta/" + wildcards.comparison_name,
        julia = config["julia"]
    shell:
        "{params.julia} {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o}"



rule whippet_delta_ME:
    input:
        lambda wildcards : expand("Whippet/Quant/{sample}.psi.ME.gz", sample= whippet_delta[wildcards.comparison_name]["A"].split(",")),
        lambda wildcards : expand("Whippet/Quant/{sample}.psi.ME.gz", sample= whippet_delta[wildcards.comparison_name]["B"].split(","))
    output:
        "Whippet/Delta/{comparison_name}.ME.diff.gz"
    params:
        bin = config["whippet_bin_folder"],
        a = lambda wildcards : ",".join(expand("Whippet/Quant/{sample}.psi.ME.gz", sample= whippet_delta[wildcards.comparison_name]["A"].split(","))),
        b = lambda wildcards : ",".join(expand("Whippet/Quant/{sample}.psi.ME.gz", sample= whippet_delta[wildcards.comparison_name]["B"].split(","))),
        o = lambda wildcards : "Whippet/Delta/" + wildcards.comparison_name + ".ME",
        julia = config["julia"]
    shell:
        "{params.julia} {params.bin}/whippet-delta.jl -a {params.a} -b {params.b} -o {params.o} "



## Differential inclusion without Whippet (delta_method : "microexonator").
## src/me_delta.py applies the whippet-delta.jl model to MicroExonator's own
## per-sample PSI and read counts. Comparisons come from the same YAML file as
## the Whippet route (config key whippet_delta).

comparison_names = whippet_delta.keys()


def delta_group_tables(comparison_name, group):
    return [downstream_PSI(sample) for sample in whippet_delta[comparison_name][group].split(",")]


rule differential_inclusion:
    input:
        expand("Delta/{comparison_name}.diff.ME.microexons", comparison_name=comparison_names)


rule microexonator_delta:
    input:
        A = lambda wildcards : delta_group_tables(wildcards.comparison_name, "A"),
        B = lambda wildcards : delta_group_tables(wildcards.comparison_name, "B"),
        microexons = FILTERED_ME_OUTPUT
    output:
        "Delta/{comparison_name}.diff.ME.microexons"
    params:
        a = lambda wildcards, input : ",".join(input.A),
        b = lambda wildcards, input : ",".join(input.B),
        gtf = "--gtf " + config["Gene_anontation_GTF"] if "Gene_anontation_GTF" in config else "",
        min_reads = config.get("delta_min_reads", 5),
        min_samples = config.get("delta_min_samples", 1),
        size = config.get("delta_empirical_size", 1000),
        seed = config.get("delta_seed", 123456)
    conda:
        "../envs/core.yaml"
    shell:
        "python3 src/me_delta.py -a {params.a} -b {params.b} --microexons {input.microexons} {params.gtf} "
        "--min-reads {params.min_reads} --min-samples {params.min_samples} --size {params.size} --seed {params.seed} > {output}"

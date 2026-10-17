"""Optional PSI-Sigma module on the shared HISAT2 alignment (src/run_psisigma.py).

PSI-Sigma v2.3 and the Perl modules its statistics need are installed once
(umbrella_psisigma_software, login node: it downloads from GitHub and CPAN).
Intron retention uses PSI-Sigma's default irmode 0, which depends only on the
annotation, so each run keeps its junction counts and its coverage at the
retained-intron positions; comparisons read those caches, never a BAM.
"""

if str(config.get("psisigma_irmode", 0)) != "0":
    raise WorkflowError("psisigma_irmode must be 0: other PSI-Sigma intron-retention modes depend on each "
                        "comparison's junctions and cannot reuse per-run caches")

PSISIGMA_SOFTWARE = "umbrella/modules/software/psisigma/bfcb1423f1ec61ac48b6f92250505ec5855f4a53"
PSISIGMA_INSTALLED = PSISIGMA_SOFTWARE + "/installed.json"
PSISIGMA_PERL5LIB = PSISIGMA_SOFTWARE + "/perl5/lib/perl5"
PSISIGMA_REFERENCE = {name: MODULE_REFERENCE["psisigma"] + "/" + name
                      for name in ("annotation.psisigma.gtf", "annotation.psisigma.gtf.mapping.txt",
                                   "ir_universe.tsv.gz", "genome.fa.fai", "reference.json")}

localrules: umbrella_psisigma_software


rule umbrella_psisigma_software:
    output:
        protected(PSISIGMA_INSTALLED)
    params:
        root=PSISIGMA_SOFTWARE
    log:
        PSISIGMA_SOFTWARE + ".install.log"
    conda:
        "../envs/umbrella-psisigma.yaml"
    shell:
        "python3 src/install_psisigma.py --root {params.root} --record {output} > {log} 2>&1"


rule umbrella_psisigma_reference:
    input:
        gtf=UMBRELLA_GTF["annotation"],
        genome=UMBRELLA_REFERENCE["genome_fasta"],
        software=PSISIGMA_INSTALLED
    output:
        [protected(path) for path in PSISIGMA_REFERENCE.values()]
    params:
        out=MODULE_REFERENCE["psisigma"],
        software=PSISIGMA_SOFTWARE
    threads: 4
    log:
        MODULE_REFERENCE["psisigma"] + "/reference.stdout.log"
    conda:
        "../envs/umbrella-psisigma.yaml"
    shell:
        "python3 src/run_psisigma.py reference --gtf {input.gtf} --genome {input.genome} "
        "--psisigma {params.software} --out-dir {params.out} --threads {threads} > {log} 2>&1"


if MODULE_PRODUCERS:
    rule umbrella_psisigma_prep:
        input:
            bam=module_ingest("psisigma", lambda w: umbrella_quant_run(w).work_dir + "/align/aligned.bam"),
            bai=module_ingest("psisigma", lambda w: umbrella_quant_run(w).work_dir + "/align/aligned.bam.bai"),
            universe=module_ingest("psisigma", PSISIGMA_REFERENCE["ir_universe.tsv.gz"]),
            fai=module_ingest("psisigma", PSISIGMA_REFERENCE["genome.fa.fai"]),
            genome=module_ingest("psisigma", UMBRELLA_REFERENCE["genome_fasta"]),
            software=module_ingest("psisigma", PSISIGMA_INSTALLED)
        output:
            [protected(MODULE_CACHE.replace("{tool}", "psisigma") + "/" + name) for name in CACHE_FILES["psisigma"]]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace("{tool}", "psisigma"),
            software=PSISIGMA_SOFTWARE,
            layout=lambda w: umbrella_quant_run(w).layout,
            junctions="" if MODULE_SETTINGS["psisigma"]["canonical_only"] else "--all-junctions",
            identity=lambda w: module_identity(w, "psisigma")
        threads: 3
        log:
            "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.psisigma.{cache_id}.log"
        benchmark:
            "umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.psisigma.{cache_id}.tsv"
        conda:
            "../envs/umbrella-psisigma.yaml"
        shell:
            "python3 src/run_psisigma.py prep --bam {input.bam} --universe {input.universe} "
            "--psisigma {params.software} --genome {input.genome} --fai {input.fai} --out-dir {params.out} "
            "--run-id {wildcards.run_id} --layout {params.layout} --threads {threads} {params.junctions} "
            "--identity {params.identity:q} > {log} 2>&1"


rule umbrella_psisigma_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        caches=module_comparison_caches("psisigma"),
        reference=list(PSISIGMA_REFERENCE.values()),
        software=PSISIGMA_INSTALLED
    output:
        results=MODULE_ANALYSIS.replace("{tool}", "psisigma") + "/results.tsv.gz",
        native=MODULE_ANALYSIS.replace("{tool}", "psisigma") + "/native.tar.gz",
        status=MODULE_ANALYSIS.replace("{tool}", "psisigma") + "/status.json"
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace("{tool}", "psisigma"),
        caches=lambda w: module_cache_arguments(w, "psisigma"),
        reference=MODULE_REFERENCE["psisigma"],
        software=PSISIGMA_SOFTWARE,
        perl5lib=PSISIGMA_PERL5LIB,
        settings=_json.dumps(MODULE_SETTINGS["psisigma"], sort_keys=True),
        ids=lambda w: module_ids(w, "psisigma", SOFTWARE["psisigma"])
    threads: 4
    benchmark:
        MODULE_ANALYSIS.replace("{tool}", "psisigma") + "/benchmark.tsv"
    conda:
        "../envs/umbrella-psisigma.yaml"
    shell:
        "python3 src/run_psisigma.py compare --preflight {input.preflight} {params.caches} "
        "--reference {params.reference} --psisigma {params.software} --perl5lib {params.perl5lib} "
        "--out-dir {params.out} --threads {threads} --settings {params.settings:q} --ids {params.ids:q}"


rule prepare_psisigma:
    input: module_target("psisigma", module_prepare_targets)


rule quant_psisigma:
    input: module_target("psisigma", module_quant_targets)

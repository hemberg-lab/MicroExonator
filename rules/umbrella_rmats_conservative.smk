"""Optional independent rMATS profile: no novel splice-site discovery."""

if MODULE_PRODUCERS:
    rule umbrella_rmats_conservative_prep:
        input:
            bam=module_ingest('rmats_conservative', lambda w: umbrella_quant_run(w).work_dir + '/align/aligned.bam'),
            bai=module_ingest('rmats_conservative', lambda w: umbrella_quant_run(w).work_dir + '/align/aligned.bam.bai'),
            length=module_ingest('rmats_conservative', lambda w: umbrella_quant_run(w).work_dir + '/align/read_length.txt'),
            gtf=module_ingest('rmats_conservative', UMBRELLA_RMATS_GTF)
        output:
            [protected(MODULE_CACHE.replace('{tool}', 'rmats_conservative') + '/' + name)
             for name in CACHE_FILES['rmats_conservative']]
        wildcard_constraints:
            **MODULE_WILDCARDS
        params:
            out=MODULE_CACHE.replace('{tool}', 'rmats_conservative'),
            layout=lambda w: umbrella_quant_run(w).layout,
            ids=lambda w: module_identity(w, 'rmats_conservative')
        threads: 4
        conda:
            '../envs/umbrella-alignment.yaml'
        benchmark:
            'umbrella/modules/benchmarks/{reference_id}/{project_id}/{batch_id}/{run_id}.rmats_conservative.{cache_id}.tsv'
        shell:
            'python3 src/run_rmats_profile.py prep --bam {input.bam:q} --gtf {input.gtf:q} '
            '--read-length $(cat {input.length:q}) --layout {params.layout} --run-id {wildcards.run_id:q} '
            '--out-dir {params.out:q} --threads {threads} --ids {params.ids:q}'


rule umbrella_rmats_conservative_compare:
    input:
        preflight=UMBRELLA_COMPARISON_ROOT + '/preflight.json',
        caches=module_comparison_caches('rmats_conservative'),
        prep=module_comparison_caches('rmats_conservative', 'run.rmats'),
        gtf=UMBRELLA_RMATS_GTF
    output:
        results=MODULE_ANALYSIS.replace('{tool}', 'rmats_conservative') + '/results.tsv.gz',
        native=MODULE_ANALYSIS.replace('{tool}', 'rmats_conservative') + '/native.tar.gz',
        status=MODULE_ANALYSIS.replace('{tool}', 'rmats_conservative') + '/status.json'
    wildcard_constraints:
        **MODULE_WILDCARDS
    params:
        out=MODULE_ANALYSIS.replace('{tool}', 'rmats_conservative'),
        caches=lambda w: module_cache_arguments(w, 'rmats_conservative'),
        ids=lambda w: module_ids(w, 'rmats_conservative', SOFTWARE['rmats_conservative'])
    threads: 4
    conda:
        '../envs/umbrella-alignment.yaml'
    benchmark:
        MODULE_ANALYSIS.replace('{tool}', 'rmats_conservative') + '/benchmark.tsv'
    shell:
        'python3 src/run_rmats_profile.py compare --preflight {input.preflight:q} {params.caches} '
        '--gtf {input.gtf:q} --out-dir {params.out:q} --threads {threads} --ids {params.ids:q}'


rule prepare_rmats_conservative:
    input: module_target('rmats_conservative', module_prepare_targets)


rule quant_rmats_conservative:
    input: module_target('rmats_conservative', module_quant_targets)

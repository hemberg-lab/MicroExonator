"""Standalone cached-only PCA: intentionally no FASTQ or matrix producers."""
from src.expression_pca_cache import cached_comparisons
from shlex import quote

PCA_ROOT = config.get('umbrella_pca_root', 'comparisons')
PCA_ROUTES = config.get('umbrella_pca_routes', ['featurecounts'])
PCA_JOBS = cached_comparisons(PCA_ROOT, PCA_ROUTES, config.get('umbrella_pca_comparisons'))

rule expression_pca:
    input:
        [PCA_ROOT + '/' + key.rsplit('/', 1)[0] + '/pca_' + key.rsplit('/', 1)[1] + extension
         for key in PCA_JOBS for extension in ('.pdf', '.scores.tsv.gz')]

rule cached_expression_pca:
    input:
        script='src/expression_pca.R',
        preflight=lambda w: PCA_JOBS['/'.join((w.reference_id, w.project_id, w.comparison_id, w.route))]['preflight'],
        matrices=lambda w: PCA_JOBS['/'.join((w.reference_id, w.project_id, w.comparison_id, w.route))]['matrices']
    output:
        pdf=PCA_ROOT + '/{reference_id}/{project_id}/{comparison_id}/pca_{route}.pdf',
        scores=PCA_ROOT + '/{reference_id}/{project_id}/{comparison_id}/pca_{route}.scores.tsv.gz'
    wildcard_constraints:
        route='featurecounts|tximport'
    params:
        prefix=PCA_ROOT + '/{reference_id}/{project_id}/{comparison_id}/pca_{route}',
        min_group=config.get('umbrella_pca_min_group', 5),
        inputs=lambda w, input: ('--counts ' + quote(input.matrices[0]) if w.route == 'featurecounts'
            else '--counts {} --tpm {} --length {}'.format(*[quote(p) for p in input.matrices]))
    threads: 1
    resources:
        mem_mb=4000
    log:
        PCA_ROOT + '/{reference_id}/{project_id}/{comparison_id}/logs/pca_{route}.log'
    conda:
        'envs/umbrella-pca.yaml'
    shell:
        'Rscript {input.script:q} --route {wildcards.route:q} --preflight {input.preflight:q} '
        '{params.inputs} --out {params.prefix:q} --min-group {params.min_group:q} > {log:q} 2>&1'

"""Optional cached-only gene expression QC; no read or count producers."""
from src.expression_pca_cache import cached_comparisons, qc_routes, retained_qc_inputs
from shlex import quote

QC_ROOT = config.get('umbrella_expression_qc_root', 'comparisons')
QC_ROUTES = qc_routes(config.get('umbrella_expression_qc_sources', []))
QC_JOBS = cached_comparisons(QC_ROOT, QC_ROUTES, config.get('umbrella_expression_qc_comparisons')) if QC_ROUTES else {}
for job in QC_JOBS.values():
    job['qc'] = retained_qc_inputs(QC_ROOT, job['preflight'])

def qc_job(w):
    return QC_JOBS['/'.join((w.reference_id, w.project_id, w.comparison_id, w.route))]

rule gene_expression_qc:
    input:
        [QC_ROOT + '/' + key.rsplit('/', 1)[0] + '/gene_expression_qc_' + key.rsplit('/', 1)[1] + ext
         for key in QC_JOBS for ext in ('.pdf', '.scores.tsv.gz', '.summary.tsv.gz')]

rule cached_gene_expression_qc:
    input:
        scripts=['src/expression_pca.R', 'src/gene_expression_qc.R'],
        preflight=lambda w: qc_job(w)['preflight'],
        matrices=lambda w: qc_job(w)['matrices'],
        qc=lambda w: qc_job(w)['qc']
    output:
        pdf=QC_ROOT + '/{reference_id}/{project_id}/{comparison_id}/gene_expression_qc_{route}.pdf',
        scores=QC_ROOT + '/{reference_id}/{project_id}/{comparison_id}/gene_expression_qc_{route}.scores.tsv.gz',
        summary=QC_ROOT + '/{reference_id}/{project_id}/{comparison_id}/gene_expression_qc_{route}.summary.tsv.gz'
    wildcard_constraints:
        route='featurecounts|tximport'
    params:
        prefix=QC_ROOT + '/{reference_id}/{project_id}/{comparison_id}/gene_expression_qc_{route}',
        cooks=str(config.get('umbrella_expression_qc_cooks', True)).lower(),
        inputs=lambda w, input: ('--counts ' + quote(input.matrices[0]) if w.route == 'featurecounts'
            else '--counts {} --tpm {} --length {}'.format(*[quote(p) for p in input.matrices])),
        qc=lambda w, input: quote(__import__('json').dumps(list(input.qc)))
    threads: 1
    resources:
        mem_mb=4000
    log:
        QC_ROOT + '/{reference_id}/{project_id}/{comparison_id}/logs/gene_expression_qc_{route}.log'
    conda:
        'envs/umbrella-pca.yaml'
    shell:
        'Rscript src/expression_pca.R --route {wildcards.route:q} --preflight {input.preflight:q} '
        '{params.inputs} --out {params.prefix:q} --qc-report true --cooks {params.cooks:q} '
        '--qc-files {params.qc} > {log:q} 2>&1'

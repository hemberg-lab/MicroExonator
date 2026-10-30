"""Expression PCA per comparison inside quant_umbrella (umbrella_pca_routes).

Same script, outputs and conda env as the standalone cached ExpressionPCA.smk, but
fed by the comparison joins of this run, so one call quantifies and plots. Every
gene-expression route the run quantifies gets its PCA: tximport when salmon is
on, featurecounts when hisat2 is on. `umbrella_pca_routes` overrides that:

  umbrella_pca_routes: [tximport]        # only the Salmon PCA
  umbrella_pca_routes: []                # no PCA
  umbrella_pca_min_group: 5              # within-group distance flags; at least 5
"""
import shlex

from src.expression_pca_cache import ROUTES as UMBRELLA_PCA_KINDS

UMBRELLA_PCA_NEEDS = {"featurecounts": "hisat2", "tximport": "salmon"}
if config.get("umbrella_pca_routes") is None:
    UMBRELLA_PCA_ROUTES = [route for route, tool in sorted(UMBRELLA_PCA_NEEDS.items())
                           if UMBRELLA_OPTIONAL[tool]]
else:
    UMBRELLA_PCA_ROUTES = list(config["umbrella_pca_routes"])
_pca_unknown = sorted(set(UMBRELLA_PCA_ROUTES) - set(UMBRELLA_PCA_NEEDS))
if _pca_unknown or len(UMBRELLA_PCA_ROUTES) != len(set(UMBRELLA_PCA_ROUTES)):
    raise WorkflowError("umbrella_pca_routes: unknown or repeated route {}; choose from {}".format(
        ", ".join(_pca_unknown or UMBRELLA_PCA_ROUTES), ", ".join(UMBRELLA_PCA_NEEDS)))
_pca_off = ["{} needs {}".format(route, UMBRELLA_PCA_NEEDS[route]) for route in UMBRELLA_PCA_ROUTES
            if not UMBRELLA_OPTIONAL[UMBRELLA_PCA_NEEDS[route]]]
if _pca_off:
    raise WorkflowError("umbrella_pca_routes: {}; switch the tool on in umbrella_optional".format(
        "; ".join(_pca_off)))


rule umbrella_expression_pca:
    input:
        script="src/expression_pca.R",
        preflight=UMBRELLA_COMPARISON_ROOT + "/preflight.json",
        matrices=lambda w: [UMBRELLA_COMPARISON_ROOT.format(**w) + "/joined/{}.tsv.gz".format(kind)
                            for kind in UMBRELLA_PCA_KINDS[w.route]]
    output:
        pdf=UMBRELLA_COMPARISON_ROOT + "/pca_{route}.pdf",
        scores=UMBRELLA_COMPARISON_ROOT + "/pca_{route}.scores.tsv.gz"
    wildcard_constraints:
        route="featurecounts|tximport"
    params:
        prefix=UMBRELLA_COMPARISON_ROOT + "/pca_{route}",
        min_group=config.get("umbrella_pca_min_group", 5),
        inputs=lambda w, input: ("--counts " + shlex.quote(input.matrices[0]) if w.route == "featurecounts"
            else "--counts {} --tpm {} --length {}".format(*[shlex.quote(p) for p in input.matrices]))
    threads: 1
    log:
        UMBRELLA_COMPARISON_ROOT + "/logs/pca_{route}.log"
    conda:
        "../envs/umbrella-pca.yaml"
    shell:
        "Rscript {input.script:q} --route {wildcards.route:q} --preflight {input.preflight:q} "
        "{params.inputs} --out {params.prefix:q} --min-group {params.min_group:q} > {log:q} 2>&1"


UMBRELLA_PCA_TARGETS = [
    "comparisons/{}/{}/{}/pca_{}{}".format(result["reference_id"], result["project_id"], comparison_id,
                                          route, extension)
    for comparison_id, result in sorted(UMBRELLA_COMPARISONS.items())
    for route in UMBRELLA_PCA_ROUTES for extension in (".pdf", ".scores.tsv.gz")]

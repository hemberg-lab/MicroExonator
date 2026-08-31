"""Snakemake preflight for the historical Gaussian-mixture filter."""

from src.legacy_mixture import (
    load_legacy_mixture_scores,
    validate_legacy_mixture_scores,
)


scores = load_legacy_mixture_scores(
    snakemake.input["ME_coverage"],
    snakemake.input["ME_matches_file"],
    int(snakemake.params["min_number_files_detected"]),
)
validate_legacy_mixture_scores(scores)

with open(snakemake.output["validation"], "w") as output:
    output.write("observations\t{}\n".format(len(scores)))
    output.write("unique_scores\t{}\n".format(len(set(scores))))

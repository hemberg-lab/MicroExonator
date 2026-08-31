"""Snakemake adapter for the final robustness-filtered ME list."""

from src.robustness_filter import (
    collect_robust_microexons,
    write_robust_microexons,
)


detection_files = (
    list(snakemake.input["bulk_se"])
    + list(snakemake.input["bulk_pe"])
    + list(snakemake.input["single_cell"])
)
robust_microexons = collect_robust_microexons(
    detection_files,
    min_detected_samples=int(snakemake.params["min_detected_samples"]),
)
write_robust_microexons(
    snakemake.input["ME_centric"],
    snakemake.output["detected_list"],
    robust_microexons,
)
